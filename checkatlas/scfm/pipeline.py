"""Top-level scFM QC pipeline orchestration.

This module ties together the three layers:

  1. ``metrics.scfm.run.cal_scfm``   — collects raw metric values
  2. ``scfm.diagnostics.diagnose``    — produces nine verdicts
  3. ``scfm.composite.compute_all``   — produces FMF, BF, PR
  4. ``scfm.report.write_all``       — writes MultiQC-friendly TSVs

It is intentionally a *thin* orchestrator. All real logic lives in
the sub-modules. The function ``run_scfm_pipeline`` is the public
entry point called from ``__main__.py``.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .. import check
from ..metrics import cluster, dimred
from ..metrics import annot as annot_mod
from ..metrics.scfm import run as scfm_run
from ..utils import folders
from . import diagnostics, report
from .config import SCFMConfig

logger = logging.getLogger("checkatlas")


def _auto_pick_scfm_embedding(adata: AnnData, hint: str) -> str:
    """Pick the scFM embedding to evaluate.

    If the user supplied ``--scfm_embedding`` we trust them. Otherwise
    we use the column detector (Becavin comment 1) to auto-detect
    any embedding whose name matches an scFM-specific pattern.
    """
    if hint and hint in adata.obsm:
        return hint
    if hint:
        logger.debug("scfm_embedding %r not in adata.obsm; falling back to auto-detect", hint)
    from ..utils.col_detector import CheckAtlasColumnDetector

    detector = CheckAtlasColumnDetector(adata)
    matched = detector.analyze_obsm_semantics()
    if not matched:
        return hint or ""
    scfm_patterns = (
        "geneformer",
        "scgpt",
        "uce",
        "scfoundation",
        "sccello",
        "scbert",
        "scfm",
        "foundation",
    )
    for key, _score in matched:
        key_l = key.lower()
        if any(p in key_l for p in scfm_patterns):
            return key
    return matched[0][0]


def _load_or_build_context(adata: AnnData, atlas_name: str, outdir: str, args) -> Optional[Any]:
    """Build (or load) a PreprocessContext for the scfm pipeline.

    This is a *thin* wrapper around ``atlas.preprocess_atlas``. It
    re-uses the kNN / distance-matrix cache that the existing
    pipeline builds, satisfying Becavin comment 2 (no re-processing
    of adata).
    """
    from ..atlas import preprocess_atlas

    try:
        adata_proc = preprocess_atlas(
            {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_PATH_KEY: getattr(adata, "filename", None) or "",
                check.ATLAS_TYPE_KEY: "AnnData",
                check.ATLAS_EXTENSION_KEY: ".h5ad",
            },
            args,
        )
    except Exception as exc:
        logger.debug("scfm: preprocess_atlas failed: %s", exc)
        return None
    return adata_proc, _try_load_context(atlas_name, outdir, adata_proc)


def _try_load_context(atlas_name: str, outdir: str, adata: AnnData) -> Optional[Any]:
    """Try to load a PreprocessContext from the temp folder."""
    from ..metrics._preprocess_context import load_context, make_preprocess_fingerprint

    try:
        temp_parent = folders.get_folder(outdir, folders.TEMP)
        fp = make_preprocess_fingerprint(
            adata,
            embedding_keys=list(adata.obsm.keys()),
            cluster_label_keys=[],
            batch_keys=[],
            k_neighbors=90,
            source_path=getattr(adata, "filename", None),
        )
        return load_context(atlas_name, temp_parent, fp)
    except Exception as exc:
        logger.debug("scfm: load_context failed: %s", exc)
        return None


def run_scfm_pipeline(
    adata: AnnData,
    config: SCFMConfig,
    args=None,
    *,
    outdir: Optional[str] = None,
) -> dict[str, Any]:
    """Run the full scfm QC pipeline on one AnnData.

    Parameters
    ----------
    adata : AnnData
    config : SCFMConfig
    args : argparse.Namespace, optional
        Used for ``args.path`` and ``args.n_jobs``.
    outdir : str, optional
        Output directory. Defaults to
        ``checkatlas_files/scfm/<atlas_name>/`` under the current
        path.

    Returns
    -------
    dict
        ``{"verdicts": [...], "composite": {...}, "out_paths": {...}}``.
    """
    if outdir is None:
        base = getattr(args, "path", ".") if args is not None else "."
        outdir = os.path.join(
            folders.get_folder(base, folders.SCFM), config.atlas_name
        )
    folders.checkatlas_folders(os.path.dirname(outdir) or ".")
    os.makedirs(outdir, exist_ok=True)

    n_jobs = getattr(args, "n_jobs", -1) if args is not None else -1

    scfm_emb = _auto_pick_scfm_embedding(adata, config.scfm_embedding)
    config_eff = SCFMConfig(
        **{**config.__dict__, "scfm_embedding": scfm_emb}
    )
    if not scfm_emb:
        logger.warning(
            "scfm: no scFM embedding detected in adata.obsm; "
            "diagnostics will be skipped (n/a)"
        )

    rows: list[dict] = []

    existing_metric_lists: list = []
    for metric_name in (
        "isolated_f1_score",
        "adj_rand_index",
        "average_silhouette_width",
        "kbet",
        "pcr",
        "graph_connectivity",
    ):
        if metric_name in annot_mod.__all__:
            existing_metric_lists.append(metric_name)
    for metric_name in ("silhouette", "davies_bouldin", "calinski_harabasz"):
        if metric_name in cluster.__all__:
            existing_metric_lists.append(metric_name)
    for metric_name in ("dCor", "kruskal_stress", "spearman_rho"):
        if metric_name in dimred.__all__:
            existing_metric_lists.append(metric_name)

    _, ctx = (
        _load_or_build_context(adata, config.atlas_name, outdir, args)
        if (args is not None and getattr(args, "path", None))
        else (adata, None)
    )

    try:
        annot_df = annot_mod.cal_annot(
            adata,
            atlas_name=config.atlas_name,
            metric_list=existing_metric_lists,
            file_dir=os.path.join(outdir, "annot_cache"),
            n_jobs=n_jobs,
            verbose=False,
            preprocess_context=ctx,
        )
        for _, r in annot_df.iterrows():
            rows.append(
                {
                    "Atlas Name": config.atlas_name,
                    "Task": "annot",
                    "Metric Name": str(r.get("Metric Name", "")),
                    "Reference/Input 1": str(r.get("Reference/Input 1", "")),
                    "Prediction/Input 2": str(r.get("Prediction/Input 2", "")),
                    "Embedding": str(r.get("Reference/Input 1", "")),
                    "Value": r.get("Value", np.nan),
                    "Time (s)": r.get("Time (s)", 0.0),
                }
            )
    except Exception as exc:
        logger.debug("scfm: cal_annot failed: %s", exc)

    try:
        clust_df = metrics_mod_cal_cluster(
            adata,
            atlas_name=config.atlas_name,
            metric_list=["silhouette"],
            file_dir=os.path.join(outdir, "cluster_cache"),
            n_jobs=n_jobs,
            verbose=False,
            preprocess_context=ctx,
        )
        for _, r in clust_df.iterrows():
            rows.append(
                {
                    "Atlas Name": config.atlas_name,
                    "Task": "cluster",
                    "Metric Name": str(r.get("Metric Name", "")),
                    "Embedding": str(r.get("Embedding", "")),
                    "Label Key": str(r.get("Label Key", "")),
                    "Value": r.get("Value", np.nan),
                    "Time (s)": r.get("Time (s)", 0.0),
                }
            )
    except Exception as exc:
        logger.debug("scfm: cal_cluster failed: %s", exc)

    try:
        scfm_df = scfm_run.cal_scfm(
            adata,
            scfm_embedding=scfm_emb,
            baseline_embeddings=config.baseline_embeddings,
            ref_label=config.ref_label,
            pred_label=config.predicted_label,
            batch_key=config.batch_key,
            domain_key=config.domain_key,
            patient_key=config.patient_key,
            outcome_key=config.outcome_key,
            atlas_name=config.atlas_name,
            scaling_fractions=config.scaling_fractions,
            n_seeds=config.n_seeds,
            noise_sigma=config.noise_sigma,
            min_domain_size=config.min_domain_size,
            n_jobs=n_jobs,
            preprocess_context=ctx,
        )
        for _, r in scfm_df.iterrows():
            row = {
                "Atlas Name": config.atlas_name,
                "Task": str(r.get("Task", "")),
                "Metric Name": str(r.get("Metric Name", "")),
                "Embedding": str(r.get("Embedding", "")),
                "Value": r.get("Value", np.nan),
                "Time (s)": r.get("Time (s)", 0.0),
            }
            for extra_col in ("Fraction", "N_Cells", "Train_Domain", "Test_Domain"):
                if extra_col in r.index:
                    row[extra_col] = r[extra_col]
            rows.append(row)
    except Exception as exc:
        logger.warning("scfm: cal_scfm failed: %s", exc)
        scfm_df = pd.DataFrame()

    metrics_df = pd.DataFrame(rows) if rows else pd.DataFrame()
    verdicts = diagnostics.diagnose(metrics_df, config_eff)

    config_dict = {
        "atlas": config.atlas_name,
        "scfm_embedding": scfm_emb,
        "baseline_embeddings": ",".join(config.baseline_embeddings),
        "predicted_label": config.predicted_label or "",
        "ref_label": config.ref_label or "",
        "batch_key": config.batch_key or "",
        "domain_key": config.domain_key or "",
        "patient_key": config.patient_key or "",
        "outcome_key": config.outcome_key or "",
        "scaling_fractions": ",".join(str(f) for f in config.scaling_fractions),
        "n_seeds": config.n_seeds,
        "n_cells": adata.n_obs,
        "n_embeddings": len(adata.obsm),
    }
    out_paths = report.write_all(
        atlas_name=config.atlas_name,
        verdicts=verdicts,
        metrics_df=metrics_df,
        outdir=outdir,
        weights_path=config.weights_path,
        thresholds_path=config.thresholds_path,
        config_dict=config_dict,
    )
    return {
        "verdicts": verdicts,
        "metrics_df": metrics_df,
        "out_paths": out_paths,
        "outdir": outdir,
    }


def metrics_mod_cal_cluster(*args, **kwargs):
    """Thin re-export of ``metrics.cal_cluster`` (avoids a top-level import cycle)."""
    from ..metrics import metrics

    return metrics.cal_cluster(*args, **kwargs)
