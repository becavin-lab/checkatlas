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


def _load_or_build_context(adata: AnnData, atlas_name: str, outdir: str, args) -> tuple[AnnData, Optional[Any]]:
    """Load (or build) a PreprocessContext for the scfm pipeline.

    This is a *thin* wrapper around ``atlas.preprocess_atlas``. It
    re-uses the kNN / distance-matrix cache that the existing
    pipeline builds, satisfying Becavin comment 2 (no re-processing
    of adata).

    The function ALWAYS attempts to load or build a context — it is
    no longer gated on ``args.path`` being set. The previous behaviour
    silently returned ``ctx=None`` for notebook and programmatic
    callers, which forced the downstream scfm metric modules to fall
    back to building their own kNN graphs. This broke the GPU path:
    JAX-accelerated kNN lives in ``metrics._neighbors`` and is only
    triggered through the preprocess pipeline.
    """
    from ..atlas import preprocess_atlas
    from ..metrics._preprocess_context import (
        load_context,
        make_preprocess_fingerprint,
    )

    # Step 1: try to load a cached context from the temp folder.
    temp_parent = folders.get_folder(outdir, folders.TEMP)
    try:
        fp = make_preprocess_fingerprint(
            adata,
            embedding_keys=list(adata.obsm.keys()),
            cluster_label_keys=[],
            batch_keys=[],
            k_neighbors=90,
            source_path=getattr(adata, "filename", None),
        )
        cached = load_context(atlas_name, temp_parent, fp)
        if cached is not None:
            logger.info(
                "scfm: PreprocessContext cache hit for %s — "
                "re-using precomputed kNN / distance matrices",
                atlas_name,
            )
            return adata, cached
    except Exception as exc:
        logger.debug("scfm: load_context failed: %s", exc)

    # Step 2: build a fresh context. ``preprocess_atlas`` reads adata
    # from disk via ``atlas_info`` so we provide the h5ad path. When
    # called from a notebook, ``adata.filename`` is set and we use it;
    # otherwise we fall back to the in-memory AnnData (the path-only
    # build will fail and we return ``(adata, None)`` so the metric
    # modules can still run, just without the precomputed cache).
    source_path = getattr(adata, "filename", None)
    if not source_path or not os.path.exists(source_path):
        logger.debug(
            "scfm: no on-disk source for %s; running without "
            "PreprocessContext (CPU-only downstream metrics)",
            atlas_name,
        )
        return adata, None

    try:
        adata_proc = preprocess_atlas(
            {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_PATH_KEY: source_path,
                check.ATLAS_TYPE_KEY: "AnnData",
                check.ATLAS_EXTENSION_KEY: ".h5ad",
            },
            args,
        )
    except Exception as exc:
        logger.debug("scfm: preprocess_atlas failed: %s", exc)
        return adata, None

    # Step 3: load the freshly-built context.
    try:
        fp = make_preprocess_fingerprint(
            adata_proc,
            embedding_keys=list(adata_proc.obsm.keys()),
            cluster_label_keys=[],
            batch_keys=[],
            k_neighbors=90,
            source_path=source_path,
        )
        ctx = load_context(atlas_name, temp_parent, fp)
    except Exception as exc:
        logger.debug("scfm: post-build load_context failed: %s", exc)
        ctx = None

    return adata_proc, ctx


def _try_load_context(atlas_name: str, outdir: str, adata: AnnData) -> Optional[Any]:
    """Try to load a PreprocessContext from the temp folder.

    Kept for backwards compatibility — ``_load_or_build_context``
    supersedes this. New code should call ``_load_or_build_context``.
    """
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

    # ALWAYS try to load or build a PreprocessContext. The previous
    # guard ``if (args is not None and getattr(args, "path", None))``
    # silently dropped the context for notebook / programmatic
    # callers, forcing the scfm metric modules to build their own kNN
    # (CPU only). Now the pipeline always attempts the cache; the
    # in-memory fallback returns ``(adata, None)`` only when no
    # on-disk source is available.
    _, ctx = _load_or_build_context(adata, config.atlas_name, outdir, args)

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


LONG_FORMAT_COLUMNS = [
    "Atlas Name",
    "Task",
    "Metric Name",
    "Embedding",
    "Reference/Input 1",
    "Prediction/Input 2",
    "Value",
    "Time (s)",
    "Fraction",
    "N_Cells",
    "Train_Domain",
    "Test_Domain",
]


def _empty_long_df() -> pd.DataFrame:
    return pd.DataFrame(columns=LONG_FORMAT_COLUMNS)


def _config_dict_for(
    config: SCFMConfig, adata: AnnData, scfm_emb: str
) -> dict:
    return {
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


def _should_run_cal_scfm(config: SCFMConfig) -> bool:
    """Return True if any of the four scfm-specific metric modules
    (scaling, stability, rare_types, cross_domain) should run on the
    live AnnData.

    The four modules are opt-in via the existing ``--scfm_*`` flags.
    ``cal_scfm`` is the only metric engine the scfm step runs from
    scratch; the per-task TSVs (cluster / annotation / dimred) are
    always loaded from disk.
    """
    return bool(
        config.scaling_fractions
        or config.n_seeds > 0
        or config.ref_label
        or config.predicted_label
        or config.domain_key
    )


def _apply_user_overrides(
    long: pd.DataFrame, kind: str, config: SCFMConfig
) -> pd.DataFrame:
    """Filter a long-format per-task DataFrame so only rows that
    match the user-specified ``--scfm_*`` keys survive.

    The per-task TSV is written by ``metrics.cal_*`` and may contain
    rows for every reference / predicted / batch / embedding column
    the column detector found. When the user explicitly named one
    column via ``--scfm_ref_label`` / ``--scfm_predicted_label`` /
    ``--scfm_batch_key`` / ``--scfm_embedding`` /
    ``--scfm_baseline_embeddings``, the scfm pipeline must consume
    *only* the rows that match those keys; the column detector's
    broader list is ignored for the scfm analysis (the underlying
    TSV is not modified — only the in-memory long table the
    diagnostic engine sees is filtered).

    The function is a no-op when the user did not set any
    ``--scfm_*`` key, so default (auto-detect) behaviour is
    preserved.
    """
    if kind == "annotation":
        keep = pd.Series(True, index=long.index)
        if config.ref_label:
            keep &= long["Reference/Input 1"].astype(str) == str(
                config.ref_label
            )
        if config.predicted_label:
            keep &= long["Prediction/Input 2"].astype(str) == str(
                config.predicted_label
            )
        if config.batch_key:
            mask = long["Metric Name"].astype(str).isin(
                {"kbet", "iLISI", "pcr", "graph_connectivity"}
            )
            keep &= (~mask) | (
                long["Prediction/Input 2"].astype(str) == str(config.batch_key)
            )
        return long[keep]
    if kind == "cluster":
        if not config.ref_label:
            return long
        keep = long["Prediction/Input 2"].astype(str) == str(config.ref_label)
        return long[keep]
    if kind == "dimred":
        allowed: set[str] = set()
        if config.scfm_embedding:
            allowed.add(str(config.scfm_embedding))
        for emb in config.baseline_embeddings or ():
            if emb:
                allowed.add(str(emb))
        if not allowed:
            return long
        keep = long["Embedding"].astype(str).isin(allowed)
        return long[keep]
    return long


def run_scfm_from_cache(
    adata: AnnData,
    config: SCFMConfig,
    args=None,
    *,
    outdir: Optional[str] = None,
) -> dict[str, Any]:
    """Run the scfm QC pipeline on one AnnData, loading per-task
    metric results from disk when they are available.

    This is the *load-then-diagnose* entry point. It does NOT
    re-compute the per-task metric engines (``cal_annot`` /
    ``cal_cluster`` / ``cal_dimred``). Those have already been run
    by the prior ``checkatlas metric_*`` (or Nextflow) invocation
    and their wide-format TSVs sit under
    ``checkatlas_files/<task>/<atlas>.tsv``.

    Behaviour:

    1. Determine the per-task TSV paths and call
       :func:`checkatlas.scfm.io.load_per_task_tsvs` to read them.
    2. Convert each wide TSV to long format with
       :func:`checkatlas.scfm.io.wide_to_long` and concatenate.
    3. If the scfm-specific metric modules are requested (scaling /
       stability / rare_types / cross_domain), call
       ``cal_scfm.cal_scfm`` on the live ``adata`` and append the
       resulting long-format rows.
    4. Call :func:`checkatlas.scfm.diagnostics.diagnose` to obtain
       the nine ``ProblemVerdict`` objects.
    5. Call :func:`checkatlas.scfm.report.write_all` to write the
       MultiQC-friendly TSVs to ``outdir``.

    Parameters
    ----------
    adata : AnnData
    config : SCFMConfig
    args : argparse.Namespace, optional
        Used for ``args.path`` and ``args.n_jobs``.
    outdir : str, optional
        Output directory. Defaults to
        ``checkatlas_files/scfm/<atlas_name>/`` under
        ``args.path``.

    Returns
    -------
    dict
        ``{"verdicts": [...], "metrics_df": ..., "out_paths": {...},
        "outdir": ...}``.
    """
    from . import io as scfm_io

    if outdir is None:
        base = getattr(args, "path", ".") if args is not None else "."
        folders.checkatlas_folders(base)
        outdir = os.path.join(
            folders.get_folder(base, folders.SCFM), config.atlas_name
        )
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

    # ── 1. Load per-task TSVs from disk ────────────────────────────
    base_path = getattr(args, "path", ".") if args is not None else "."
    tsvs = scfm_io.load_per_task_tsvs(base_path, config.atlas_name)

    long_frames: list[pd.DataFrame] = []
    for kind, df in tsvs.items():
        if df is None:
            logger.info("scfm: no %s TSV for %s", kind, config.atlas_name)
            continue
        long = scfm_io.wide_to_long(df, kind, config.atlas_name)
        if long.empty:
            continue
        long = _apply_user_overrides(long, kind, config_eff)
        if long.empty:
            logger.info(
                "scfm: %s TSV rows filtered out by user --scfm_* "
                "overrides for %s (the user named explicit ref / "
                "predicted / batch / embedding keys; the column "
                "detector's broader list is ignored for the scfm step)",
                kind,
                config.atlas_name,
            )
            continue
        long_frames.append(long)
        logger.info(
            "scfm: %s TSV -> %d long rows", kind, len(long)
        )

    # ── 2. Optionally run the four scfm-specific metric modules ───
    if _should_run_cal_scfm(config_eff):
        try:
            scfm_long = scfm_run.cal_scfm(
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
                n_seeds=config.n_seeds,
                noise_sigma=config.noise_sigma,
                min_domain_size=config.min_domain_size,
                n_jobs=n_jobs,
            )
            if scfm_long is not None and not scfm_long.empty:
                # Ensure all LONG_FORMAT_COLUMNS are present
                for col in LONG_FORMAT_COLUMNS:
                    if col not in scfm_long.columns:
                        scfm_long[col] = np.nan
                long_frames.append(scfm_long)
                logger.info(
                    "scfm: cal_scfm -> %d long rows", len(scfm_long)
                )
        except Exception as exc:
            logger.warning("scfm: cal_scfm failed: %s", exc)

    if long_frames:
        metrics_df = pd.concat(long_frames, ignore_index=True, sort=False)
        # Make sure the standard column set exists for downstream
        # consumers (verdicts / report / MultiQC).
        for col in LONG_FORMAT_COLUMNS:
            if col not in metrics_df.columns:
                metrics_df[col] = np.nan
    else:
        metrics_df = _empty_long_df()

    verdicts = diagnostics.diagnose(metrics_df, config_eff)

    config_dict = _config_dict_for(config_eff, adata, scfm_emb)
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
