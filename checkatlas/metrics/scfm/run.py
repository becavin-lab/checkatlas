"""Orchestrator for scFM-specific metrics.

This module is the **only** entry point for the scfm Layer-1
pipeline. It follows Becavin comment 2 strictly:

  * No re-loading of adata.
  * No re-detection of columns.
  * No re-computation of kNN graphs or distance matrices that
    ``preprocess_context`` already provides.
  * The new scfm modules are run **after** the existing
    ``cal_cluster`` / ``cal_annot`` / ``cal_dimred`` calls, using
    the same precomputed artefacts.

The function returns a long-format DataFrame with the same column
schema as the existing ``cal_*`` pipelines (``[Atlas Name, Metric
Name, Embedding, ..., Value, Time (s)]``) so the result can be
written into the same TSV files that MultiQC reads.
"""

from __future__ import annotations

import logging
import time
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger("checkatlas")


def cal_scfm(
    adata: AnnData,
    *,
    scfm_embedding: str,
    baseline_embeddings: Sequence[str] = (),
    ref_label: Optional[str] = None,
    pred_label: Optional[str] = None,
    batch_key: Optional[str] = None,
    domain_key: Optional[str] = None,
    patient_key: Optional[str] = None,
    outcome_key: Optional[str] = None,
    atlas_name: str = "",
    scaling_fractions: Sequence[float] = (0.01, 0.10, 0.50, 1.00),
    n_seeds: int = 5,
    noise_sigma: float = 0.10,
    min_domain_size: int = 50,
    n_jobs: int = -1,
    preprocess_context=None,
    run_scaling: bool = True,
    run_stability: bool = True,
    run_rare_types: bool = True,
    run_cross_domain: bool = True,
) -> pd.DataFrame:
    """Run all four scfm-specific metric modules.

    Parameters
    ----------
    adata : AnnData
    scfm_embedding : str
        Embedding key for the scFM under test.
    baseline_embeddings : sequence of str
        Baseline embedding keys (e.g. ``["X_pca", "X_scvi"]``).
    ref_label, pred_label, batch_key, domain_key, patient_key, outcome_key : str, optional
        Annotation / metadata columns. Modules that require missing
        columns are skipped silently.
    atlas_name : str
    scaling_fractions : sequence of float
    n_seeds : int
    noise_sigma : float
        Reserved for the gaussian-robustness function in
        ``information_loss``. Currently not used by the four modules
        kept after Becavin comment 7; kept for API symmetry.
    min_domain_size : int
    n_jobs : int
    preprocess_context : PreprocessContext, optional
        Re-uses precomputed kNN / distance matrices.

    Returns
    -------
    pd.DataFrame
        Long-format metric table.
    """
    from . import cross_domain, rare_types, scaling, stability

    all_embeddings: list[str] = []
    if scfm_embedding in adata.obsm:
        all_embeddings.append(scfm_embedding)
    for emb in baseline_embeddings:
        if emb in adata.obsm and emb not in all_embeddings:
            all_embeddings.append(emb)

    rows: list[dict] = []

    if run_scaling and all_embeddings:
        t0 = time.time()
        try:
            df = scaling.run(
                adata,
                all_embeddings,
                ref_label=ref_label,
                pred_label=pred_label,
                batch_key=batch_key,
                fractions=scaling_fractions,
                n_jobs=n_jobs,
                preprocess_context=preprocess_context,
            )
            for _, r in df.iterrows():
                rows.append(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "scfm_scaling",
                        "Metric Name": str(r["Metric"]),
                        "Embedding": str(r["Embedding"]),
                        "Fraction": float(r["Fraction"]),
                        "N_Cells": int(r["N_Cells"]),
                        "Value": float(r["Value"])
                        if pd.notna(r["Value"])
                        else np.nan,
                        "Time (s)": float(r["Time (s)"]),
                    }
                )
        except Exception as exc:
            logger.warning("cal_scfm: scaling failed: %s", exc)
        _ = time.time() - t0

    if run_stability and all_embeddings:
        try:
            df = stability.run(
                adata,
                all_embeddings,
                ref_label=ref_label,
                pred_label=pred_label,
                batch_key=batch_key,
                n_seeds=n_seeds,
                n_jobs=n_jobs,
                preprocess_context=preprocess_context,
            )
            for _, r in df.iterrows():
                rows.append(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "scfm_stability",
                        "Metric Name": f"{r['Metric']}_{r['Statistic']}",
                        "Embedding": str(r["Embedding"]),
                        "Fraction": np.nan,
                        "N_Cells": adata.n_obs,
                        "Value": float(r["Value"])
                        if pd.notna(r["Value"])
                        else np.nan,
                        "Time (s)": 0.0,
                    }
                )
        except Exception as exc:
            logger.warning("cal_scfm: stability failed: %s", exc)

    if run_rare_types and pred_label and ref_label:
        try:
            df = rare_types.run(
                adata,
                all_embeddings,
                pred_label=pred_label,
                ref_label=ref_label,
            )
            for _, r in df.iterrows():
                rows.append(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "scfm_rare_types",
                        "Metric Name": str(r["Metric"]),
                        "Embedding": str(r["Embedding"]),
                        "Fraction": np.nan,
                        "N_Cells": adata.n_obs,
                        "Value": float(r["Value"])
                        if pd.notna(r["Value"])
                        else np.nan,
                        "Time (s)": 0.0,
                    }
                )
        except Exception as exc:
            logger.warning("cal_scfm: rare_types failed: %s", exc)

    if run_cross_domain and domain_key and ref_label and scfm_embedding in adata.obsm:
        try:
            df = cross_domain.run(
                adata,
                scfm_embedding,
                ref_label=ref_label,
                domain_key=domain_key,
                min_domain_size=min_domain_size,
            )
            for _, r in df.iterrows():
                rows.append(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "scfm_cross_domain",
                        "Metric Name": "cross_domain_ari",
                        "Embedding": str(r["Embedding"]),
                        "Train_Domain": str(r["Train_Domain"]),
                        "Test_Domain": str(r["Test_Domain"]),
                        "Fraction": np.nan,
                        "N_Cells": int(r["N_Test"]),
                        "Value": float(r["ARI"])
                        if pd.notna(r["ARI"])
                        else np.nan,
                        "Time (s)": 0.0,
                    }
                )
        except Exception as exc:
            logger.warning("cal_scfm: cross_domain failed: %s", exc)

    if not rows:
        return pd.DataFrame(
            columns=[
                "Atlas Name",
                "Task",
                "Metric Name",
                "Embedding",
                "Value",
                "Time (s)",
            ]
        )
    return pd.DataFrame(rows)
