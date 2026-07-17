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

For large atlases (>300k cells), core-set subsampling is applied
via the same ``_stratified_subsample`` protocol used by the main
pipeline. The subsample indices from ``preprocess_context`` are
reused when available; otherwise a uniform subsample to
``_CORESET_SIZE`` (50k) is used.

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

_LARGE_ATLAS_THRESHOLD = 300_000
_CORESET_SIZE = 50_000
_CORESET_MIN_PER_CLUSTER = 100
_CORESET_SEED = 42


def _resolve_cell_limit(
    adata: AnnData,
    preprocess_context,
) -> tuple[AnnData, Optional[np.ndarray], int]:
    """Apply core-set subsampling for large atlases (>300k cells).

    Returns ``(adata, subsample_indices, n_cells)`` where
    ``adata`` is unchanged (the atlas is never modified in-place),
    ``subsample_indices`` is the array of cell indices to use (or
    ``None`` if subsampling is not needed), and ``n_cells`` is the
    number of cells after subsampling.

    When ``preprocess_context`` carries precomputed
    ``subsample_indices`` those are reused directly.
    """
    n_cells = adata.n_obs
    if n_cells <= _LARGE_ATLAS_THRESHOLD:
        return adata, None, n_cells

    if (
        preprocess_context is not None
        and hasattr(preprocess_context, "subsample_indices")
        and preprocess_context.subsample_indices is not None
    ):
        logger.info(
            "scfm: reusing precomputed core-set subsample (%d cells) from preprocess_context",
            len(preprocess_context.subsample_indices),
        )
        return adata, preprocess_context.subsample_indices, len(
            preprocess_context.subsample_indices
        )

    rng = np.random.default_rng(_CORESET_SEED)
    indices = rng.choice(n_cells, size=_CORESET_SIZE, replace=False)
    indices = np.sort(indices)
    logger.info(
        "scfm: large atlas (%d cells) — applying uniform core-set subsample to %d cells",
        n_cells,
        _CORESET_SIZE,
    )
    return adata, indices, _CORESET_SIZE


def _subsample_array(arr: np.ndarray, indices: Optional[np.ndarray]) -> np.ndarray:
    """Return ``arr[indices]`` if indices are given, else ``arr`` unchanged."""
    if indices is None:
        return arr
    return arr[indices]


def _subsample_embedding(
    adata: AnnData, embedding_key: str, indices: Optional[np.ndarray]
) -> np.ndarray:
    """Return the embedding array, optionally subsampled."""
    X = adata.obsm[embedding_key]
    return _subsample_array(X, indices)


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

    For atlases larger than 300k cells, core-set subsampling is
    applied so the scaling and stability modules operate on a
    manageable subset. The subsample indices from
    ``preprocess_context`` are reused when available; otherwise a
    uniform random sample of 50k cells is drawn. The rare_types
    and cross_domain modules always run on the full atlas (they
    are not O(N^2) operations).

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
    n_seeds : int
        Number of perturbation seeds for the stability metric.
    noise_sigma : float
        Noise level for the stability embedding-perturbation
        (sigma in units of per-feature std). Default 0.10.
    min_domain_size : int
    n_jobs : int
    preprocess_context : PreprocessContext, optional
        Re-uses precomputed kNN / distance matrices. Built or loaded
        by ``scfm.pipeline._load_or_build_context``.

    Returns
    -------
    pd.DataFrame
        Long-format metric table.
    """
    from . import cross_domain, rare_types, scaling, stability

    # Apply core-set subsampling for large atlases (>300k cells)
    _, subset_indices, n_eff = _resolve_cell_limit(adata, preprocess_context)

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
                n_jobs=n_jobs,
                preprocess_context=preprocess_context,
                subset_indices=subset_indices,
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
                sigma_scale=noise_sigma,
                n_jobs=n_jobs,
                preprocess_context=preprocess_context,
                subset_indices=subset_indices,
            )
            for _, r in df.iterrows():
                rows.append(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "scfm_stability",
                        "Metric Name": f"{r['Metric']}_{r['Statistic']}",
                        "Embedding": str(r["Embedding"]),
                        "Fraction": np.nan,
                        "N_Cells": n_eff,
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
                row = {
                    "Atlas Name": atlas_name,
                    "Task": "scfm_cross_domain",
                    "Metric Name": "cross_domain_ari",
                    "Embedding": str(r["Embedding"]),
                    "Fraction": np.nan,
                    # Per the user contract, N_Cells in the output
                    # TSV is always the full-atlas size — the
                    # per-domain N_Test lives in a separate column
                    # (or is dropped) so downstream consumers
                    # always see a consistent "size of the
                    # evaluated atlas".
                    "N_Cells": adata.n_obs,
                    "Value": float(r["ARI"])
                    if pd.notna(r["ARI"])
                    else np.nan,
                    "Time (s)": 0.0,
                }
                if "Test_Domain" in r.index:
                    row["Test_Domain"] = str(r["Test_Domain"])
                if "Train_Domain" in r.index:
                    row["Train_Domain"] = str(r["Train_Domain"])
                if "N_Test" in r.index:
                    row["N_Test_Domain"] = int(r["N_Test"])
                rows.append(row)
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
