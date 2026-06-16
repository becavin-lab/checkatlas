"""Scaling-law saturation metric (Problem 3).

Reproduces the DenAdel et al. (2025) and Wang (2025) finding that scFM
embedding quality plateaus very quickly as the cell-count grows.

The implementation does NOT re-process adata: it relies on
``preprocess_context`` for precomputed kNN / distance matrices and
re-uses the existing silhouette / ARI / iLISI implementations.

The function returns a long-format DataFrame compatible with the rest
of the CheckAtlas metric pipeline (one row per ``(embedding, fraction,
metric_subset_item)`` triple), plus a derived ``plateau_ratio``.
"""

from __future__ import annotations

import logging
import os
import time
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger("checkatlas")


def _silhouette_with_precomputed(
    adata: AnnData,
    emb: str,
    label: Optional[str],
    precomputed_dists: Optional[dict],
    n_jobs: int,
) -> float:
    """ASW using precomputed distances if available (no re-computation)."""
    from ..cluster import silhouette

    X = adata.obsm[emb]
    if label is not None and label in adata.obs.columns:
        labels = np.asarray(adata.obs[label])
    else:
        return float("nan")
    pdists = precomputed_dists.get(emb) if precomputed_dists else None
    try:
        return float(
            silhouette.run(
                X,
                labels,
                n_jobs=n_jobs,
                verbose=False,
                precomputed_dists=pdists,
            )
        )
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("scaling.silhouette failed for %s/%s: %s", emb, label, exc)
        return float("nan")


def _ari_with_precomputed(
    adata: AnnData,
    emb: str,
    ref_label: Optional[str],
    pred_label: Optional[str],
) -> float:
    """Adjusted Rand Index between two obs columns, no re-computation."""
    from sklearn.metrics import adjusted_rand_score

    if ref_label is None or pred_label is None:
        return float("nan")
    if ref_label not in adata.obs or pred_label not in adata.obs:
        return float("nan")
    try:
        return float(
            adjusted_rand_score(
                adata.obs[ref_label].astype(str),
                adata.obs[pred_label].astype(str),
            )
        )
    except Exception as exc:  # pragma: no cover
        logger.debug("scaling.ari failed: %s", exc)
        return float("nan")


def _ilisi_with_precomputed(
    adata: AnnData,
    emb: str,
    batch_key: Optional[str],
    emb_nn: dict,
) -> float:
    """iLISI using precomputed kNN (no re-computation)."""
    from ..annot import lisi

    if batch_key is None or batch_key not in adata.obs.columns:
        return float("nan")
    nn = emb_nn.get(emb)
    if nn is None:
        return float("nan")
    try:
        return float(
            lisi.run_with_neighbors(
                nn, adata.obs[batch_key], verbose=False
            )
        )
    except Exception as exc:  # pragma: no cover
        logger.debug("scaling.ilisi failed for %s/%s: %s", emb, batch_key, exc)
        return float("nan")


def run(
    adata: AnnData,
    embeddings: Sequence[str],
    *,
    ref_label: Optional[str] = None,
    pred_label: Optional[str] = None,
    batch_key: Optional[str] = None,
    fractions: Sequence[float] = (0.01, 0.10, 0.50, 1.00),
    metric_subset: Sequence[str] = ("silhouette", "ari", "ilisi"),
    seed: int = 42,
    n_jobs: int = -1,
    preprocess_context=None,
) -> pd.DataFrame:
    """Re-evaluate ``metric_subset`` on subsampled adata at each fraction.

    Parameters
    ----------
    adata : AnnData
        Input atlas.
    embeddings : sequence of str
        Embedding keys to evaluate (e.g. ``["X_geneformer", "X_pca"]``).
    ref_label, pred_label : str, optional
        Annotation columns for ARI. Either may be None.
    batch_key : str, optional
        Batch column for iLISI.
    fractions : sequence of float
        Subsample fractions. Default ``(0.01, 0.10, 0.50, 1.00)``.
    metric_subset : sequence of str
        Subset of ``{"silhouette", "ari", "ilisi"}``.
    seed : int
        RNG seed for the subsample.
    n_jobs : int
        Parallel jobs.
    preprocess_context : PreprocessContext, optional
        If provided, precomputed kNN graphs (``emb_nn``) and distance
        matrices (``precomputed_dists``) are reused and the adata is
        NEVER re-processed.

    Returns
    -------
    pd.DataFrame
        Long-format with columns
        ``[Embedding, Metric, Fraction, N_Cells, Value, Time (s)]`` and
        a derived ``plateau_ratio`` row per (embedding, metric).
    """
    if not embeddings:
        return pd.DataFrame(
            columns=[
                "Embedding",
                "Metric",
                "Fraction",
                "N_Cells",
                "Value",
                "Time (s)",
            ]
        )

    n_total = adata.n_obs
    rng = np.random.default_rng(seed)

    emb_nn: dict = {}
    precomputed_dists: dict = {}
    if preprocess_context is not None:
        if hasattr(preprocess_context, "knn_paths"):
            for emb, npz_path in preprocess_context.knn_paths.items():
                try:
                    from .._cache import load_knn
                    from .._neighbors import NeighborResults

                    loaded = load_knn(
                        os.path.dirname(npz_path),
                        os.path.splitext(os.path.basename(npz_path))[0],
                    )
                    if loaded is not None:
                        emb_nn[emb] = NeighborResults(
                            indices=loaded[0], distances=loaded[1]
                        )
                except Exception:
                    pass
        if hasattr(preprocess_context, "cluster_dir"):
            import os as _os

            def _safe(s: str) -> str:
                return s.replace("/", "_").replace(" ", "_")

            for emb in embeddings:
                tri_path = _os.path.join(
                    preprocess_context.cluster_dir, f"dist_{_safe(emb)}.tri"
                )
                npy_path = tri_path.replace(".tri", ".npy")
                if _os.path.exists(tri_path):
                    try:
                        from .._triangular import TriangularMatrix

                        precomputed_dists[emb] = TriangularMatrix(
                            n=n_total, filepath=tri_path, mode="r"
                        )
                    except Exception:
                        pass
                elif _os.path.exists(npy_path):
                    try:
                        precomputed_dists[emb] = np.load(npy_path)
                    except Exception:
                        pass

    rows: list[dict] = []
    for emb in embeddings:
        if emb not in adata.obsm:
            logger.debug("scaling: embedding %s not in adata.obsm", emb)
            continue
        for fraction in sorted(set(float(f) for f in fractions)):
            t0 = time.time()
            n_sub = max(50, int(n_total * fraction))
            n_sub = min(n_sub, n_total)
            if n_sub >= n_total:
                sub = adata
            else:
                idx = rng.choice(n_total, size=n_sub, replace=False)
                sub = adata[idx].copy()
            for metric in metric_subset:
                if metric == "silhouette":
                    val = _silhouette_with_precomputed(
                        sub, emb, ref_label or pred_label, precomputed_dists, n_jobs
                    )
                elif metric == "ari":
                    val = _ari_with_precomputed(sub, emb, ref_label, pred_label)
                elif metric == "ilisi":
                    val = _ilisi_with_precomputed(sub, emb, batch_key, emb_nn)
                else:
                    continue
                rows.append(
                    {
                        "Embedding": emb,
                        "Metric": metric,
                        "Fraction": fraction,
                        "N_Cells": n_sub,
                        "Value": val,
                        "Time (s)": round(time.time() - t0, 3),
                    }
                )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    plateau_rows: list[dict] = []
    for (emb, metric), group in df.groupby(["Embedding", "Metric"]):
        sub = group.dropna(subset=["Value"])
        if sub.empty:
            continue
        sub = sub.sort_values("Fraction")
        full = sub[sub["Fraction"] >= 0.999]["Value"]
        small = sub[sub["Fraction"] <= 0.011]["Value"]
        if full.empty or small.empty:
            continue
        full_v = float(full.iloc[0])
        small_v = float(small.iloc[0])
        ratio = float(full_v) / float(small_v) if abs(small_v) > 1e-9 else float("nan")
        plateau_rows.append(
            {
                "Embedding": emb,
                "Metric": metric,
                "Fraction": -1.0,
                "N_Cells": n_total,
                "Value": ratio,
                "Time (s)": 0.0,
            }
        )
    if plateau_rows:
        df = pd.concat([df, pd.DataFrame(plateau_rows)], ignore_index=True)
    return df
