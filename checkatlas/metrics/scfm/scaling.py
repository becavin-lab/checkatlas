"""Scaling-law saturation metric (Problem 3).

Reproduces the DenAdel et al. (2025) and Wang (2025) finding that
scFM embedding quality plateaus very quickly as the pretraining
cell-count grows.

**Contract (per the user instruction): the entire input atlas is
processed — no subsampling of the atlas. The scaling curve is a
property of the scFM pretraining corpus, not of the input atlas,
so the metric engine evaluates the *full atlas* once per (embedding,
metric) pair and emits a single value.**

The ``plateau_ratio`` row in the output is therefore populated
from a *user-supplied* pair of metric values (one from a small
atlas, one from a large atlas), not from internal subsampling. The
two values are passed via the ``small_value`` / ``large_value``
arguments of :func:`run`. When the user does not provide a pair,
the plateau row is omitted — the per-embedding metric on the full
atlas is still the primary result.

The implementation does NOT re-process adata: it relies on
``preprocess_context`` for precomputed kNN / distance matrices and
re-uses the existing silhouette / ARI / iLISI implementations.
When the context is missing, the function falls back to
pynndescent kNN (CPU) — but never re-loads the atlas, never
re-detects columns, and never re-computes distance matrices
that already exist in the cache.

The function returns a long-format DataFrame compatible with the
rest of the CheckAtlas metric pipeline (one row per ``(embedding,
metric)`` pair), plus an optional ``plateau_ratio`` row when the
user supplies small / large atlas values.
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


def _silhouette_value(
    X: np.ndarray,
    labels: np.ndarray,
    precomputed_dists=None,
    n_jobs: int = -1,
) -> float:
    """ASW on the full atlas using precomputed distances if available."""
    from ..cluster import silhouette

    try:
        return float(
            silhouette.run(
                X,
                labels,
                n_jobs=n_jobs,
                verbose=False,
                precomputed_dists=precomputed_dists,
            )
        )
    except Exception as exc:  # pragma: no cover
        logger.debug("scaling.silhouette failed: %s", exc)
        return float("nan")


def _ari_value(ref: np.ndarray, pred: np.ndarray) -> float:
    """Adjusted Rand Index on the full atlas — no re-computation."""
    from sklearn.metrics import adjusted_rand_score

    try:
        return float(adjusted_rand_score(ref.astype(str), pred.astype(str)))
    except Exception as exc:  # pragma: no cover
        logger.debug("scaling.ari failed: %s", exc)
        return float("nan")


def _ilisi_value(nn, labels: np.ndarray) -> float:
    """iLISI on the full atlas using precomputed kNN."""
    from ..annot import lisi

    try:
        return float(lisi.run_with_neighbors(nn, labels, verbose=False))
    except Exception as exc:  # pragma: no cover
        logger.debug("scaling.ilisi failed: %s", exc)
        return float("nan")


def _load_precomputed(
    preprocess_context,
    embeddings: Sequence[str],
    n_total: int,
) -> tuple[dict, dict]:
    """Load precomputed kNN graphs and distance matrices from
    ``preprocess_context``.

    Returns
    -------
    (emb_nn, precomputed_dists)
        Two dicts keyed by embedding name. Empty if the context
        has no precomputed artefacts for an embedding.
    """
    emb_nn: dict = {}
    precomputed_dists: dict = {}
    if preprocess_context is None:
        return emb_nn, precomputed_dists

    if hasattr(preprocess_context, "knn_paths"):
        from .._cache import load_knn
        from .._neighbors import NeighborResults

        for emb, npz_path in preprocess_context.knn_paths.items():
            try:
                loaded = load_knn(
                    os.path.dirname(npz_path),
                    os.path.splitext(os.path.basename(npz_path))[0],
                )
                if loaded is not None:
                    emb_nn[emb] = NeighborResults(
                        indices=loaded[0], distances=loaded[1]
                    )
            except Exception:
                continue

    if hasattr(preprocess_context, "cluster_dir"):

        def _safe(s: str) -> str:
            return s.replace("/", "_").replace(" ", "_")

        for emb in embeddings:
            tri_path = os.path.join(
                preprocess_context.cluster_dir, f"dist_{_safe(emb)}.tri"
            )
            npy_path = tri_path.replace(".tri", ".npy")
            if os.path.exists(tri_path):
                try:
                    from .._triangular import TriangularMatrix

                    precomputed_dists[emb] = TriangularMatrix(
                        n=n_total, filepath=tri_path, mode="r"
                    )
                except Exception:
                    continue
            elif os.path.exists(npy_path):
                try:
                    precomputed_dists[emb] = np.load(npy_path)
                except Exception:
                    continue

    return emb_nn, precomputed_dists


def _resolve_labels(
    adata: AnnData,
    ref_label: Optional[str],
    pred_label: Optional[str],
) -> Optional[np.ndarray]:
    """Return the obs column to use for ASW / Dunn-style metrics.

    Prefers the ref label, then the pred label. Returns ``None`` if
    neither is available so the caller can skip the silhouette row.
    """
    for key in (ref_label, pred_label):
        if key and key in adata.obs.columns:
            return np.asarray(adata.obs[key])
    return None


def run(
    adata: AnnData,
    embeddings: Sequence[str],
    *,
    ref_label: Optional[str] = None,
    pred_label: Optional[str] = None,
    batch_key: Optional[str] = None,
    metric_subset: Sequence[str] = ("silhouette", "ari", "ilisi"),
    n_jobs: int = -1,
    preprocess_context=None,
    plateau_pairs: Optional[dict] = None,
) -> pd.DataFrame:
    """Evaluate ``metric_subset`` on the **full atlas** for each embedding.

    The full input ``adata`` is consumed as-is — no subsampling, no
    re-processing, no re-detection of columns. Precomputed kNN graphs
    and distance matrices from ``preprocess_context`` are reused
    without modification. When the context is missing, the
    silhouette fallback falls back to ``pynndescent`` (CPU) for the
    distance matrix.

    Parameters
    ----------
    adata : AnnData
        Input atlas. Used as-is.
    embeddings : sequence of str
        Embedding keys to evaluate.
    ref_label, pred_label, batch_key : str, optional
        Annotation / batch columns. Silhouette requires one of
        ``ref_label`` / ``pred_label``; ARI requires both; iLISI
        requires ``batch_key``.
    metric_subset : sequence of str
        Subset of ``{"silhouette", "ari", "ilisi"}``.
    n_jobs : int
        Parallel jobs.
    preprocess_context : PreprocessContext, optional
        Source of precomputed kNN / distance matrices.
    plateau_pairs : dict, optional
        ``{embedding: {metric: (small_value, large_value)}}`` for
        external scaling-curve computation. When provided, a
        ``plateau_ratio = large / small`` row is appended per
        ``(embedding, metric)``.

    Returns
    -------
    pd.DataFrame
        Long-format with columns
        ``[Embedding, Metric, Fraction, N_Cells, Value, Time (s)]``.
        One row per ``(embedding, metric)`` pair on the full atlas
        (``Fraction=1.0, N_Cells=adata.n_obs``). When ``plateau_pairs``
        is provided, one extra row per pair is appended
        (``Fraction=-1.0``).
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
    emb_nn, precomputed_dists = _load_precomputed(
        preprocess_context, embeddings, n_total
    )

    labels = _resolve_labels(adata, ref_label, pred_label)
    ref_arr = (
        np.asarray(adata.obs[ref_label].astype(str))
        if ref_label and ref_label in adata.obs.columns
        else None
    )
    pred_arr = (
        np.asarray(adata.obs[pred_label].astype(str))
        if pred_label and pred_label in adata.obs.columns
        else None
    )
    batch_arr = (
        np.asarray(adata.obs[batch_key])
        if batch_key and batch_key in adata.obs.columns
        else None
    )

    rows: list[dict] = []
    for emb in embeddings:
        if emb not in adata.obsm:
            logger.debug("scaling: embedding %s not in adata.obsm", emb)
            continue
        X_full = adata.obsm[emb]
        for metric in metric_subset:
            t0 = time.time()
            if metric == "silhouette":
                if labels is None:
                    continue
                pdists = precomputed_dists.get(emb)
                val = _silhouette_value(X_full, labels, pdists, n_jobs)
            elif metric == "ari":
                if ref_arr is None or pred_arr is None:
                    continue
                val = _ari_value(ref_arr, pred_arr)
            elif metric == "ilisi":
                if batch_arr is None or emb not in emb_nn:
                    continue
                val = _ilisi_value(emb_nn[emb], batch_arr)
            else:
                continue
            rows.append(
                {
                    "Embedding": emb,
                    "Metric": metric,
                    "Fraction": 1.0,
                    "N_Cells": n_total,
                    "Value": val,
                    "Time (s)": round(time.time() - t0, 3),
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    if plateau_pairs:
        plateau_rows: list[dict] = []
        for emb, by_metric in plateau_pairs.items():
            for metric, (small_v, large_v) in by_metric.items():
                if abs(small_v) <= 1e-9:
                    ratio = float("nan")
                else:
                    try:
                        ratio = float(large_v) / float(small_v)
                    except (TypeError, ValueError):
                        ratio = float("nan")
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
