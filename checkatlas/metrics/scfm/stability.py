"""Seed-stability metric (Problem 8).

Reproduces the Liu 2024 and Xu 2026 finding that scFMs are unstable
across random subsamples (high coefficient of variation).

The implementation does NOT re-process adata: it relies on the
``preprocess_context`` for precomputed kNN / distance matrices and
re-uses the existing silhouette / ARI / iLISI implementations.
"""

from __future__ import annotations

import logging
import time
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score

logger = logging.getLogger("checkatlas")


def _silhouette_value(
    X: np.ndarray,
    labels: np.ndarray,
    precomputed_dists=None,
    n_jobs: int = -1,
) -> float:
    from ..cluster import silhouette

    try:
        return float(
            silhouette.run(
                X, labels, n_jobs=n_jobs, verbose=False, precomputed_dists=precomputed_dists
            )
        )
    except Exception:
        return float("nan")


def _ari_value(ref: np.ndarray, pred: np.ndarray) -> float:
    try:
        return float(adjusted_rand_score(ref.astype(str), pred.astype(str)))
    except Exception:
        return float("nan")


def _ilisi_value(nn, labels: np.ndarray) -> float:
    from ..annot import lisi

    try:
        return float(lisi.run_with_neighbors(nn, labels, verbose=False))
    except Exception:
        return float("nan")


def _kbet_value(nn, labels: np.ndarray) -> float:
    from ..annot import kbet

    try:
        k = min(25, nn.n_neighbors)
        sub = nn.subset_neighbors(k)
        return float(kbet.run_with_neighbors(sub, labels, verbose=False))
    except Exception:
        return float("nan")


def run(
    adata: AnnData,
    embeddings: Sequence[str],
    *,
    ref_label: Optional[str] = None,
    pred_label: Optional[str] = None,
    batch_key: Optional[str] = None,
    n_seeds: int = 5,
    subsample_frac: float = 0.90,
    base_seed: int = 0,
    metrics: Sequence[str] = ("silhouette", "ari", "ilisi", "kbet"),
    n_jobs: int = -1,
    preprocess_context=None,
) -> pd.DataFrame:
    """Re-evaluate ``metrics`` on ``n_seeds`` random subsamples.

    For each embedding, the function computes ``mean``, ``std``, and
    ``CV = std / |mean|`` for every requested metric. A high ``CV``
    indicates the embedding is sensitive to the subsample.

    Parameters
    ----------
    adata : AnnData
    embeddings : sequence of str
    ref_label, pred_label, batch_key : str, optional
        Required for ARI / iLISI / kBET respectively.
    n_seeds : int
        Number of subsamples to draw.
    subsample_frac : float
        Fraction of cells in each subsample (0 < f <= 1).
    metrics : sequence of str
        Subset of ``{"silhouette", "ari", "ilisi", "kbet"}``.

    Returns
    -------
    pd.DataFrame
        Long format with columns
        ``[Embedding, Metric, Statistic, Value]`` where ``Statistic``
        is one of ``mean``, ``std``, ``cv``.
    """
    if not embeddings or n_seeds <= 0:
        return pd.DataFrame(
            columns=["Embedding", "Metric", "Statistic", "Value"]
        )

    n_total = adata.n_obs
    n_sub = max(50, int(n_total * subsample_frac))
    n_sub = min(n_sub, n_total)

    emb_nn: dict = {}
    precomputed_dists: dict = {}
    if preprocess_context is not None and hasattr(preprocess_context, "knn_paths"):
        import os as _os

        from .._cache import load_knn
        from .._neighbors import NeighborResults

        for emb, npz_path in preprocess_context.knn_paths.items():
            try:
                loaded = load_knn(
                    _os.path.dirname(npz_path),
                    _os.path.splitext(_os.path.basename(npz_path))[0],
                )
                if loaded is not None:
                    emb_nn[emb] = NeighborResults(
                        indices=loaded[0], distances=loaded[1]
                    )
            except Exception:
                pass
        if hasattr(preprocess_context, "cluster_dir"):
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

    rng = np.random.default_rng(base_seed)
    seeds = [int(rng.integers(0, 1 << 31)) for _ in range(n_seeds)]

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
            continue
        X_full = adata.obsm[emb]
        per_metric_values: dict[str, list[float]] = {m: [] for m in metrics}
        for s in seeds:
            t0 = time.time()
            sub_rng = np.random.default_rng(s)
            if n_sub >= n_total:
                idx = np.arange(n_total)
            else:
                idx = sub_rng.choice(n_total, size=n_sub, replace=False)
            X_sub = X_full[idx]
            sub_ref = ref_arr[idx] if ref_arr is not None else None
            sub_pred = pred_arr[idx] if pred_arr is not None else None
            sub_batch = batch_arr[idx] if batch_arr is not None else None

            for metric in metrics:
                if metric == "silhouette":
                    if sub_ref is None and sub_pred is None:
                        continue
                    labels = sub_ref if sub_ref is not None else sub_pred
                    pdists = precomputed_dists.get(emb)
                    if pdists is not None and hasattr(pdists, "_filepath"):
                        # TriangularMatrix: cannot subsample the full
                        # distance matrix cheaply, so we recompute on
                        # the subsampled X (this is per-seed and
                        # necessary — it is NOT re-processing of the
                        # original adata).
                        pdists = None
                    val = _silhouette_value(X_sub, labels, pdists, n_jobs)
                elif metric == "ari":
                    if sub_ref is None or sub_pred is None:
                        continue
                    val = _ari_value(sub_ref, sub_pred)
                elif metric == "ilisi":
                    if sub_batch is None or emb not in emb_nn:
                        continue
                    full_nn = emb_nn[emb]
                    sub_nn = full_nn.subset_neighbors(
                        min(90, full_nn.n_neighbors)
                    )
                    val = _ilisi_value(sub_nn, sub_batch)
                elif metric == "kbet":
                    if sub_batch is None or emb not in emb_nn:
                        continue
                    full_nn = emb_nn[emb]
                    sub_nn = full_nn.subset_neighbors(
                        min(25, full_nn.n_neighbors)
                    )
                    val = _kbet_value(sub_nn, sub_batch)
                else:
                    continue
                per_metric_values[metric].append(val)
            _ = time.time() - t0

        for metric, values in per_metric_values.items():
            if not values:
                continue
            arr = np.asarray(values, dtype=float)
            arr = arr[~np.isnan(arr)]
            if arr.size == 0:
                continue
            mean = float(np.mean(arr))
            std = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
            cv = float(std / abs(mean)) if abs(mean) > 1e-9 else 0.0
            for stat, val in (("mean", mean), ("std", std), ("cv", cv)):
                rows.append(
                    {
                        "Embedding": emb,
                        "Metric": metric,
                        "Statistic": stat,
                        "Value": val,
                    }
                )
    return pd.DataFrame(rows)
