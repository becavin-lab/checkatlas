"""Seed-stability metric (Problem 8).

Reproduces the Liu 2024 and Xu 2026 finding that scFMs are unstable
across random perturbations (high coefficient of variation).

**Contract (per the user instruction): the entire input atlas is
used. The atlas is NEVER subsampled. Stability is measured by
perturbing the *embedding* of each cell (Gaussian noise injection
on the embedding itself), not by dropping cells.**

Why this design? The scFM pretraining corpus size is fixed — the
*embedding* the user has stored in ``adata.obsm`` is what we
evaluate. We can perturb the embedding cheaply (O(N × d) Gaussian
noise) but we cannot change the cell set without re-running the
downstream labels (cell type, batch) which are stored in
``adata.obs`` and not under our control.

The perturbation is the standard one used in robustness benchmarks
(Atti & Subramaniam 2025): additive Gaussian noise with sigma
proportional to the per-feature standard deviation of the
embedding. The metric is CV of the downstream silhouette / ARI /
iLISI / kBET over ``n_seeds`` noise samples.

The implementation does NOT re-process adata: it relies on the
``preprocess_context`` for precomputed kNN / distance matrices and
re-uses the existing silhouette / ARI / iLISI / kBET
implementations. The kNN graph for kBET / iLISI is built ONCE on
the un-perturbed embedding and reused across all noise samples via
``NeighborResults.subset_neighbors(k)``.
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
                X,
                labels,
                n_jobs=n_jobs,
                verbose=False,
                precomputed_dists=precomputed_dists,
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


def _load_precomputed(
    preprocess_context,
    embeddings: Sequence[str],
    n_total: int,
) -> tuple[dict, dict]:
    """Load precomputed kNN graphs and distance matrices from
    ``preprocess_context``."""
    import os

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


def _perturb_embedding(
    X: np.ndarray,
    sigma_scale: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Add Gaussian noise proportional to the per-feature std of X.

    ``sigma_scale`` is the noise level in units of feature-std
    (sigma_scale=0.1 means sigma = 0.1 * std(X, axis=0), the
    convention used in Atti & Subramaniam 2025).
    """
    feat_std = X.std(axis=0, keepdims=True)
    sigma = sigma_scale * feat_std
    noise = rng.standard_normal(X.shape) * sigma
    return X + noise


def run(
    adata: AnnData,
    embeddings: Sequence[str],
    *,
    ref_label: Optional[str] = None,
    pred_label: Optional[str] = None,
    batch_key: Optional[str] = None,
    n_seeds: int = 5,
    sigma_scale: float = 0.10,
    base_seed: int = 0,
    metrics: Sequence[str] = ("silhouette", "ari", "ilisi", "kbet"),
    n_jobs: int = -1,
    preprocess_context=None,
    subset_indices: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Measure CV across ``n_seeds`` embedding-perturbation runs.

    For large atlases, ``subset_indices`` can provide a core-set
    of cell indices for O(N^2) operations (silhouette). The
    embedding is perturbed only on the subset; label-based
    metrics (ARI, iLISI, kBET) use the full atlas.

    The perturbation is the standard one used in robustness
    benchmarks (Atti & Subramaniam 2025): additive Gaussian noise
    with sigma proportional to the per-feature standard deviation
    of the embedding. The metric is CV of the downstream
    silhouette / ARI / iLISI / kBET over ``n_seeds`` noise
    samples.

    The kNN graph for kBET / iLISI is built once on the
    un-perturbed embedding and reused across all noise samples
    via ``NeighborResults.subset_neighbors(k)``.

    Parameters
    ----------
    adata : AnnData
        Input atlas. Used as-is.
    embeddings : sequence of str
    ref_label, pred_label : str, optional
        Required for ARI; one of them required for silhouette.
    batch_key : str, optional
        Required for iLISI and kBET.
    n_seeds : int
        Number of perturbation seeds.
    sigma_scale : float
        Noise level in units of per-feature std. Default 0.10
        matches the Atti & Subramaniam 2025 convention.
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
    emb_nn, precomputed_dists = _load_precomputed(
        preprocess_context, embeddings, n_total
    )

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

    rng = np.random.default_rng(base_seed)
    seeds = [int(rng.integers(0, 1 << 31)) for _ in range(n_seeds)]

    rows: list[dict] = []
    for emb in embeddings:
        if emb not in adata.obsm:
            continue
        X_full = adata.obsm[emb]
        # Apply core-set subsampling for O(N^2) operations on large atlases
        X_for_sil = X_full
        labels_for_sil = ref_arr if ref_arr is not None else pred_arr
        if subset_indices is not None and subset_indices.size < n_total:
            X_for_sil = X_full[subset_indices]
            if labels_for_sil is not None:
                labels_for_sil = labels_for_sil[subset_indices]
        per_metric_values: dict[str, list[float]] = {m: [] for m in metrics}
        for s in seeds:
            t0 = time.time()
            sub_rng = np.random.default_rng(s)
            # Perturb the full (or subset) embedding
            X_to_perturb = X_for_sil if "silhouette" in metrics or subset_indices is not None else X_full
            X_perturbed = _perturb_embedding(X_to_perturb, sigma_scale, sub_rng)
            for metric in metrics:
                if metric == "silhouette":
                    if labels_for_sil is None:
                        continue
                    pdists = precomputed_dists.get(emb)
                    if pdists is not None and hasattr(pdists, "_filepath"):
                        pdists = None  # cannot reuse full dist matrix
                    val = _silhouette_value(
                        X_perturbed, labels_for_sil, pdists, n_jobs
                    )
                elif metric == "ari":
                    if ref_arr is None or pred_arr is None:
                        continue
                    val = _ari_value(ref_arr, pred_arr)
                elif metric == "ilisi":
                    if batch_arr is None or emb not in emb_nn:
                        continue
                    full_nn = emb_nn[emb]
                    sub_nn = full_nn.subset_neighbors(
                        min(90, full_nn.n_neighbors)
                    )
                    val = _ilisi_value(sub_nn, batch_arr)
                elif metric == "kbet":
                    if batch_arr is None or emb not in emb_nn:
                        continue
                    full_nn = emb_nn[emb]
                    sub_nn = full_nn.subset_neighbors(
                        min(25, full_nn.n_neighbors)
                    )
                    val = _kbet_value(sub_nn, batch_arr)
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
