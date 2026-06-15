from typing import Callable, Optional

import numpy as np
from scipy.sparse import issparse
from sklearn.metrics import pairwise_distances

# ── GPU cdist import (optional, falls back to sklearn) ────────────
_gpu_cdist: Optional[Callable[..., np.ndarray]] = None
_HAS_JAX: bool = False
try:
    from .._jax_utils import _JAX_AVAILABLE as _HAS_JAX
    from .._jax_utils import cdist as _gpu_cdist
except ImportError:
    pass


def run(X, labels, n_jobs=-1, verbose=True, sample_size=None, precomputed_dists=None):
    """
    Calculate the Average Silhouette Width for clustering quality.

    Silhouette measures how well each point fits within its assigned
    cluster relative to its nearest neighbouring cluster.
    Values range from -1 (misassigned) to +1 (well assigned).

    **Large‑data strategy**: intra‑cluster distances *a(i)* are computed
    exactly via bounded‑memory chunked pairwise blocks.  Inter‑cluster
    distances *b(i)* use the squared‑Euclidean centroid‑expansion
    formula — an O(N·K·d) approximation that avoids the full O(N²)
    distance matrix entirely while preserving cluster ordering.

    When *precomputed_dists* is provided the full path is skipped and
    the N×N distance matrix is passed directly to sklearn's
    ``silhouette_score(metric="precomputed")`` for an exact, O(1)
    computation per call.

    `Average Silhouette Width readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/silhouette/>`__

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix.  Sparse inputs are densified.
    labels : array-like of shape (n_samples,)
        Cluster labels.
    n_jobs : int, default=-1
        Parallel jobs (-1 = all cores).
    verbose : bool, default=True
        Print progress.
    sample_size : int, optional
        **Ignored** — always processes every sample.
    precomputed_dists : ndarray or TriangularMatrix, optional
        Precomputed N×N pairwise distance matrix.  When provided the
        raw *X* is not used.

    Returns
    -------
    float
        Mean silhouette coefficient.  Range [-1, 1].
    """
    if precomputed_dists is not None:
        if issparse(precomputed_dists):
            precomputed_dists = precomputed_dists.toarray()
        labels_arr = np.asarray(labels)
        n_total = precomputed_dists.shape[0]
        if n_total != len(labels_arr):
            raise ValueError(
                f"precomputed_dists and labels must have same length. "
                f"Got dists: {n_total}, labels: {len(labels_arr)}"
            )
        if len(np.unique(labels_arr)) < 2:
            return 0.0
        if hasattr(precomputed_dists, "to_dense"):
            # TriangularMatrix — too large to densify (O(N²) memory).
            # Fall through to the X-based centroid-expansion algorithm
            # which is O(N·K·d) and avoids the full distance matrix.
            if verbose:
                print("TriangularMatrix detected — using X-based silhouette...")
        else:
            if verbose:
                print("Using precomputed distances for Silhouette...")
            from sklearn.metrics import silhouette_score

            return float(silhouette_score(precomputed_dists, labels_arr, metric="precomputed"))

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    labels = np.asarray(labels)
    n_total = len(X)

    if n_total != len(labels):
        raise ValueError(
            f"X and labels must have same length. "
            f"Got X: {n_total}, labels: {len(labels)}"
        )

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    if n_clusters < 2:
        return 0.0

    if verbose:
        print(
            f"Computing Average Silhouette Width "
            f"({n_total:,} samples, {n_clusters} clusters, n_jobs={n_jobs})..."
        )

    # ── Step 1: a(i) — exact intra-cluster mean Euclidean distance ─
    a_values = np.zeros(n_total)
    cluster_sizes = {}

    for lbl in unique_labels:
        idx = np.where(labels == lbl)[0]
        nk = len(idx)
        cluster_sizes[lbl] = nk
        if nk < 2:
            continue
        Xk = X[idx]
        if verbose and nk > 5000:
            print(f"  Intra-cluster '{lbl}' ({nk:,} samples)…")

        _exact_intra_mean(Xk, a_values, idx, n_jobs=n_jobs, chunk=2000)

    # ── Step 2: b(i) via centroid-expansion formula ───────────────
    # For point i and cluster C:
    #   mean_sq_dist(i, C) = ||i||² + mean(||c||²) - 2·i·centroid
    #   b_approx(i) = min_{C≠cl(i)} sqrt(mean_sq_dist(i, C))
    #
    # Pre-compute per-cluster statistics.
    norms_sq = np.einsum("ij,ij->i", X, X)  # ||x||² per point

    cluster_centroids = {}
    cluster_mean_norms = {}
    for lbl in unique_labels:
        idx = np.where(labels == lbl)[0]
        cluster_centroids[lbl] = X[idx].mean(axis=0)
        cluster_mean_norms[lbl] = norms_sq[idx].mean()

    b_values = np.full(n_total, np.inf)

    # For each cluster, compute b(i) against every *other* cluster in
    # one vectorised block (O(|Ck|·K·d)).
    for lbl_i in unique_labels:
        idx_i = np.where(labels == lbl_i)[0]
        Xi = X[idx_i]
        ni = len(idx_i)

        # All data for points in this cluster
        norms_i = norms_sq[idx_i]  # shape (ni,)

        # For each other cluster, compute mean_sq_dist block
        for lbl_j in unique_labels:
            if lbl_i == lbl_j:
                continue
            centroid_j = cluster_centroids[lbl_j]  # shape (d,)
            mean_norm_j = cluster_mean_norms[lbl_j]

            # mean_sq_dist(i, Cj) = ||i||² + mean(||c||²) - 2·i·centroid_j
            # shape: (ni,)
            mean_sq = norms_i + mean_norm_j - 2.0 * np.dot(Xi, centroid_j)
            # Clamp near-zero negatives (floating-point noise)
            mean_sq = np.maximum(mean_sq, 0.0)
            mean_dist = np.sqrt(mean_sq)

            b_values[idx_i] = np.minimum(b_values[idx_i], mean_dist)

    if verbose:
        print(f"  Completed b(i) approximation via centroid expansion.")

    # ── Step 3: per-point silhouette ──────────────────────────────
    s = _silhouette_from_apoint(a_values, b_values)

    if verbose:
        print(f"  ASW = {s:.6f}")

    return s


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════


def _exact_intra_mean(Xk, a_out, global_idx, n_jobs=1, chunk=2000):
    """
    Compute a(i) = mean Euclidean distance from each point to all
    *other* points in the same cluster.  Uses chunked upper-triangular
    pairwise blocks so the full |Ck|×|Ck| matrix is never stored.

    When JAX is available, GPU cdist replaces sklearn pairwise_distances
    for ~100× faster intra-cluster distance computation.
    """
    n = len(Xk)
    counts = np.zeros(n)  # per-point count of other points seen
    _use_gpu = _HAS_JAX and _gpu_cdist is not None

    for r in range(0, n, chunk):
        re = min(r + chunk, n)
        Xr = Xk[r:re]
        for c in range(r, n, chunk):
            ce = min(c + chunk, n)
            Xc = Xk[c:ce]
            if _use_gpu:
                block = np.array(_gpu_cdist(Xr, Xc), copy=True)  # JAX returns read-only
            else:
                block = pairwise_distances(Xr, Xc, n_jobs=n_jobs)
            # block shape: (re-r) × (ce-c)

            if r == c:
                # diagonal — zero out self-distances
                dlen = block.shape[0]
                for d in range(dlen):
                    if d < block.shape[1]:
                        block[d, d] = np.nan
                with np.errstate(all="ignore"):
                    row_sum = np.nansum(block, axis=1)
                row_cnt = (ce - c) - 1  # exclude self
            else:
                row_sum = block.sum(axis=1)
                row_cnt = ce - c

            a_out[global_idx[r:re]] += row_sum
            counts[r:re] += row_cnt

    # Normalise: a(i) = sum / (n-1)
    for i in range(n):
        if counts[i] > 0:
            a_out[global_idx[i]] /= counts[i]


def _silhouette_from_apoint(a, b):
    """Per-point silhouette from precomputed a(i) and b(i)."""
    denom = np.maximum(a, b)
    mask = denom > 0
    s = np.zeros(len(a))
    s[mask] = (b[mask] - a[mask]) / denom[mask]
    return float(np.mean(s))


# Keep old name for backward-compat references
_silhouette_per_point = _silhouette_from_apoint
