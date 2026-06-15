import numpy as np
from scipy.sparse import issparse
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors


def run(X, labels, n_jobs=-1, verbose=True, max_samples=None, precomputed_dists=None):
    """
    Calculate the Dunn Index for clustering quality evaluation.

    The Dunn Index is the ratio of the minimum inter-cluster distance to
    the maximum intra-cluster distance.  Higher values indicate better,
    more compact and well-separated clusters.

    On large datasets the computation uses spatial indexing (ball-tree
    nearest-neighbour search) for inter-cluster distances and a chunked
    pairwise strategy for intra-cluster diameters — the full O(n²)
    distance matrix is **never** materialised in memory.

    When *precomputed_dists* is provided the full path is skipped and
    intra‑/inter‑cluster distances are computed by slicing the
    precomputed N×N matrix — O(1) per metric call.

    `Dunn Index readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/dunn/>`__

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix.  Sparse matrices are densified automatically.
    labels : array-like of shape (n_samples,)
        Cluster labels.
    n_jobs : int, default=-1
        Number of parallel jobs for distance computations (-1 = all cores).
    verbose : bool, default=True
        Print progress.
    max_samples : int, optional
        **Ignored** — kept for backward compatibility with older callers.
        The function always processes the full dataset.
    precomputed_dists : ndarray or TriangularMatrix, optional
        Precomputed N×N pairwise distance matrix.  When provided the
        raw *X* is not used — distances are sliced by cluster index.

    Returns
    -------
    float
        Dunn Index.  Range [0, ∞); higher is better.
    """
    if precomputed_dists is not None:
        if issparse(precomputed_dists):
            precomputed_dists = precomputed_dists.toarray()
        if hasattr(precomputed_dists, "_get_block"):
            if verbose:
                print("Using precomputed distances for Dunn Index...")
                print("  TriangularMatrix detected — switching to X-based path...")
            precomputed_dists = None
        else:
            if verbose:
                print("Using precomputed distances for Dunn Index...")
            if precomputed_dists.shape[0] != len(labels):
                raise ValueError(
                    f"precomputed_dists and labels must have same length. "
                    f"Got dists: {precomputed_dists.shape[0]}, labels: {len(labels)}"
                )
            labels_arr = np.asarray(labels)
            unique_lbl = np.unique(labels_arr)
            n_cl = len(unique_lbl)
            if n_cl < 2:
                return 0.0
            cluster_idx = {}
            for lbl in unique_lbl:
                cluster_idx[lbl] = np.where(labels_arr == lbl)[0]
            min_inter = np.inf
            for i in range(n_cl):
                li = unique_lbl[i]
                idx_i = cluster_idx[li]
                for j in range(i + 1, n_cl):
                    lj = unique_lbl[j]
                    idx_j = cluster_idx[lj]
                    d = float(precomputed_dists[np.ix_(idx_i, idx_j)].min())
                    if d < min_inter:
                        min_inter = d
            max_intra = 0.0
            for lbl in unique_lbl:
                idx = cluster_idx[lbl]
                if len(idx) < 2:
                    continue
                dk = float(precomputed_dists[np.ix_(idx, idx)].max())
                if dk > max_intra:
                    max_intra = dk
            if max_intra == 0:
                return 0.0
            dunn_val = min_inter / max_intra
            if verbose:
                print(f"  Dunn Index = {dunn_val:.6f}")
            return dunn_val

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    labels = np.asarray(labels)

    if len(X) != len(labels):
        raise ValueError(
            f"X and labels must have same length. "
            f"Got X: {len(X)}, labels: {len(labels)}"
        )

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    if n_clusters < 2:
        return 0.0

    # ── Per-cluster index arrays ──────────────────────────────────
    cluster_idx = {}
    for lbl in unique_labels:
        cluster_idx[lbl] = np.where(labels == lbl)[0]

    # ── Inter-cluster: minimum distance between clusters ──────────
    # Pre-fit one ball-tree per cluster (22 fits) and query across
    # all C²/2 = 231 pairs.  Each fit is O(|Ck|·log|Ck|·d), each query
    # is O(|Ci|·log|Cj|·d).  Total ~30 s for 85k cells / 50 dims.
    inter_pairs = (n_clusters * (n_clusters - 1)) // 2
    if verbose:
        print(f"Computing Dunn Index " f"({len(X):,} samples, {n_clusters} clusters)...")
        print(f"  Inter-cluster distances ({inter_pairs} pairs)...")

    trees = {}
    for lbl in unique_labels:
        trees[lbl] = NearestNeighbors(n_neighbors=1, algorithm="auto").fit(
            X[cluster_idx[lbl]]
        )

    min_inter = np.inf
    for i in range(n_clusters):
        li = unique_labels[i]
        Xi = X[cluster_idx[li]]
        for j in range(i + 1, n_clusters):
            lj = unique_labels[j]
            Xj = X[cluster_idx[lj]]

            if len(Xi) <= len(Xj):
                dists, _ = trees[lj].kneighbors(Xi)
            else:
                dists, _ = trees[li].kneighbors(Xj)

            d = float(dists.min())
            if d < min_inter:
                min_inter = d

    # ── Intra-cluster: maximum diameter per cluster ───────────────
    # Chunked pairwise to avoid storing the full |Ck|×|Ck| matrix.
    if verbose:
        print(f"  Intra-cluster diameters ({n_clusters} clusters)...")

    max_intra = 0.0

    # Threshold below which we compute the full intra-cluster matrix
    # in one shot (fast, bounded memory); above it we fall back to
    # chunked upper-triangular blocks.
    _DIRECT_THRESHOLD = 10000

    for lbl in unique_labels:
        idx = cluster_idx[lbl]
        nk = len(idx)
        if nk < 2:
            continue

        Xk = X[idx]

        if nk <= _DIRECT_THRESHOLD:
            dk = pairwise_distances(Xk, n_jobs=n_jobs).max()
        else:
            dk = _chunked_diameter(Xk, n_jobs=n_jobs)

        if dk > max_intra:
            max_intra = dk

    if max_intra == 0:
        return 0.0

    dunn_index = min_inter / max_intra

    if verbose:
        print(f"  Dunn Index = {dunn_index:.6f}")

    return dunn_index


# ── Chunked upper-triangular pairwise max ────────────────────────────


def _chunked_diameter(X, n_jobs=1, chunk_size=2000):
    """Maximum pairwise distance in *X*, computed in memory-bounded
    upper-triangular blocks so the full |X|×|X| matrix is never stored.

    Uses GPU-accelerated cdist when JAX is available (~10× faster).
    """
    from .._jax_utils import _JAX_AVAILABLE
    from .._jax_utils import cdist as gpu_cdist

    n = len(X)
    max_d = 0.0
    _use_gpu = _JAX_AVAILABLE

    for i in range(0, n, chunk_size):
        end_i = min(i + chunk_size, n)
        Xi = X[i:end_i]

        for j in range(i, n, chunk_size):
            end_j = min(j + chunk_size, n)
            Xj = X[j:end_j]

            if _use_gpu:
                block = np.array(gpu_cdist(Xi, Xj), copy=True)  # JAX returns read-only
            else:
                block = pairwise_distances(Xi, Xj, n_jobs=n_jobs)
            bmax = float(block.max())
            if bmax > max_d:
                max_d = bmax

    return max_d
