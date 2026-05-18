import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse


def run(X, labels, n_jobs=-1, verbose=True, max_samples=None):
    """
    Calculate the Dunn Index for clustering quality evaluation.

    The Dunn Index is the ratio of the minimum inter-cluster distance to
    the maximum intra-cluster distance.  Higher values indicate better,
    more compact and well-separated clusters.

    On large datasets the computation uses spatial indexing (ball-tree
    nearest-neighbour search) for inter-cluster distances and a chunked
    pairwise strategy for intra-cluster diameters — the full O(n²)
    distance matrix is **never** materialised in memory.

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

    Returns
    -------
    float
        Dunn Index.  Range [0, ∞); higher is better.
    """
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
    # Uses NearestNeighbors (ball-tree for low-dim data) so each
    # cross-cluster query is O(|Ci|・log|Cj|) instead of O(|Ci|・|Cj|).
    inter_pairs = (n_clusters * (n_clusters - 1)) // 2
    if verbose:
        print(
            f"Computing Dunn Index "
            f"({len(X):,} samples, {n_clusters} clusters)..."
        )
        print(f"  Inter-cluster distances ({inter_pairs} pairs)...")

    min_inter = np.inf
    for i in range(n_clusters):
        li = unique_labels[i]
        Xi = X[cluster_idx[li]]
        for j in range(i + 1, n_clusters):
            lj = unique_labels[j]
            Xj = X[cluster_idx[lj]]

            # Fit ball-tree on one cluster, query the other for the
            # single nearest-neighbour distance.  For low‑dimensional
            # data ball‑tree is lightning fast → single-threaded.
            if len(Xi) <= len(Xj):
                nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
                nn.fit(Xj)
                dists, _ = nn.kneighbors(Xi)
            else:
                nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
                nn.fit(Xi)
                dists, _ = nn.kneighbors(Xj)

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
    upper-triangular blocks so the full |X|×|X| matrix is never stored."""
    n = len(X)
    max_d = 0.0

    for i in range(0, n, chunk_size):
        end_i = min(i + chunk_size, n)
        Xi = X[i:end_i]

        for j in range(i, n, chunk_size):
            end_j = min(j + chunk_size, n)
            Xj = X[j:end_j]

            block = pairwise_distances(Xi, Xj, n_jobs=n_jobs)
            bmax = float(block.max())
            if bmax > max_d:
                max_d = bmax

    return max_d
