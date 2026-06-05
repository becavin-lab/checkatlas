import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree


def run(X, labels, n_jobs=-1, verbose=True, max_samples=None, k=10):
    """
    Calculate the Density-Based Clustering Validation Index (DBCVI).

    DBCVI evaluates clustering quality based on density.  It builds a
    mutual-reachability-distance graph, constructs a Minimum Spanning
    Tree per cluster, and compares intra-cluster density (DSC) against
    inter-cluster separation (DSPC).

    **Large‑data strategy**: when *N* > 10 000 the full O(N²) distance
    matrix is replaced by an approximate kNN graph (``pynndescent``)
    that is both memory‑friendly and fast while preserving the DBCVI
    semantics — the MST on the kNN‑MRD graph is equivalent to the MST
    on the full MRD graph for all edges that matter.

    Higher DBCVI values indicate better clustering.
    Range: [-1, 1]

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix.-----------
    labels : array-like of shape (n_samples,)
        Cluster labels.
    n_jobs : int, default=-1
        Parallel jobs for kNN computation.
    verbose : bool, default=True
        Print progress.
    max_samples : int, optional
        Maximum samples for exact DBCVI before switching to
        the approximate kNN‑graph path.  ``None`` disables
        subsampling entirely — the kNN path is always used for
        large atlases.
    k : int, default=10
        Number of neighbours for core‑distance computation.

    Returns
    -------
    float
        DBCVI score.
    """
    X = np.asarray(X, dtype=np.float64)
    labels = np.asarray(labels)

    mask = labels != -1
    if np.sum(mask) == 0:
        return 0.0

    X_clean = X[mask]
    labels_clean = labels[mask]
    unique_labels = np.unique(labels_clean)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    n = len(X_clean)
    k_core = min(k, n - 1)
    if k_core < 1:
        return 0.0

    # ── path selection ────────────────────────────────────────────
    use_knn_path = n > 10000
    if max_samples is not None and n > max_samples:
        use_knn_path = True

    try:
        if use_knn_path:
            dbcvi = _dbcvi_knn_approx(
                X_clean,
                labels_clean,
                unique_labels,
                n_clusters,
                k_core,
                n_jobs,
                verbose,
            )
        else:
            dbcvi = _dbcvi_exact(
                X_clean,
                labels_clean,
                unique_labels,
                n_clusters,
                k_core,
                n_jobs,
                verbose,
            )
        if verbose:
            print(f"DBCVI: {dbcvi:.4f}")
        return dbcvi
    except Exception as e:
        print(f"Error calculating DBCVI: {e}")
        return np.nan


# ═══════════════════════════════════════════════════════════════════════
# Exact path (N ≤ 10 000) — full O(N²) distance matrix
# ═══════════════════════════════════════════════════════════════════════


def _dbcvi_exact(X, labels, unique_labels, n_clusters, k_core, n_jobs, verbose):
    from scipy.spatial.distance import pdist, squareform
    from sklearn.neighbors import NearestNeighbors

    if verbose:
        print(f"Computing core distances (k={k_core})...")
    nbrs = NearestNeighbors(n_neighbors=k_core + 1, n_jobs=n_jobs).fit(X)
    distances, _ = nbrs.kneighbors(X)
    core_dists = distances[:, -1]

    if verbose:
        print("Computing full pairwise distances...")
    dists = squareform(pdist(X))

    n = len(X)
    cdm = np.tile(core_dists, (n, 1))
    mrd = np.maximum(dists, np.maximum(cdm, cdm.T))

    return _compute_dbcvi(mrd, labels, unique_labels, n_clusters, verbose)


# ═══════════════════════════════════════════════════════════════════════
# Approximate path (N > 10 000) — kNN‑graph based, sparse MRD
# ═══════════════════════════════════════════════════════════════════════


def _dbcvi_knn_approx(X, labels, unique_labels, n_clusters, k_core, n_jobs, verbose):
    n = len(X)
    # kNN for MRD graph — use ~100 neighbours to capture cluster borders
    knn_k = min(100, n - 1)

    if verbose:
        print(f"Building kNN graph (k={knn_k}) for approximate DBCVI...")

    try:
        from pynndescent import NNDescent

        index = NNDescent(
            X,
            n_neighbors=knn_k + 1,
            metric="euclidean",
            n_jobs=n_jobs,
            random_state=42,
        )
        knn_idx, knn_dist = index.neighbor_graph
    except ImportError:
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=knn_k + 1, n_jobs=n_jobs).fit(X)
        knn_dist, knn_idx = nbrs.kneighbors(X)

    # Drop self-neighbour (column 0)
    knn_idx = knn_idx[:, 1:]
    knn_dist = knn_dist[:, 1:]

    # Compute core distances from the same kNN graph (using k_core)
    core_dists = knn_dist[:, k_core - 1].copy()

    # ── Build sparse MRD matrix (only on kNN edges) ───────────────
    rows = np.repeat(np.arange(n), knn_k)
    cols = knn_idx.ravel()
    dists = knn_dist.ravel()

    # MRD(ij) = max(core_dist_i, core_dist_j, dist_ij)
    mrd_vals = np.maximum(dists, np.maximum(core_dists[rows], core_dists[cols]))

    # Symmetrise: also add edge (j, i) for each (i, j)
    rows_sym = np.concatenate([rows, cols])
    cols_sym = np.concatenate([cols, rows])
    vals_sym = np.concatenate([mrd_vals, mrd_vals])

    mrd_sparse = csr_matrix((vals_sym, (rows_sym, cols_sym)), shape=(n, n))

    if verbose:
        nnz = mrd_sparse.nnz
        print(f"  Sparse MRD: {nnz:,} edges ({nnz/n**2*100:.2f}% of full)")

    return _compute_dbcvi_sparse(mrd_sparse, labels, unique_labels, n_clusters, verbose)


# ═══════════════════════════════════════════════════════════════════════
# Shared DBCVI computation (works on dense or sparse MRD)
# ═══════════════════════════════════════════════════════════════════════


def _compute_dbcvi(mrd, labels, unique_labels, n_clusters, verbose):
    """DBCVI from a dense MRD matrix."""
    cluster_indices = {}
    for lbl in unique_labels:
        cluster_indices[lbl] = np.where(labels == lbl)[0]

    dsc = _dsc_dense(mrd, unique_labels, cluster_indices)
    dspc = _dspc_dense(mrd, unique_labels, cluster_indices, n_clusters)
    return _validity_index(unique_labels, cluster_indices, dsc, dspc)


def _compute_dbcvi_sparse(mrd_sparse, labels, unique_labels, n_clusters, verbose):
    """DBCVI from a sparse MRD matrix."""
    cluster_indices = {}
    for lbl in unique_labels:
        cluster_indices[lbl] = np.where(labels == lbl)[0]

    dsc = _dsc_sparse(mrd_sparse, unique_labels, cluster_indices)
    dspc = _dspc_sparse(mrd_sparse, unique_labels, cluster_indices, n_clusters)
    return _validity_index(unique_labels, cluster_indices, dsc, dspc)


# ── DSC (Density Sparseness of a Cluster) ──────────────────────────


def _dsc_dense(mrd, unique_labels, cluster_indices):
    """DSC = max MST edge weight per cluster (dense)."""
    dsc = {}
    for lbl in unique_labels:
        idx = cluster_indices[lbl]
        if len(idx) < 2:
            dsc[lbl] = np.max(mrd) if np.max(mrd) > 0 else 1.0
            continue
        sub = mrd[np.ix_(idx, idx)]
        mst = minimum_spanning_tree(sub)
        dsc[lbl] = float(np.max(mst.data)) if mst.nnz > 0 else 0.0
    return dsc


def _dsc_sparse(mrd_sparse, unique_labels, cluster_indices):
    """DSC = max MST edge weight per cluster (sparse sub‑extraction)."""
    dsc = {}
    global_max = mrd_sparse.data.max() if mrd_sparse.nnz > 0 else 1.0
    for lbl in unique_labels:
        idx = cluster_indices[lbl]
        if len(idx) < 2:
            dsc[lbl] = global_max
            continue
        sub = mrd_sparse[idx, :][:, idx]
        mst = minimum_spanning_tree(sub)
        dsc[lbl] = float(np.max(mst.data)) if mst.nnz > 0 else 0.0
    return dsc


# ── DSPC (Density Separation between Pairs of Clusters) ─────────────


def _dspc_dense(mrd, unique_labels, cluster_indices, n_clusters):
    """Min MRD between each cluster pair (dense)."""
    dspc = {}
    for i in range(n_clusters):
        l1 = unique_labels[i]
        idx1 = cluster_indices[l1]
        for j in range(i + 1, n_clusters):
            l2 = unique_labels[j]
            idx2 = cluster_indices[l2]
            sep = float(np.min(mrd[np.ix_(idx1, idx2)]))
            dspc[(l1, l2)] = sep
            dspc[(l2, l1)] = sep
    return dspc


def _dspc_sparse(mrd_sparse, unique_labels, cluster_indices, n_clusters):
    """Min MRD between each cluster pair (sparse)."""
    dspc = {}
    for i in range(n_clusters):
        l1 = unique_labels[i]
        idx1 = cluster_indices[l1]
        for j in range(i + 1, n_clusters):
            l2 = unique_labels[j]
            idx2 = cluster_indices[l2]
            sub = mrd_sparse[idx1, :][:, idx2]
            sep = float(np.min(sub.data)) if sub.nnz > 0 else np.inf
            dspc[(l1, l2)] = sep
            dspc[(l2, l1)] = sep
    return dspc


# ── Validity index (shared) ─────────────────────────────────────────


def _validity_index(unique_labels, cluster_indices, dsc, dspc):
    validity_indices = []
    weights = []

    for lbl in unique_labels:
        min_sep = np.inf
        for other in unique_labels:
            if lbl == other:
                continue
            sep = dspc.get((lbl, other), np.inf)
            if sep < min_sep:
                min_sep = sep

        if min_sep == np.inf:
            min_sep = 0.0

        dens = dsc[lbl]
        denom = max(min_sep, dens)
        vc = (min_sep - dens) / denom if denom > 0 else 0.0
        validity_indices.append(vc)
        weights.append(len(cluster_indices[lbl]))

    return float(np.average(validity_indices, weights=weights))
