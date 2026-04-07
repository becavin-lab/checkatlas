import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse


def run(X, labels, n_jobs=-1, verbose=True):
    """
    Calculate the Dunn Index for clustering quality evaluation.
    
    The Dunn Index is an internal clustering validation metric that measures
    the ratio of the minimum inter-cluster distance to the maximum intra-cluster
    distance. Higher values indicate better clustering (well-separated, compact clusters).
    
    `Dunn Index readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/dunn/>`__

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points. Can be sparse or dense.
    :param labels: array-like of shape (n_samples,)
        Cluster labels for each sample
    :param n_jobs: int, default=-1
        Number of parallel jobs for pairwise distance computation.
        -1 uses all available cores.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return: float
        The Dunn Index score. Higher values indicate better clustering.
        Range: [0, infinity), where higher is better.
    
    Notes
    -----
    The Dunn Index is defined as:
    
    .. math::
        DI = \\frac{\\min_{i \\neq j} \\delta(C_i, C_j)}{\\max_k \\Delta(C_k)}
    
    where:
    - δ(C_i, C_j) is the inter-cluster distance (minimum distance between clusters)
    - Δ(C_k) is the intra-cluster distance (maximum diameter of a cluster)
    
    """
    # Convert to numpy arrays, handle sparse matrices
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
    
    labels = np.asarray(labels)
    
    # Get unique cluster labels
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    # Handle edge cases
    if n_clusters < 2:
        return 0.0
    
    if len(X) != len(labels):
        raise ValueError(f"X and labels must have the same length. "
                        f"Got X: {len(X)}, labels: {len(labels)}")
    
    if verbose:
        print(f"Computing Dunn Index ({len(X):,} samples, {n_clusters} clusters)...")
    
    # Compute pairwise distances with parallelism
    if verbose:
        print(f"  Computing pairwise distances (n_jobs={n_jobs})...")
    distances = pairwise_distances(X, n_jobs=n_jobs)
    
    # Pre-compute cluster masks and indices for vectorized operations
    cluster_masks = {}
    cluster_indices = {}
    for label in unique_labels:
        mask = labels == label
        cluster_masks[label] = mask
        cluster_indices[label] = np.where(mask)[0]
    
    # Vectorized inter-cluster minimum distances
    if verbose:
        print(f"  Computing inter-cluster distances ({n_clusters * (n_clusters - 1) // 2} pairs)...")
    
    min_inter_cluster_dist = np.inf
    for i in range(n_clusters):
        idx_i = cluster_indices[unique_labels[i]]
        for j in range(i + 1, n_clusters):
            idx_j = cluster_indices[unique_labels[j]]
            # Use np.ix_ for efficient submatrix extraction
            min_dist = distances[np.ix_(idx_i, idx_j)].min()
            if min_dist < min_inter_cluster_dist:
                min_inter_cluster_dist = min_dist
    
    # Vectorized intra-cluster maximum distances
    if verbose:
        print(f"  Computing intra-cluster diameters ({n_clusters} clusters)...")
    
    max_intra_cluster_dist = 0.0
    for label in unique_labels:
        idx = cluster_indices[label]
        if len(idx) < 2:
            continue
        max_dist = distances[np.ix_(idx, idx)].max()
        if max_dist > max_intra_cluster_dist:
            max_intra_cluster_dist = max_dist
    
    # Avoid division by zero
    if max_intra_cluster_dist == 0:
        return 0.0
    
    dunn_index = min_inter_cluster_dist / max_intra_cluster_dist
    
    if verbose:
        print(f"  Dunn Index = {dunn_index:.6f}")
    
    return dunn_index
