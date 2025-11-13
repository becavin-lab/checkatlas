import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse


def run(X, labels):
    """
    Calculate the Dunn Index for clustering quality evaluation.
    
    The Dunn Index is an internal clustering validation metric that measures
    the ratio of the minimum inter-cluster distance to the maximum intra-cluster
    distance. Higher values indicate better clustering (well-separated, compact clusters).
    
    `Dunn Index readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/dunn/>`__

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points
        Can be sparse or dense matrix
    :param labels: array-like of shape (n_samples,)
        Cluster labels for each sample
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
        # Cannot compute Dunn Index with less than 2 clusters
        return 0.0
    
    if len(X) != len(labels):
        raise ValueError(f"X and labels must have the same length. "
                        f"Got X: {len(X)}, labels: {len(labels)}")
    
    # Compute pairwise distances between all points
    distances = pairwise_distances(X)
    
    # Calculate minimum inter-cluster distance
    min_inter_cluster_dist = np.inf
    
    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            # Get indices of points in cluster i and j
            cluster_i_indices = np.where(labels == unique_labels[i])[0]
            cluster_j_indices = np.where(labels == unique_labels[j])[0]
            
            # Get distances between points in cluster i and cluster j
            inter_cluster_dists = distances[np.ix_(cluster_i_indices, cluster_j_indices)]
            
            # Find minimum distance between these two clusters
            min_dist = np.min(inter_cluster_dists)
            min_inter_cluster_dist = min(min_inter_cluster_dist, min_dist)
    
    # Calculate maximum intra-cluster distance (diameter)
    max_intra_cluster_dist = 0.0
    
    for label in unique_labels:
        # Get indices of points in this cluster
        cluster_indices = np.where(labels == label)[0]
        
        # Skip if cluster has only one point
        if len(cluster_indices) < 2:
            continue
        
        # Get distances within this cluster
        intra_cluster_dists = distances[np.ix_(cluster_indices, cluster_indices)]
        
        # Find maximum distance within this cluster (diameter)
        max_dist = np.max(intra_cluster_dists)
        max_intra_cluster_dist = max(max_intra_cluster_dist, max_dist)
    
    # Avoid division by zero
    if max_intra_cluster_dist == 0:
        return 0.0
    
    # Calculate Dunn Index
    dunn_index = min_inter_cluster_dist / max_intra_cluster_dist
    
    return dunn_index
