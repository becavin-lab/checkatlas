import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
import warnings


def run(X, labels, n_jobs=-1, verbose=True, max_samples=5000, k=10):
    """
    Calculate the Density-Based Clustering Validation Index (DBCVI).
    
    DBCVI evaluates clustering quality based on density. It defines density
    using the "mutual reachability distance" and constructs a Minimum Spanning Tree (MST)
    to measure density within clusters (DSC) and separation between clusters (DSPC).
    
    Higher DBCVI values indicate better clustering.
    Range: [-1, 1]
    
    Reference:
    Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and Sander, J., 2014.
    Density-based clustering validation.
    In Proceedings of the 2014 SIAM International Conference on Data Mining (pp. 839-847).
    
    :param X: array-like of shape (n_samples, n_features)
        Feature matrix.
    :param labels: array-like of shape (n_samples,)
        Cluster labels.
    :param n_jobs: int, optional (default=-1)
        Number of parallel jobs for kNN computation. -1 uses all processors.
    :param verbose: bool, optional (default=True)
        Whether to print progress.
    :param max_samples: int, optional (default=5000)
        Maximum number of samples for full DBCVI computation. Subsamples if exceeded.
    :param k: int, optional (default=10)
        Number of neighbors for core distance computation.
    :return: float
        DBCVI score.
    """
    X = np.asarray(X)
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    # Filter out noise points if any (often labeled as -1)
    mask = labels != -1
    if np.sum(mask) == 0:
        return 0.0
        
    X_clean = X[mask]
    labels_clean = labels[mask]
    unique_labels = np.unique(labels_clean)
    n_clusters = len(unique_labels)
    
    if n_clusters < 2:
        return 0.0
        
    try:
        from sklearn.neighbors import NearestNeighbors
        
        # Subsample if too large (MST is O(N²))
        if len(X_clean) > max_samples:
            if verbose:
                print(f"Subsampling {max_samples} points for DBCVI (from {len(X_clean)})...")
            np.random.seed(42)
            indices = np.random.choice(len(X_clean), max_samples, replace=False)
            X_clean = X_clean[indices]
            labels_clean = labels_clean[indices]
            unique_labels = np.unique(labels_clean)
            n_clusters = len(unique_labels)
            if n_clusters < 2:
                return 0.0
        
        # Calculate core distances using kNN (parallelized)
        k_actual = min(k, len(X_clean) - 1)
        if k_actual < 1:
            return 0.0
        
        if verbose: print(f"Computing core distances (k={k_actual})...")
        nbrs = NearestNeighbors(n_neighbors=k_actual + 1, n_jobs=n_jobs).fit(X_clean)
        distances, _ = nbrs.kneighbors(X_clean)
        core_dists = distances[:, -1]  # k-th nearest neighbor distance
            
        # Compute pairwise Euclidean distances
        if verbose: print("Computing pairwise distances...")
        dists = squareform(pdist(X_clean))
        
        # Mutual Reachability Distance: max(core_dist(a), core_dist(b), dist(a,b))
        n = len(X_clean)
        core_dists_matrix_1 = np.tile(core_dists, (n, 1))
        core_dists_matrix_2 = core_dists_matrix_1.T
        mrd = np.maximum(dists, np.maximum(core_dists_matrix_1, core_dists_matrix_2))
        
        # Calculate DSC (Density Sparseness of a Cluster) for each cluster
        if verbose: print("Computing cluster density (DSC via MST)...")
        dsc = {}
        cluster_indices = {}
        
        for label in unique_labels:
            idx = np.where(labels_clean == label)[0]
            cluster_indices[label] = idx
            
            if len(idx) < 2:
                dsc[label] = np.max(mrd) if np.max(mrd) > 0 else 1.0
                continue
                
            # MST of cluster subgraph
            sub_mrd = mrd[np.ix_(idx, idx)]
            mst = minimum_spanning_tree(sub_mrd)
            
            # DSC = max weight in the cluster's MST
            if mst.nnz > 0:
                dsc[label] = np.max(mst.data)
            else:
                dsc[label] = 0.0

        # Calculate DSPC (Density Separation between Pairs of Clusters)
        if verbose: print("Computing cluster separation (DSPC)...")
        dspc = {}
        for i in range(n_clusters):
            l1 = unique_labels[i]
            for j in range(i + 1, n_clusters):
                l2 = unique_labels[j]
                
                idx1 = cluster_indices[l1]
                idx2 = cluster_indices[l2]
                
                # Min MRD between the two clusters
                sub_mrd = mrd[np.ix_(idx1, idx2)]
                sep = np.min(sub_mrd)
                dspc[(l1, l2)] = sep
                dspc[(l2, l1)] = sep
        
        # Calculate Validity Index per cluster
        if verbose: print("Computing per-cluster validity indices...")
        validity_indices = []
        weights = []
        
        for label in unique_labels:
            min_sep = np.inf
            for other_label in unique_labels:
                if label == other_label:
                    continue
                sep = dspc.get((label, other_label), np.inf)
                if sep < min_sep:
                    min_sep = sep
            
            if min_sep == np.inf:
                min_sep = 0.0
            
            dens_sparseness = dsc[label]
            
            if max(min_sep, dens_sparseness) == 0:
                vc = 0.0
            else:
                vc = (min_sep - dens_sparseness) / max(min_sep, dens_sparseness)
            
            validity_indices.append(vc)
            weights.append(len(cluster_indices[label]))
            
        # Weighted average
        dbcvi = np.average(validity_indices, weights=weights)
        
        if verbose: print(f"DBCVI: {dbcvi:.4f}")
        return dbcvi

    except Exception as e:
        print(f"Error calculating DBCVI: {e}")
        return np.nan
