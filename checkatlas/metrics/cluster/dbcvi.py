import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
import warnings

def run(X, labels):
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
        
    # Calculate Mutual Reachability Distance (MRD)
    # We need k-nearest neighbor distance for "core distance"
    # For simplicity and speed without hdbscan, we can approximate or use a fixed k.
    # Original paper suggests k. Let's use k=min(10, n_samples-1).
    
    try:
        from sklearn.neighbors import NearestNeighbors
        k = min(10, len(X_clean) - 1)
        if k < 1:
            return 0.0
            
        nbrs = NearestNeighbors(n_neighbors=k+1).fit(X_clean)
        distances, _ = nbrs.kneighbors(X_clean)
        core_dists = distances[:, -1] # k-th nearest neighbor distance
        
        # Calculate all pairwise Euclidean distances
        # For large datasets, this is expensive (O(N^2)). 
        # We might need to subsample if N is very large.
        if len(X_clean) > 5000:
            warnings.warn("Dataset too large for full DBCVI calculation. Subsampling to 5000 points.")
            indices = np.random.choice(len(X_clean), 5000, replace=False)
            X_clean = X_clean[indices]
            labels_clean = labels_clean[indices]
            core_dists = core_dists[indices]
            
        dists = squareform(pdist(X_clean))
        
        # Mutual Reachability Distance: max(core_dist(a), core_dist(b), dist(a,b))
        # Efficient calculation using broadcasting
        n = len(X_clean)
        core_dists_matrix_1 = np.tile(core_dists, (n, 1))
        core_dists_matrix_2 = core_dists_matrix_1.T
        mrd = np.maximum(dists, np.maximum(core_dists_matrix_1, core_dists_matrix_2))
        
        # Calculate DSC (Density Sparseness of a Cluster) for each cluster
        dsc = {}
        cluster_indices = {}
        
        for label in unique_labels:
            indices = np.where(labels_clean == label)[0]
            cluster_indices[label] = indices
            
            if len(indices) < 2:
                dsc[label] = 0.0 # Or max possible? 0 density means infinite sparseness?
                # Actually DSC is "sparseness", so lower is denser (better).
                # If 1 point, density is undefined or infinite.
                # Let's handle it gracefully.
                dsc[label] = np.max(mrd) if np.max(mrd) > 0 else 1.0
                continue
                
            # Extract subgraph for this cluster
            sub_mrd = mrd[np.ix_(indices, indices)]
            
            # MST of the cluster
            # minimum_spanning_tree expects CSR or dense matrix
            mst = minimum_spanning_tree(sub_mrd)
            
            # DSC is the maximum weight in the MST (internal sparseness)
            # "The density sparseness of a cluster C is defined as the maximum edge weight of its MST"
            # Wait, the paper says: DSC(C) = max weight of edges in MST_MRD(C)
            if mst.nnz > 0:
                dsc[label] = np.max(mst.data)
            else:
                dsc[label] = 0.0

        # Calculate DSPC (Density Separation Pair of Clusters)
        # DSPC(Ci, Cj) = min weight of edge separating Ci and Cj in the global MST?
        # No, "minimum reachability distance between points of Ci and points of Cj"
        # "DSPC(Ci, Cj) = min(MRD(x, y)) for x in Ci, y in Cj"
        # But we also need to consider if they are "connected" via other clusters?
        # The paper defines DSPC as the minimum reachability distance between the two clusters.
        
        dspc = {}
        for i in range(n_clusters):
            l1 = unique_labels[i]
            for j in range(i + 1, n_clusters):
                l2 = unique_labels[j]
                
                idx1 = cluster_indices[l1]
                idx2 = cluster_indices[l2]
                
                # Submatrix of MRD between C1 and C2
                sub_mrd = mrd[np.ix_(idx1, idx2)]
                
                sep = np.min(sub_mrd)
                dspc[(l1, l2)] = sep
                dspc[(l2, l1)] = sep
        
        # Calculate Validity Index for each cluster
        # VC(C) = (min_j(DSPC(C, Cj)) - DSC(C)) / max(min_j(DSPC(C, Cj)), DSC(C))
        validity_indices = []
        weights = []
        
        for label in unique_labels:
            # Find min separation to any other cluster
            min_sep = np.inf
            for other_label in unique_labels:
                if label == other_label:
                    continue
                sep = dspc.get((label, other_label), np.inf)
                if sep < min_sep:
                    min_sep = sep
            
            if min_sep == np.inf:
                min_sep = 0.0 # Should not happen if >1 cluster
            
            dens_sparseness = dsc[label]
            
            if max(min_sep, dens_sparseness) == 0:
                vc = 0.0
            else:
                vc = (min_sep - dens_sparseness) / max(min_sep, dens_sparseness)
            
            validity_indices.append(vc)
            weights.append(len(cluster_indices[label]))
            
        # Weighted average
        dbcvi = np.average(validity_indices, weights=weights)
        
        return dbcvi

    except Exception as e:
        print(f"Error calculating DBCVI: {e}")
        return np.nan
