import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr

def run(adata, low_dim_key='X_umap', high_dim_key='X', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True, log_transform=False,
        precomputed_high_knn_dists=None, precomputed_low_knn_dists=None,
        precomputed_high_knn=None, precomputed_low_knn=None): # Added for API consistency
    """
    Density Preservation
    Measures the preservation of local density by calculating the Pearson correlation between the radii of k-nearest neighbors in high-dimensional and low-dimensional spaces.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs.
        verbose (bool): Whether to print progress.
        log_transform (bool): Whether to log-transform the radii before correlation.
        precomputed_high_knn_dists (np.ndarray): Precomputed k-NN distances for high-dim data.
        precomputed_low_knn_dists (np.ndarray): Precomputed k-NN distances for low-dim data.
        precomputed_high_knn (np.ndarray): Precomputed k-NN indices (not used directly, but kept for API consistency).
        precomputed_low_knn (np.ndarray): Precomputed k-NN indices (not used directly, but kept for API consistency).

    Returns:
        float: The Density Preservation score (Pearson Correlation).

    Interpretation:
        Range -1 to 1.
        Higher is better (1 means perfect preservation of relative density).
    """

    # 1. Use Precomputed Neighbors if available
    if precomputed_high_knn_dists is not None and precomputed_low_knn_dists is not None:
        if verbose: print("Using precomputed k-NN graphs...")
        dists_high = precomputed_high_knn_dists
        dists_low = precomputed_low_knn_dists
    else:
        # Fallback to local computation
        if verbose: print("Precomputed k-NN not provided. Calculating locally...")

        # Check keys
        if low_dim_key not in adata.obsm.keys():
            if verbose: print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        # Prepare Data
        # Subsampling
        n_obs = adata.n_obs
        if n_samples is not None and n_samples < n_obs:
            if verbose: print(f"Subsampling {n_samples} cells...")
            np.random.seed(seed)
            indices = np.random.choice(n_obs, n_samples, replace=False)
        else:
            indices = np.arange(n_obs)

        # High Dim Data
        if high_dim_key == 'X':
            high_dim_data = adata.X[indices]
            if hasattr(high_dim_data, "toarray"): high_dim_data = high_dim_data.toarray()
        else:
            if high_dim_key not in adata.obsm.keys():
                if verbose: print(f"Calculating {high_dim_key}...")
                sc.tl.pca(adata, random_state=seed)
            high_dim_data = adata.obsm[high_dim_key][indices]

        low_dim_data = adata.obsm[low_dim_key][indices]

        # Fit Nearest Neighbors
        # Query k+1 (to exclude self)
        query_k = k_neighbors + 1
        
        if verbose: print(f"Computing local radii (k={k_neighbors}) in High-Dim...")
        nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
        dists_high, _ = nbrs_high.kneighbors(high_dim_data)
        
        if verbose: print(f"Computing local radii (k={k_neighbors}) in Low-Dim...")
        nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
        dists_low, _ = nbrs_low.kneighbors(low_dim_data)

    # 3. Compute Radii
    # Take the distance to the furthest neighbor in the set (the k-th neighbor)
    radii_high = dists_high[:, -1]
    radii_low = dists_low[:, -1]

    if log_transform:
        radii_high = np.log1p(radii_high)
        radii_low = np.log1p(radii_low)

    # 4. Compute Pearson Correlation
    if verbose: print("Calculating correlation of radii...")
    correlation, p_value = pearsonr(radii_high, radii_low)
    
    return correlation