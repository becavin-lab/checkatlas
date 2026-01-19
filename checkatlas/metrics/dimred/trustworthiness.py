import numpy as np
import scanpy as sc
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors

def run(adata, low_dim_key='X_umap', high_dim_key='X', k_neighbors=30, 
        n_samples=5000, seed=42, verbose=True, n_jobs=-1,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Computes the Trustworthiness metric.
    Measures preservation of local neighborhood (penalizes false neighbors).
    Compares High-Dim (adata.X) vs Low-Dim (X_umap).
    
    Note: sklearn.manifold.trustworthiness does NOT accept precomputed distance matrices
    or n_jobs. It computes distances internally. Subsampling is the primary
    optimization here.
    """

    # 1. Check keys
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    # 2. Subsampling (Crucial for performance)
    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        if verbose: print(f"Subsampling {n_samples} cells...")
        np.random.seed(seed)
        indices = np.random.choice(n_obs, n_samples, replace=False)
    else:
        indices = np.arange(n_obs)

    # 3. Prepare Data
    if high_dim_key == 'X':
        high_dim_data = adata.X[indices]
        if hasattr(high_dim_data, "toarray"): high_dim_data = high_dim_data.toarray()
    else:
        if high_dim_key not in adata.obsm.keys():
            if verbose: print(f"Calculating {high_dim_key}...")
            sc.tl.pca(adata, random_state=seed)
        high_dim_data = adata.obsm[high_dim_key][indices]

    low_dim_data = adata.obsm[low_dim_key][indices]

    # 4. Calculate Trustworthiness
    # Only pairwise euclidean is supported by simple generic function
    if verbose: print(f"Calculating Trustworthiness (k={k_neighbors})...")
    
    # We could implement a parallel version using joblib and exact formula 
    # if we wanted to use precomputed distances, but sklearn's optimized C implementation
    # is usually fast enough on subsampled data (5000 samples).
    # Re-implementing the rank/penalty logic in python might be slower even with parallelism.
    
    score = trustworthiness(
        high_dim_data, 
        low_dim_data, 
        n_neighbors=k_neighbors, 
        metric='euclidean'
    )

    return score
