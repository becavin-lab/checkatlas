import numpy as np
import scanpy as sc
from sklearn.manifold import trustworthiness

def run(adata, low_dim_key='X_umap', high_dim_key='X', k_neighbors=30, 
        n_samples=5000, seed=42, verbose=True, n_jobs=-1,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Computes the Continuity metric.
    Measures preservation of original neighbors (penalizes broken trajectories).
    Compares High-Dim (adata.X) vs Low-Dim (X_umap).
    
    Implemented as Trustworthiness with SWAPPED inputs.
    Note: sklearn.manifold.trustworthiness does NOT accept precomputed distance matrices.
    """

    # 1. Check keys
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    # 2. Subsampling
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

    # 4. Calculate Continuity
    if verbose: print(f"Calculating Continuity (Swapped Trustworthiness, k={k_neighbors})...")
    
    # SWAPPED: We check if Low-Dim neighbors are preserved in High-Dim (Wait, no?)
    # Continuity checks if High-Dim neighbors are preserved in Low-Dim.
    # Trustworthiness: High -> Low (False Neighbors: Low neighbor NOT High neighbor)
    # Continuity: Low -> High (Missing neighbor: High neighbor NOT Low neighbor)
    # So we treat Low as Reference and High as Embedding? 
    # Yes, sklearn documentation says trustworthiness(X, X_embedded) penalizes rank errors.
    # To check Continuity, we penalize rank errors in High-Dim relative to Low-Dim??
    # Actually the standard trick is to swap them.
    
    score = trustworthiness(
        low_dim_data,   # <--- Treated as Reference (The space we define neighborhood in?)
        high_dim_data,  # <--- Treated as Embedding
        n_neighbors=k_neighbors, 
        metric='euclidean'
    )

    return score