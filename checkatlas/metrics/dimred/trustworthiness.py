import numpy as np
import scanpy as sc
from sklearn.manifold import trustworthiness
import logging

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', k_neighbors=30, 
        n_samples=5000, seed=42, verbose=True):
    """
    Computes the Trustworthiness metric.
    
    This metric measures the extent to which the low-dimensional embedding preserves 
    the local neighborhood of the high-dimensional space. Specifically, it penalizes 
    "false neighbors" (points that are close in UMAP but far in PCA).

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm for reference (default: 'X_pca').
    k_neighbors : int, optional (default=30)
        Number of neighbors to check.
    n_samples : int or None
        Number of cells to subsample. Computing rank-based trustworthiness 
        on full datasets is O(N^2) or O(N log N) and very slow. 
        5000 is usually a good trade-off.
    seed : int
        Random seed.
    verbose : bool
        Print status updates.

    Returns:
    --------
    float
        Trustworthiness score (0 to 1).
        1.0 = Highly trustworthy (no false neighbors).
    """

    # 1. Check keys / Compute if missing
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    if high_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{high_dim_key}' not found. Calculating PCA...")
        sc.tl.pca(adata, random_state=seed) 

    # 2. Subsampling (Highly Recommended)
    # sklearn.manifold.trustworthiness calculates pairwise distances internally.
    # For 27k cells, this can be slow. Subsampling to 5k makes it instant.
    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        if verbose:
            print(f"Subsampling {n_samples} cells out of {n_obs} for efficiency...")
        np.random.seed(seed)
        indices = np.random.choice(n_obs, n_samples, replace=False)
        high_dim_data = adata.obsm[high_dim_key][indices]
        low_dim_data = adata.obsm[low_dim_key][indices]
    else:
        if verbose:
            print(f"Using all {n_obs} cells (Warning: May be slow)...")
        high_dim_data = adata.obsm[high_dim_key]
        low_dim_data = adata.obsm[low_dim_key]

    # 3. Calculate Trustworthiness
    # sklearn signature: trustworthiness(X, X_embedded, n_neighbors=5, metric='euclidean')
    # Note: X is high-dim (Ground Truth), X_embedded is low-dim.
    
    if verbose:
        print(f"Calculating Trustworthiness (k={k_neighbors})...")

    score = trustworthiness(
        high_dim_data, 
        low_dim_data, 
        n_neighbors=k_neighbors, 
        metric='euclidean'
    )

    return score
