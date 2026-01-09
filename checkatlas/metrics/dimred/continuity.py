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
    Computes the Continuity metric.
    
    Continuity measures the preservation of original neighbors (Recall).
    It penalizes "missing neighbors" -- points that are close in High-Dim (PCA)
    but end up far apart in Low-Dim (UMAP).
    
    A low Continuity score implies the embedding has "torn" the manifold, 
    breaking continuous trajectories into artificial separate clusters.

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
        Number of cells to subsample for speed.
    seed : int
        Random seed.
    verbose : bool
        Print status updates.

    Returns:
    --------
    float
        Continuity score (0 to 1).
        1.0 = Perfect continuity (manifold is intact).
        < 0.7 = Significant tearing of the data structure.
    """

    # 1. Check keys / Compute if missing
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    if high_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{high_dim_key}' not found. Calculating PCA...")
        sc.tl.pca(adata, random_state=seed) 

    # 2. Subsampling
    # Like Trustworthiness, this is O(N^2) or O(N log N). Subsampling is required.
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

    # 3. Calculate Continuity
    # TRICK: We use the trustworthiness function but SWAP the inputs.
    # We treat Low-Dim as the "Ground Truth" and High-Dim as the "Embedding".
    # This checks: "Are the neighbors in High-Dim preserved in Low-Dim?"
    
    if verbose:
        print(f"Calculating Continuity (Swapped Trustworthiness, k={k_neighbors})...")

    # NOTICE THE SWAP: trustworthiness(Low, High) instead of (High, Low)
    score = trustworthiness(
        low_dim_data,   # <--- Treated as Reference
        high_dim_data,  # <--- Treated as Embedding
        n_neighbors=k_neighbors, 
        metric='euclidean'
    )

    return score