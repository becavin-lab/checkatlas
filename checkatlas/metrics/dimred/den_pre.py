import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr
import logging

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True, log_transform=False):
    """
    Computes the Density Preservation metric (Pearson correlation of local radii).
    
    It assesses whether the heterogeneity (crowdedness) of cell states is preserved.
    It calculates the distance to the k-th nearest neighbor in High-Dim and Low-Dim
    and computes their correlation.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm for reference (default: 'X_pca').
    k_neighbors : int, optional (default=30)
        The neighbor index to define the local radius. 
        We measure the distance to the 30th neighbor.
    n_samples : int or None
        Number of cells to subsample for speed.
    seed : int
        Random seed.
    n_jobs : int
        Parallel jobs for neighbor search.
    verbose : bool
        Show progress.

    Returns:
    --------
    float
        Density Preservation Score (Pearson Correlation, -1 to 1).
        Higher is better (relative density patterns are preserved).
    """

    # 1. Check keys / Compute if missing
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    if high_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{high_dim_key}' not found. Calculating PCA...")
        sc.tl.pca(adata, random_state=seed) 

    # 2. Subsampling
    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        if verbose:
            print(f"Subsampling {n_samples} cells out of {n_obs}...")
        np.random.seed(seed)
        indices = np.random.choice(n_obs, n_samples, replace=False)
        high_dim_data = adata.obsm[high_dim_key][indices]
        low_dim_data = adata.obsm[low_dim_key][indices]
    else:
        if verbose:
            print(f"Using all {n_obs} cells...")
        high_dim_data = adata.obsm[high_dim_key]
        low_dim_data = adata.obsm[low_dim_key]

    # 3. Compute Radii (Distance to k-th neighbor)
    # We query k+1 neighbors because the 1st neighbor is the point itself (dist=0).
    # We want the distance to the last neighbor in the list.
    query_k = k_neighbors + 1
    
    if verbose: print(f"Computing local radii (k={k_neighbors}) in High-Dim...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
    dists_high, _ = nbrs_high.kneighbors(high_dim_data)
    
    # Take the distance to the furthest neighbor in the set (the k-th neighbor)
    # This represents the "radius" of the local neighborhood.
    radii_high = dists_high[:, -1]

    if verbose: print(f"Computing local radii (k={k_neighbors}) in Low-Dim...")
    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
    dists_low, _ = nbrs_low.kneighbors(low_dim_data)
    
    radii_low = dists_low[:, -1]

    if log_transform:
        radii_high = np.log1p(radii_high)  # log(1 + x) avoids log(0)
        radii_low = np.log1p(radii_low)

    # 4. Compute Pearson Correlation
    # We check if larger radii in High-Dim correspond to larger radii in Low-Dim.
    if verbose: print("Calculating correlation of radii...")
    
    # Note: Log-transforming distances can sometimes improve this metric 
    # because volumes grow exponentially, but standard correlation is the baseline.
    correlation, p_value = pearsonr(radii_high, radii_low)

    # if verbose:
    #     print(f"Density Preservation (Pearson): {correlation:.4f}")
    #     logger.debug(f"Density Preservation : {correlation}")

    return correlation