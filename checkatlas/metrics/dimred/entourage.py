import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import logging

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)  # Ensure debug prints show up if needed

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True):
    """
    Computes the Entourage metric (also known as Neighborhood Preservation or KNN Overlap).
    
    It measures what fraction of a cell's nearest neighbors in the high-dimensional space (PCA)
    remain among its nearest neighbors in the low-dimensional embedding (UMAP).

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm for reference (default: 'X_pca').
    k_neighbors : int, optional (default=30)
        The size of the local neighborhood. 
        Note: 30 is standard for scRNA-seq (matches Scanpy/Seurat defaults).
    n_samples : int or None
        Number of cells to subsample for speed.
    seed : int
        Random seed.
    n_jobs : int
        Number of parallel jobs for neighbor search (-1 = all cores).
    verbose : bool
        Show progress bars.

    Returns:
    --------
    float
        Entourage Score (0 to 1).
        1.0 = Perfect preservation of local neighborhoods.
        ~0.0 = Neighbors are completely scrambled.
    """

    # 1. Check keys / Compute if missing
    if low_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    if high_dim_key not in adata.obsm.keys():
        if verbose: print(f"Key '{high_dim_key}' not found. Calculating PCA...")
        sc.tl.pca(adata, random_state=seed) 

    # 2. Subsampling
    # While KNN is faster than pairwise distances, subsampling is still useful for 
    # immediate feedback on massive datasets.
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

    n_cells = high_dim_data.shape[0]

    # 3. Fit Nearest Neighbors
    # We query k + 1 because the closest neighbor is always the point itself (dist=0).
    # We will slice that one off later.
    query_k = k_neighbors + 1
    
    if verbose: print(f"Finding {k_neighbors} nearest neighbors in High-Dim (PCA)...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs, algorithm='auto').fit(high_dim_data)
    _, indices_high = nbrs_high.kneighbors(high_dim_data)

    if verbose: print(f"Finding {k_neighbors} nearest neighbors in Low-Dim (UMAP)...")
    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs, algorithm='auto').fit(low_dim_data)
    _, indices_low = nbrs_low.kneighbors(low_dim_data)

    # 4. Calculate Overlap
    # Slice [:, 1:] to remove the self-match (index 0)
    indices_high = indices_high[:, 1:]
    indices_low = indices_low[:, 1:]

    total_intersection = 0
    
    # We iterate through cells to calculate intersection size
    # Python sets are optimized for intersection operations
    iterator = range(n_cells)

    # for i in iterator:
    #     # Convert row to set for O(1) lookups
    #     set_high = set(indices_high[i])
    #     set_low = set(indices_low[i])
        
    #     # Calculate intersection size
    #     intersection_size = len(set_high.intersection(set_low))
    #     total_intersection += intersection_size

    # Vectorized overlap calculation (replaces the for loop)
    # Convert to sets row-wise and compute intersection sizes
    # Using NumPy broadcasting for faster computation

    if verbose:
        print("Calculate final entourage score...")
    
    intersection_counts = np.array([
        len(np.intersect1d(indices_high[i], indices_low[i])) 
        for i in range(n_cells)
    ])
    total_intersection = intersection_counts.sum()

    # 5. Compute Final Score
    # Formula: (Sum of overlaps) / (Total number of neighbors checked)
    # Total neighbors checked = Number of cells * k per cell
    entourage_score = total_intersection / (n_cells * k_neighbors)


    return entourage_score