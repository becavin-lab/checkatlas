import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
from joblib import Parallel, delayed
import logging

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True):
    """
    Computes the Average Jaccard Distance.
    
    This measures the dissimilarity between the k-NN graph in High-Dim vs Low-Dim.
    It is equivalent to (1 - coKNN).
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding.
    high_dim_key : str
        Key in adata.obsm for reference.
    k_neighbors : int, optional (default=30)
        Size of the neighborhood.
    n_samples : int or None
        Number of cells to subsample.
    seed : int
        Random seed.
    verbose : bool
        Show progress.

    Returns:
    --------
    float
        Average Jaccard Distance (0 to 1).
        0.0 = Perfect match (Neighborhoods are identical).
        1.0 = Complete mismatch.
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

    n_cells = high_dim_data.shape[0]

    # 3. Fit Nearest Neighbors
    # Query k+1 (to exclude self)
    query_k = k_neighbors + 1
    
    if verbose: print(f"Finding neighbors (k={k_neighbors})...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
    _, indices_high = nbrs_high.kneighbors(high_dim_data)

    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
    _, indices_low = nbrs_low.kneighbors(low_dim_data)

    # 4. Calculate Jaccard Distance with parallel processing
    # Slice off the first column (self-match)
    indices_high = indices_high[:, 1:]
    indices_low = indices_low[:, 1:]

    def _jaccard_distance(high_row, low_row, k):
        set_high, set_low = set(high_row), set(low_row)
        intersection = len(set_high & set_low)
        union = 2 * k - intersection
        j_index = intersection / union if union > 0 else 0
        return 1.0 - j_index

    jaccard_distances = Parallel(n_jobs=n_jobs)(
        delayed(_jaccard_distance)(indices_high[i], indices_low[i], k_neighbors)
        for i in tqdm(range(n_cells), desc="Calculating Jaccard Dist", disable=not verbose)
    )

    # 5. Average
    avg_jaccard_dist = np.mean(jaccard_distances)

    return avg_jaccard_dist

# Example Usage
# j_dist = run(adata, low_dim_key='X_umap', high_dim_key='X_pca')