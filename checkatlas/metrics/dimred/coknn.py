import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import logging
from joblib import Parallel, delayed

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True):
    """
    Computes the coKNN metric (Jaccard Index of neighborhoods).
    
    Unlike Entourage (which calculates Recall), coKNN calculates the Jaccard similarity
    between the High-Dim and Low-Dim neighborhoods. It is a stricter metric that
    measures exact set similarity.

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
        coKNN Score (0 to 1).
        Average Jaccard similarity of neighborhoods.
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
    # Query k+1 because index 0 is the cell itself
    query_k = k_neighbors + 1
    
    if verbose: print(f"Finding {k_neighbors} nearest neighbors in High-Dim (PCA)...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs, algorithm='auto').fit(high_dim_data)
    _, indices_high = nbrs_high.kneighbors(high_dim_data)

    if verbose: print(f"Finding {k_neighbors} nearest neighbors in Low-Dim (UMAP)...")
    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs, algorithm='auto').fit(low_dim_data)
    _, indices_low = nbrs_low.kneighbors(low_dim_data)

    # 4. Calculate Jaccard Similarity (coKNN)
    indices_high = indices_high[:, 1:]
    indices_low = indices_low[:, 1:]

    jaccard_sum = 0
    
    iterator = range(n_cells)
    if verbose:
        iterator = tqdm(iterator, desc="Calculating coKNN (Jaccard)")

    def _jaccard(high_row, low_row, k):
        set_high, set_low = set(high_row), set(low_row)
        intersection = len(set_high & set_low)
        union = 2 * k - intersection
        return intersection / union if union > 0 else 0
    
    jaccard_scores = Parallel(n_jobs=n_jobs)(
        delayed(_jaccard)(indices_high[i], indices_low[i], k_neighbors)
        for i in tqdm(range(n_cells), desc="Calculating coKNN", disable=not verbose)
    )
    coknn_score = np.mean(jaccard_scores)

    # if verbose:
    #     print(f"coKNN Score: {coknn_score:.4f}")
    #     logger.debug(f"coKNN : {coknn_score}")

    return coknn_score