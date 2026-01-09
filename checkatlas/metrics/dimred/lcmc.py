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
    Computes the LCMC (Local Continuity Meta-Criterion).
    
    LCMC measures the overlap of local neighborhoods (like Entourage), 
    but it is normalized against the baseline probability of overlapping by chance.
    
    Formula: LCMC = (Overlap - Baseline) / (1 - Baseline)
    Where Baseline approx = k / (N - 1)

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding.
    high_dim_key : str
        Key in adata.obsm for reference.
    k_neighbors : int, optional (default=30)
        The size of the local neighborhood.
    n_samples : int or None
        Number of cells to subsample.
    seed : int
        Random seed.
    n_jobs : int
        Parallel jobs.
    verbose : bool
        Show progress.

    Returns:
    --------
    float
        LCMC Score (usually 0 to 1).
        1.0 = Perfect.
        0.0 = Random performance.
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
    # Query k+1 to account for self-neighbor
    query_k = k_neighbors + 1
    
    if verbose: print(f"Finding {k_neighbors} neighbors (High-Dim)...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
    _, indices_high = nbrs_high.kneighbors(high_dim_data)

    if verbose: print(f"Finding {k_neighbors} neighbors (Low-Dim)...")
    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
    _, indices_low = nbrs_low.kneighbors(low_dim_data)

    # 4. Calculate Raw Overlap with parallel processing
    indices_high = indices_high[:, 1:]
    indices_low = indices_low[:, 1:]

    def _overlap(high_row, low_row):
        return len(set(high_row) & set(low_row))

    overlaps = Parallel(n_jobs=n_jobs)(
        delayed(_overlap)(indices_high[i], indices_low[i])
        for i in tqdm(range(n_cells), desc="Calculating Overlap", disable=not verbose)
    )
    total_overlap = sum(overlaps)

    # Average overlap fraction (This is the Entourage Score)
    # Average number of shared neighbors / k
    mean_overlap = (total_overlap / n_cells) / k_neighbors

    # 5. Apply LCMC Adjustment
    # Baseline probability of picking a neighbor by chance
    # Note: The pool of potential neighbors is N-1 (excluding self)
    baseline = k_neighbors / (n_cells - 1)

    # LCMC Formula
    lcmc_score = (mean_overlap - baseline) / (1 - baseline)

    return lcmc_score