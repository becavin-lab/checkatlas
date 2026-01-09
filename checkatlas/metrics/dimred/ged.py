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
    Computes the Normalized Graph Edit Distance (GED) on k-NN graphs.
    
    It calculates the fraction of edges that do not match between the High-Dim 
    and Low-Dim graphs (Symmetric Difference).
    
    GED = (Edges unique to High + Edges unique to Low) / (Total Edges in both)

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm for embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm for reference (default: 'X_pca').
    k_neighbors : int, optional (default=30)
        The number of edges per node (k).
    n_samples : int or None
        Number of cells to subsample.
    seed : int
        Random seed.
    n_jobs : int
        Parallel jobs for neighbor search.
    verbose : bool
        Show progress.

    Returns:
    --------
    float
        Normalized GED Score (0 to 1).
        0.0 = Graphs are identical (Perfect preservation).
        1.0 = Graphs are completely disjoint.
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
    # Query k+1 (exclude self)
    query_k = k_neighbors + 1
    
    if verbose: print(f"Building High-Dim Graph (k={k_neighbors})...")
    nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
    _, indices_high = nbrs_high.kneighbors(high_dim_data)

    if verbose: print(f"Building Low-Dim Graph (k={k_neighbors})...")
    nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
    _, indices_low = nbrs_low.kneighbors(low_dim_data)

    # 4. Calculate Graph Edit Distance with parallel processing
    # Slice off self-match
    indices_high = indices_high[:, 1:]
    indices_low = indices_low[:, 1:]

    total_edges_possible = n_cells * k_neighbors * 2  # Total edges in Graph A + Graph B

    def _symmetric_diff(high_row, low_row):
        return len(set(high_row).symmetric_difference(set(low_row)))

    sym_diffs = Parallel(n_jobs=n_jobs)(
        delayed(_symmetric_diff)(indices_high[i], indices_low[i])
        for i in tqdm(range(n_cells), desc="Calculating Graph Edits", disable=not verbose)
    )
    total_symmetric_difference = sum(sym_diffs)

    # 5. Normalize
    # We divide by the sum of cardinalities (|E_H| + |E_L|)
    ged_score = total_symmetric_difference / total_edges_possible

    return ged_score

# Example Usage
# ged = run(adata, low_dim_key='X_umap', high_dim_key='X_pca')