import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
from joblib import Parallel, delayed

def run(adata, low_dim_key='X_umap', high_dim_key='X', k_neighbors=30, 
        n_samples=5000, seed=42, n_jobs=-1, verbose=True,
        precomputed_high_knn=None, precomputed_low_knn=None):
    """
    Graph Edit Distance (GED)
    Measures the normalized graph edit distance between the k-NN graphs of the high-dimensional and low-dimensional data.
    It counts the number of edge insertions and deletions needed to transform one graph to the other.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs.
        verbose (bool): Whether to print progress.
        precomputed_high_knn (np.ndarray): Precomputed k-NN indices for high-dim data.
        precomputed_low_knn (np.ndarray): Precomputed k-NN indices for low-dim data.

    Returns:
        float: The Normalized Graph Edit Distance.

    Interpretation:
        Range 0 to 1.
        Lower is better (0 means identical graphs, 1 means completely different graphs).
    """

    # 1. Use Precomputed Neighbors if available
    if precomputed_high_knn is not None and precomputed_low_knn is not None:
        if verbose: print("Using precomputed k-NN graphs...")
        indices_high = precomputed_high_knn
        indices_low = precomputed_low_knn
        n_cells = indices_high.shape[0]
        
        # Precomputed (from metrics.py) includes self. Slice it.
        indices_high = indices_high[:, 1:]
        indices_low = indices_low[:, 1:]
    else:
        # Fallback to local computation
        if verbose: print("Precomputed k-NN not provided. Calculating locally...")

        # Check keys
        if low_dim_key not in adata.obsm.keys():
            if verbose: print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        # Prepare Data
        # Subsampling
        n_obs = adata.n_obs
        if n_samples is not None and n_samples < n_obs:
            if verbose: print(f"Subsampling {n_samples} cells...")
            np.random.seed(seed)
            indices = np.random.choice(n_obs, n_samples, replace=False)
        else:
            indices = np.arange(n_obs)

        # High Dim Data
        if high_dim_key == 'X':
            high_dim_data = adata.X[indices]
            if hasattr(high_dim_data, "toarray"): high_dim_data = high_dim_data.toarray()
        else:
            if high_dim_key not in adata.obsm.keys():
                if verbose: print(f"Calculating {high_dim_key}...")
                sc.tl.pca(adata, random_state=seed)
            high_dim_data = adata.obsm[high_dim_key][indices]

        low_dim_data = adata.obsm[low_dim_key][indices]
        n_cells = high_dim_data.shape[0]

        # Fit Nearest Neighbors
        query_k = k_neighbors + 1
        
        if verbose: print(f"Building High-Dim Graph (k={k_neighbors})...")
        nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(high_dim_data)
        _, indices_high = nbrs_high.kneighbors(high_dim_data)

        if verbose: print(f"Building Low-Dim Graph (k={k_neighbors})...")
        nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(low_dim_data)
        _, indices_low = nbrs_low.kneighbors(low_dim_data)

        # Slice off self-match
        indices_high = indices_high[:, 1:]
        indices_low = indices_low[:, 1:]

    # 4. Calculate Graph Edit Distance (Parallel)
    total_edges_possible = n_cells * k_neighbors * 2  # Total edges in Graph A + Graph B

    def _symmetric_diff(high_row, low_row):
        return len(set(high_row).symmetric_difference(set(low_row)))

    sym_diffs = Parallel(n_jobs=n_jobs)(
        delayed(_symmetric_diff)(indices_high[i], indices_low[i])
        for i in tqdm(range(n_cells), desc="Calculating Graph Edits", disable=not verbose)
    )
    total_symmetric_difference = sum(sym_diffs)

    # 5. Normalize
    ged_score = total_symmetric_difference / total_edges_possible
    return ged_score