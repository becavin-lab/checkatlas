import numpy as np
import scanpy as sc
from sklearn.metrics import pairwise_distances
from tqdm import tqdm
import gc

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=None, seed=42, 
        n_jobs=-1, verbose=True, batch_size=1000, sample_pairs=1000000,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Distance Correlation (dCor) - Memory-Optimized
    Measures the dependence between the pairwise distance matrices of the high-dimensional and low-dimensional data.
    It assesses how well the global structure (distances) is preserved.

    Uses a sampling-based approach for large datasets to avoid memory issues.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        n_samples (int): Number of samples to subsample for calculation. None = use all.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs.
        verbose (bool): Whether to print progress.
        batch_size (int): Batch size for chunked processing.
        sample_pairs (int): Number of cells to sample for dCor (smaller subset for double-centering).
                           dCor requires O(n^2) memory for centering, so we sample cells, not pairs.
        precomputed_high_dists (np.ndarray or memmap): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray or memmap): Precomputed low-dim distance matrix.

    Returns:
        float: The Distance Correlation score.

    Interpretation:
        Range 0 to 1.
        Higher is better (1 means perfect correlation of distances, 0 means no correlation).
    """

    # 1. Use Precomputed Distances if available
    if precomputed_high_dists is not None and precomputed_low_dists is not None:
        if verbose: print("Using precomputed distance matrices...")
        high_dim_dists = precomputed_high_dists
        low_dim_dists = precomputed_low_dists
        n_cells = high_dim_dists.shape[0]
    else:
        # Fallback to local computation
        if verbose: print("Precomputed distances not provided. Calculating locally...")
        
        # Check keys
        if low_dim_key not in adata.obsm.keys():
            if verbose: print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        # Prepare Data
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

        # Compute Distances
        high_dim_dists = _compute_distances_chunked(high_dim_data, n_jobs, batch_size, "High-Dim Dists", verbose)
        low_dim_dists = _compute_distances_chunked(low_dim_data, n_jobs, batch_size, "Low-Dim Dists", verbose)

    # 2. For dCor, we need to subsample CELLS (not pairs) because double-centering requires full submatrix
    # dCor on 5000 cells needs ~200MB for centering. Limit to sample_pairs (interpreted as n_cells here).
    max_cells_for_dcor = min(5000, sample_pairs // 1000) if sample_pairs else 5000
    
    if n_cells > max_cells_for_dcor:
        if verbose: 
            print(f"Subsampling {max_cells_for_dcor} cells for dCor (from {n_cells})...")
        np.random.seed(seed)
        sample_idx = np.random.choice(n_cells, max_cells_for_dcor, replace=False)
        
        # Extract submatrix
        A = np.zeros((max_cells_for_dcor, max_cells_for_dcor), dtype=np.float32)
        B = np.zeros((max_cells_for_dcor, max_cells_for_dcor), dtype=np.float32)
        
        for i in tqdm(range(max_cells_for_dcor), desc="Extracting submatrix", disable=not verbose):
            for j in range(max_cells_for_dcor):
                A[i, j] = high_dim_dists[sample_idx[i], sample_idx[j]]
                B[i, j] = low_dim_dists[sample_idx[i], sample_idx[j]]
    else:
        # Use full matrices (copy to ensure contiguous for speed)
        A = np.array(high_dim_dists, dtype=np.float32)
        B = np.array(low_dim_dists, dtype=np.float32)
    
    # 3. Compute Distance Correlation
    if verbose: print("Calculating Distance Correlation (Double Centering)...")
    dcor_score = _distance_correlation(A, B)
    
    # Cleanup
    del A, B
    gc.collect()
    
    return dcor_score


def _distance_correlation(A, B):
    """Computes the Distance Correlation between two distance matrices A and B."""
    # 1. Double Centering (Optimized with numpy broadcasting)
    A_centered = _double_center(A)
    B_centered = _double_center(B)

    # 2. Compute components
    dCov_sq = np.mean(A_centered * B_centered)
    dVar_A_sq = np.mean(A_centered * A_centered)
    dVar_B_sq = np.mean(B_centered * B_centered)

    if dVar_A_sq > 0 and dVar_B_sq > 0:
        dCor = np.sqrt(dCov_sq / (np.sqrt(dVar_A_sq) * np.sqrt(dVar_B_sq)))
    else:
        dCor = 0.0
    return dCor


def _double_center(D):
    """Performs double centering efficiently."""
    row_means = np.mean(D, axis=1, keepdims=True)
    col_means = np.mean(D, axis=0, keepdims=True)
    grand_mean = np.mean(D)
    return D - row_means - col_means + grand_mean


def _compute_distances_chunked(data, n_jobs=-1, batch_size=1000, desc="Computing distances", verbose=True):
    """Compute pairwise distances in chunks."""
    n_samples = data.shape[0]
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    for i in tqdm(range(0, n_samples, batch_size), desc=desc, disable=not verbose):
        end = min(i + batch_size, n_samples)
        distances[i:end, :] = pairwise_distances(data[i:end], data, n_jobs=n_jobs)
    
    np.fill_diagonal(distances, 0)
    return distances