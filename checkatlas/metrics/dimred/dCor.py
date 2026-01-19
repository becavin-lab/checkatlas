import numpy as np
import scanpy as sc
from sklearn.metrics import pairwise_distances
from tqdm import tqdm
from joblib import Parallel, delayed

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=5000, seed=42, 
        n_jobs=-1, verbose=True, batch_size=2000,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Distance Correlation (dCor)
    Measures the dependence between the pairwise distance matrices of the high-dimensional and low-dimensional data.
    It assesses how well the global structure (distances) is preserved.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        n_samples (int): Number of samples to subsample for calculation.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs.
        verbose (bool): Whether to print progress.
        batch_size (int): Batch size for distance computation.
        precomputed_high_dists (np.ndarray): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray): Precomputed low-dim distance matrix.

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
    else:
        # Fallback to local computation
        if verbose: print("Precomputed distances not provided. Calculating locally...")
        
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

        # Compute Distances
        high_dim_dists = _compute_distances_parallel(high_dim_data, n_jobs, batch_size, "High-Dim Dists", verbose)
        low_dim_dists = _compute_distances_parallel(low_dim_data, n_jobs, batch_size, "Low-Dim Dists", verbose)

    # 4. Compute Distance Correlation
    if verbose: print("Calculating Distance Correlation (Double Centering)...")
    dcor_score = _distance_correlation(high_dim_dists, low_dim_dists)
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


def _compute_batch_distances(data, start_idx, end_idx, n_jobs):
    """Compute distances for a single batch."""
    return pairwise_distances(data[start_idx:end_idx], data, n_jobs=n_jobs)


def _compute_distances_parallel(data, n_jobs=-1, batch_size=2000, desc="Computing distances", verbose=True):
    """Batch-wise parallel distance computation."""
    n_samples = data.shape[0]
    if n_samples <= batch_size:
        return pairwise_distances(data, n_jobs=n_jobs)
    
    batch_indices = [(i, min(i + batch_size, n_samples)) for i in range(0, n_samples, batch_size)]
    
    batch_results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_compute_batch_distances)(data, start, end, 1)
        for start, end in tqdm(batch_indices, desc=desc, disable=not verbose)
    )
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    for idx, (start, end) in enumerate(batch_indices):
        distances[start:end, :] = batch_results[idx]
            
    np.fill_diagonal(distances, 0)
    return distances