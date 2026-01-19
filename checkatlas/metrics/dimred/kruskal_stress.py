import numpy as np
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm
from joblib import Parallel, delayed

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=5000, seed=42,
        n_jobs=-1, verbose=True, batch_size=2000, 
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Kruskal Stress
    Measures the mismatch between pairwise distances in high-dimensional and low-dimensional spaces.
    It quantifies how well the global distance structure is preserved.

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
        float: The Kruskal Stress score.

    Interpretation:
        Range 0 to 1 (typically).
        Lower is better (0 means perfect preservation of distances).
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
        high_dim_dists = _compute_distances_with_progress(high_dim_data, n_jobs, batch_size, "High-Dim Dists", verbose)
        low_dim_dists = _compute_distances_with_progress(low_dim_data, n_jobs, batch_size, "Low-Dim Dists", verbose)

    # 4. Extract upper triangle indices (excluding diagonal)
    if verbose: print("Extracting pairwise distances...")
    triu_indices = np.triu_indices_from(high_dim_dists, k=1)
    H = high_dim_dists[triu_indices]
    L = low_dim_dists[triu_indices]

    # 5. Normalize Distances
    if verbose: print("Normalizing distances...")
    H_max = np.max(H)
    L_max = np.max(L)
    
    # Avoid division by zero
    if H_max > 0: H_norm = H / H_max
    else: H_norm = H
        
    if L_max > 0: L_norm = L / L_max
    else: L_norm = L

    # 6. Calculate Stress
    if verbose: print("Calculating stress score...")
    numerator = np.sum(np.square(H_norm - L_norm))
    denominator = np.sum(np.square(L_norm))
    
    if denominator == 0:
        return 0.0
        
    stress = np.sqrt(numerator / denominator)
    return stress

def _compute_distances_with_progress(data, n_jobs=-1, batch_size=2000, desc="Computing distances", verbose=True):
    """Batch-wise parallel distance computation."""
    n_samples = data.shape[0]
    if n_samples <= batch_size:
        return pairwise_distances(data, n_jobs=n_jobs)
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    n_batches = (n_samples + batch_size - 1) // batch_size
    
    # Simple loop with internal parallelism
    with tqdm(total=n_batches, desc=desc, disable=not verbose) as pbar:
        for i in range(0, n_samples, batch_size):
            end_i = min(i + batch_size, n_samples)
            distances[i:end_i, :] = pairwise_distances(data[i:end_i], data, n_jobs=n_jobs)
            pbar.update(1)
            
    np.fill_diagonal(distances, 0)
    return distances