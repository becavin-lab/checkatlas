import numpy as np
import scipy.stats
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm
from joblib import Parallel, delayed

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=5000, seed=42,
        n_jobs=-1, verbose=True, batch_size=2000,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Computes Spearman's Rho using precomputed distances if available.
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

    # 4. Extract upper triangle indices (Flattening)
    if verbose: print("Extracting pairwise distance vectors...")
    
    triu_indices = np.triu_indices_from(high_dim_dists, k=1)
    H_vec = high_dim_dists[triu_indices]
    L_vec = low_dim_dists[triu_indices]

    # 5. Compute Spearman's Rho
    if verbose: print("Computing Spearman correlation (Ranking)...")
    
    rho_result = scipy.stats.spearmanr(H_vec, L_vec)
    return rho_result.statistic


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