import numpy as np
import scipy.stats
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm
from joblib import Parallel, delayed
import gc

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=None, seed=42,
        n_jobs=-1, verbose=True, batch_size=1000, sample_pairs=1000000,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Spearman's Rho (Memory-Optimized)
    Measures the rank correlation between pairwise distances in high-dimensional and low-dimensional spaces.
    It assesses how well the relative ordering of distances is preserved.

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
        sample_pairs (int): Number of distance pairs to sample for large datasets (default: 1M).
                           Set to None to use all pairs (may cause memory issues for large N).
        precomputed_high_dists (np.ndarray or memmap): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray or memmap): Precomputed low-dim distance matrix.

    Returns:
        float: The Spearman's Rho score.

    Interpretation:
        Range -1 to 1.
        Higher is better (1 means perfect monotonic relationship of distances).
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

        # Compute Distances in chunks
        high_dim_dists = _compute_distances_chunked(high_dim_data, n_jobs, batch_size, "High-Dim Dists", verbose)
        low_dim_dists = _compute_distances_chunked(low_dim_data, n_jobs, batch_size, "Low-Dim Dists", verbose)

    # 2. Calculate Spearman correlation using memory-efficient approach
    n_pairs_total = n_cells * (n_cells - 1) // 2  # Upper triangle size
    
    if sample_pairs is not None and sample_pairs < n_pairs_total:
        # Sampling-based approach for large datasets
        if verbose: 
            print(f"Using sampling-based Spearman (sampling {sample_pairs:,} of {n_pairs_total:,} pairs)...")
        rho = _sampled_spearman(high_dim_dists, low_dim_dists, sample_pairs, seed, verbose)
    else:
        # Full computation for smaller datasets
        if verbose: 
            print(f"Computing full Spearman correlation ({n_pairs_total:,} pairs)...")
        rho = _chunked_spearman(high_dim_dists, low_dim_dists, batch_size, verbose)
    
    return rho


def _sampled_spearman(high_dists, low_dists, n_pairs, seed, verbose=True):
    """Compute Spearman correlation by sampling distance pairs."""
    n = high_dists.shape[0]
    
    np.random.seed(seed)
    
    # Generate random upper-triangle indices
    if verbose: print("Generating random pair indices...")
    
    # Efficient sampling: generate row, col pairs where row < col
    sampled_pairs = set()
    while len(sampled_pairs) < n_pairs:
        remaining = n_pairs - len(sampled_pairs)
        # Over-sample to account for duplicates
        rows = np.random.randint(0, n-1, size=int(remaining * 1.2))
        cols = np.random.randint(0, n, size=int(remaining * 1.2))
        # Ensure row < col
        mask = rows < cols
        rows, cols = rows[mask], cols[mask]
        for r, c in zip(rows, cols):
            if len(sampled_pairs) >= n_pairs:
                break
            sampled_pairs.add((r, c))
    
    # Extract sampled distances
    if verbose: print("Extracting sampled distance pairs...")
    sampled_pairs = list(sampled_pairs)
    rows = np.array([p[0] for p in sampled_pairs])
    cols = np.array([p[1] for p in sampled_pairs])
    
    # Extract in batches to handle memmap arrays efficiently
    H_vec = np.empty(len(sampled_pairs), dtype=np.float32)
    L_vec = np.empty(len(sampled_pairs), dtype=np.float32)
    
    batch_size = 100000
    for i in range(0, len(sampled_pairs), batch_size):
        end = min(i + batch_size, len(sampled_pairs))
        batch_rows = rows[i:end]
        batch_cols = cols[i:end]
        H_vec[i:end] = high_dists[batch_rows, batch_cols]
        L_vec[i:end] = low_dists[batch_rows, batch_cols]
    
    # Compute Spearman correlation
    if verbose: print("Computing Spearman correlation...")
    rho_result = scipy.stats.spearmanr(H_vec, L_vec)
    
    # Cleanup
    del H_vec, L_vec, sampled_pairs, rows, cols
    gc.collect()
    
    return rho_result.statistic


def _chunked_spearman(high_dists, low_dists, batch_size, verbose=True):
    """Compute Spearman correlation by extracting upper triangle in chunks."""
    n = high_dists.shape[0]
    
    # Collect upper triangle values row by row
    H_list = []
    L_list = []
    
    for i in tqdm(range(n-1), desc="Extracting distances", disable=not verbose):
        # For row i, upper triangle is columns i+1 to n
        H_list.append(high_dists[i, i+1:])
        L_list.append(low_dists[i, i+1:])
        
        # Periodically concatenate to avoid too many small arrays
        if len(H_list) >= batch_size:
            H_list = [np.concatenate(H_list)]
            L_list = [np.concatenate(L_list)]
    
    H_vec = np.concatenate(H_list)
    L_vec = np.concatenate(L_list)
    
    del H_list, L_list
    gc.collect()
    
    if verbose: print("Computing Spearman correlation...")
    rho_result = scipy.stats.spearmanr(H_vec, L_vec)
    
    del H_vec, L_vec
    gc.collect()
    
    return rho_result.statistic


def _compute_distances_chunked(data, n_jobs=-1, batch_size=1000, desc="Computing distances", verbose=True):
    """Compute pairwise distances in chunks."""
    n_samples = data.shape[0]
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    for i in tqdm(range(0, n_samples, batch_size), desc=desc, disable=not verbose):
        end = min(i + batch_size, n_samples)
        distances[i:end, :] = pairwise_distances(data[i:end], data, n_jobs=n_jobs)
    
    np.fill_diagonal(distances, 0)
    return distances