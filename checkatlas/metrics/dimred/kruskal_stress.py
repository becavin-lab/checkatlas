import numpy as np
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm
import gc

def run(adata, low_dim_key='X_umap', high_dim_key='X', n_samples=None, seed=42,
        n_jobs=-1, verbose=True, batch_size=1000, sample_pairs=1000000,
        precomputed_high_dists=None, precomputed_low_dists=None):
    """
    Kruskal Stress (Memory-Optimized)
    Measures the mismatch between pairwise distances in high-dimensional and low-dimensional spaces.
    It quantifies how well the global distance structure is preserved.

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
                           Set to None to use all pairs.
        precomputed_high_dists (np.ndarray or memmap): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray or memmap): Precomputed low-dim distance matrix.

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

    # 2. Calculate stress using memory-efficient approach
    n_pairs_total = n_cells * (n_cells - 1) // 2
    
    if sample_pairs is not None and sample_pairs < n_pairs_total:
        # Sampling-based approach
        if verbose: 
            print(f"Using sampling-based stress calculation (sampling {sample_pairs:,} of {n_pairs_total:,} pairs)...")
        stress = _sampled_kruskal_stress(high_dim_dists, low_dim_dists, sample_pairs, seed, verbose)
    else:
        # Full computation
        if verbose: 
            print(f"Computing full Kruskal stress ({n_pairs_total:,} pairs)...")
        stress = _chunked_kruskal_stress(high_dim_dists, low_dim_dists, batch_size, verbose)
    
    return stress


def _sampled_kruskal_stress(high_dists, low_dists, n_pairs, seed, verbose=True):
    """Compute Kruskal stress by sampling distance pairs."""
    n = high_dists.shape[0]
    
    np.random.seed(seed)
    
    # Generate random upper-triangle indices
    if verbose: print("Generating random pair indices...")
    
    sampled_pairs = set()
    while len(sampled_pairs) < n_pairs:
        remaining = n_pairs - len(sampled_pairs)
        rows = np.random.randint(0, n-1, size=int(remaining * 1.2))
        cols = np.random.randint(0, n, size=int(remaining * 1.2))
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
    
    # Extract in batches for memmap efficiency
    H = np.empty(len(sampled_pairs), dtype=np.float32)
    L = np.empty(len(sampled_pairs), dtype=np.float32)
    
    extract_batch_size = 100000
    for i in range(0, len(sampled_pairs), extract_batch_size):
        end = min(i + extract_batch_size, len(sampled_pairs))
        batch_rows = rows[i:end]
        batch_cols = cols[i:end]
        H[i:end] = high_dists[batch_rows, batch_cols]
        L[i:end] = low_dists[batch_rows, batch_cols]
    
    # Normalize
    if verbose: print("Normalizing distances...")
    H_max = np.max(H)
    L_max = np.max(L)
    
    if H_max > 0: H_norm = H / H_max
    else: H_norm = H
        
    if L_max > 0: L_norm = L / L_max
    else: L_norm = L

    # Calculate Stress
    if verbose: print("Calculating stress score...")
    numerator = np.sum(np.square(H_norm - L_norm))
    denominator = np.sum(np.square(L_norm))
    
    # Cleanup
    del H, L, H_norm, L_norm, sampled_pairs, rows, cols
    gc.collect()
    
    if denominator == 0:
        return 0.0
        
    stress = np.sqrt(numerator / denominator)
    return stress


def _chunked_kruskal_stress(high_dists, low_dists, batch_size, verbose=True):
    """Compute Kruskal stress by processing upper triangle in chunks."""
    n = high_dists.shape[0]
    
    # First pass: find max values for normalization
    if verbose: print("Finding max distances for normalization...")
    H_max = 0.0
    L_max = 0.0
    
    for i in tqdm(range(n-1), desc="Finding max", disable=not verbose):
        H_row = high_dists[i, i+1:]
        L_row = low_dists[i, i+1:]
        H_max = max(H_max, np.max(H_row))
        L_max = max(L_max, np.max(L_row))
    
    # Second pass: compute stress
    if verbose: print("Computing stress values...")
    numerator = 0.0
    denominator = 0.0
    
    for i in tqdm(range(n-1), desc="Computing stress", disable=not verbose):
        H_row = high_dists[i, i+1:]
        L_row = low_dists[i, i+1:]
        
        # Normalize
        if H_max > 0: H_norm = H_row / H_max
        else: H_norm = H_row
            
        if L_max > 0: L_norm = L_row / L_max
        else: L_norm = L_row
        
        numerator += np.sum(np.square(H_norm - L_norm))
        denominator += np.sum(np.square(L_norm))
    
    if denominator == 0:
        return 0.0
        
    stress = np.sqrt(numerator / denominator)
    return stress


def _compute_distances_chunked(data, n_jobs=-1, batch_size=1000, desc="Computing distances", verbose=True):
    """Compute pairwise distances in chunks."""
    n_samples = data.shape[0]
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    for i in tqdm(range(0, n_samples, batch_size), desc=desc, disable=not verbose):
        end = min(i + batch_size, n_samples)
        distances[i:end, :] = pairwise_distances(data[i:end], data, n_jobs=n_jobs)
    
    np.fill_diagonal(distances, 0)
    return distances