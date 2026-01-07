import numpy as np
import scipy.stats
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', n_samples=5000, seed=42,
        n_jobs=-1, verbose=True, batch_size=2000):
    """
    Computes Spearman's Rho to assess how well the global topology (relative ordering of cells)
    is preserved in the low-dimensional embedding compared to the high-dimensional space.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm to use as the embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm to use as the reference high-dim space (default: 'X_pca').
    n_samples : int or None
        Number of cells to subsample. Computing ranks for pairwise distances of 
        27k+ cells is extremely memory intensive. 5000 is usually sufficient.
    seed : int
        Random seed for reproducibility.
    n_jobs : int, optional (default=-1)
        Number of parallel jobs for distance computation.
    verbose : bool, optional (default=True)
        Whether to show progress bars.
    batch_size : int, optional (default=2000)
        Batch size for distance computation.

    Returns:
    --------
    float
        Spearman's Rho (-1 to 1). 
        Higher (+1) is better. Indicates excellent preservation of neighborhood hierarchy.
    """
    
    # 1. Check if keys exist and compute if missing
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
            print(f"Subsampling {n_samples} cells out of {n_obs} for efficiency...")
        np.random.seed(seed)
        indices = np.random.choice(n_obs, n_samples, replace=False)
        high_dim_data = adata.obsm[high_dim_key][indices]
        low_dim_data = adata.obsm[low_dim_key][indices]
    else:
        if verbose:
            print(f"Using all {n_obs} cells (Warning: High memory usage for Rank calculation)...")
        high_dim_data = adata.obsm[high_dim_key]
        low_dim_data = adata.obsm[low_dim_key]
    
    # 3. Compute Pairwise Distances
    # We reuse the helper function for parallel processing
    high_dim_dists = _compute_distances_with_progress(
        high_dim_data, 
        n_jobs=n_jobs, 
        batch_size=batch_size,
        desc="Computing high-dim distances (PCA)",
        verbose=verbose
    )
    
    low_dim_dists = _compute_distances_with_progress(
        low_dim_data, 
        n_jobs=n_jobs, 
        batch_size=batch_size,
        desc="Computing low-dim distances (UMAP)",
        verbose=verbose
    )

    # 4. Extract upper triangle indices (Flattening)
    if verbose:
        print("Extracting pairwise distance vectors...")
    
    # We only need the upper triangle (k=1 removes diagonal)
    triu_indices = np.triu_indices_from(high_dim_dists, k=1)
    
    # Flatten into 1D arrays
    # Note: For 5000 cells, this results in vectors of length ~12.5 million.
    # This fits easily in RAM (approx 100MB).
    H_vec = high_dim_dists[triu_indices]
    L_vec = low_dim_dists[triu_indices]

    # 5. Compute Spearman's Rho
    # scipy.stats.spearmanr calculates correlation on ranks
    if verbose:
        print("Computing Spearman correlation (Ranking)...")
    
    # We use scipy.stats.spearmanr. 
    # It returns a generic object with 'statistic' and 'pvalue'. We need 'statistic'.
    rho_result = scipy.stats.spearmanr(H_vec, L_vec)
    
    return rho_result.statistic


def _compute_distances_with_progress(data, n_jobs=-1, batch_size=2000, desc="Computing distances", verbose=True):
    """
    Helper function: Compute pairwise distances with batch-wise progress tracking.
    """
    n_samples = data.shape[0]
    
    if n_samples <= batch_size or not verbose:
        return pairwise_distances(data, n_jobs=n_jobs)
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    # Calculate total batches
    n_batches = int(np.ceil(n_samples / batch_size))
    
    with tqdm(total=n_batches, desc=desc, disable=not verbose) as pbar:
        for i in range(0, n_samples, batch_size):
            end_i = min(i + batch_size, n_samples)
            distances[i:end_i, :] = pairwise_distances(
                data[i:end_i], 
                data, 
                n_jobs=n_jobs
            )
            pbar.update(1)
    
    np.fill_diagonal(distances, 0)
    return distances