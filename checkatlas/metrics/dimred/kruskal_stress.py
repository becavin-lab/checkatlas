import numpy as np
from sklearn.metrics import pairwise_distances
import scanpy as sc
from tqdm import tqdm


def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', n_samples=5000, seed=42,
        n_jobs=-1, verbose=True, batch_size=2000):
    """
    Computes the Kruskal stress to assess how well global distances are preserved 
    in the low-dimensional embedding compared to the high-dimensional space.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm to use as the embedding (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm to use as the reference high-dim space (default: 'X_pca').
    n_samples : int or None
        Number of cells to subsample for calculation. 
        Using all 27k+ cells is computationally expensive (O(N^2)). 
        Set to None to use all cells (warning: high memory usage).
    seed : int
        Random seed for reproducibility.
    n_jobs : int, optional (default=-1)
        Number of parallel jobs for distance computation.
        -1 means using all available CPU cores.
    verbose : bool, optional (default=True)
        Whether to show progress bars and status messages.
    batch_size : int, optional (default=2000)
        Batch size for progress tracking during distance computation.
        Larger batch_size = faster but less frequent progress updates.

    Returns:
    --------
    float
        The Kruskal stress value (0 to 1). Lower is better preservation.
    """
    
    # 1. Check if keys exist
    if low_dim_key not in adata.obsm.keys():
        if verbose:
            print(f"Key '{low_dim_key}' not found in adata.obsm \n Calculating UMAP...")

        ## calculate the "X_umap" if not present
        sc.tl.umap(adata, n_components=2, random_state=seed)

    if high_dim_key not in adata.obsm.keys():
        if verbose:
            print(f"Key '{high_dim_key}' not found in adata.obsm \n Calculating PCA...")
        
        ## calculate the "X_pca" if not present
        sc.tl.pca(adata, random_state=seed) 

    # 2. Subsampling (Crucial for very large datasets)
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
            print(f"Using all {n_obs} cells (n_jobs={n_jobs} for parallel computation)...")
        high_dim_data = adata.obsm[high_dim_key]
        low_dim_data = adata.obsm[low_dim_key]
    
    n_points = high_dim_data.shape[0]

    # 3. Compute Pairwise Distances with parallelism and progress tracking
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

    # 4. Extract upper triangle indices (excluding diagonal)
    if verbose:
        print("Extracting pairwise distances...")
    triu_indices = np.triu_indices_from(high_dim_dists, k=1)
    H = high_dim_dists[triu_indices]
    L = low_dim_dists[triu_indices]

    # 5. Normalize Distances (CRITICAL STEP)
    if verbose:
        print("Normalizing distances...")
    H_norm = H / np.max(H)
    L_norm = L / np.max(L)

    # 6. Calculate Stress (Stress-1 formula)
    if verbose:
        print("Calculating final stress score...")
    numerator = np.sum(np.square(H_norm - L_norm))
    denominator = np.sum(np.square(L_norm))
    
    stress = np.sqrt(numerator / denominator)
    
    return stress


def _compute_distances_with_progress(data, n_jobs=-1, batch_size=2000, desc="Computing distances", verbose=True):
    """
    Compute pairwise distances with batch-wise progress tracking and parallelism.
    
    Parameters:
    -----------
    data : array-like of shape (n_samples, n_features)
        The data matrix.
    n_jobs : int
        Number of parallel jobs (-1 = all cores).
    batch_size : int
        Size of each batch for progress updates.
    desc : str
        Description for the progress bar.
    verbose : bool
        Whether to show progress bar.
        
    Returns:
    --------
    ndarray of shape (n_samples, n_samples)
        Pairwise distance matrix.
    """
    n_samples = data.shape[0]
    
    # For small datasets, just compute directly
    if n_samples <= batch_size or not verbose:
        return pairwise_distances(data, n_jobs=n_jobs)
    
    # For larger datasets, compute in batches with progress bar
    n_batches = (n_samples + batch_size - 1) // batch_size
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    with tqdm(total=n_batches, desc=desc, disable=not verbose) as pbar:
        for i in range(0, n_samples, batch_size):
            end_i = min(i + batch_size, n_samples)
            # Compute distances from batch to all points (parallelized)
            distances[i:end_i, :] = pairwise_distances(
                data[i:end_i], 
                data, 
                n_jobs=n_jobs
            )
            pbar.update(1)
    
    # Ensure diagonal is exactly 0 (floating point fix)
    np.fill_diagonal(distances, 0)
    
    return distances