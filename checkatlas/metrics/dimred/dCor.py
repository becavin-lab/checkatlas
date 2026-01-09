import numpy as np
import scanpy as sc
from sklearn.metrics import pairwise_distances
from tqdm import tqdm
import logging

# # Configure logger
# logger = logging.getLogger("checkatlas")
# logger.setLevel(logging.DEBUG)

def run(adata, low_dim_key='X_umap', high_dim_key='X_pca', n_samples=5000, seed=42, 
        n_jobs=-1, verbose=True, batch_size=2000):
    """
    Computes Distance Correlation (dCor) to assess dimensionality reduction.
    
    dCor ranges from 0 to 1. 
    It measures the dependence between the pairwise distances in High-Dim 
    and the pairwise distances in Low-Dim. 
    Unlike Pearson, dCor=0 implies true independence.

    Parameters:
    -----------
    adata : AnnData
        The annotated data matrix.
    low_dim_key : str
        Key in adata.obsm (default: 'X_umap').
    high_dim_key : str
        Key in adata.obsm (default: 'X_pca').
    n_samples : int or None
        Number of cells to subsample. Computing full distance matrices 
        and performing double-centering is O(N^2). 
        5000 is recommended.
    seed : int
        Random seed.
    n_jobs : int
        Parallel jobs for distance calc.
    verbose : bool
        Show progress.

    Returns:
    --------
    float
        Distance Correlation (0 to 1).
        Higher is better (stronger dependence between original and embedded geometry).
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
            print(f"Using all {n_obs} cells (Warning: High memory usage)...")
        high_dim_data = adata.obsm[high_dim_key]
        low_dim_data = adata.obsm[low_dim_key]

    # 3. Compute Pairwise Distance Matrices
    # We need the full matrices for double centering
    high_dim_dists = _compute_distances_with_progress(
        high_dim_data, n_jobs=n_jobs, batch_size=batch_size, 
        desc="Computing High-Dim Distances", verbose=verbose
    )
    
    low_dim_dists = _compute_distances_with_progress(
        low_dim_data, n_jobs=n_jobs, batch_size=batch_size, 
        desc="Computing Low-Dim Distances", verbose=verbose
    )

    # 4. Compute Distance Correlation
    if verbose:
        print("Calculating Distance Correlation (Double Centering)...")

    dcor_score = _distance_correlation(high_dim_dists, low_dim_dists)

    # if verbose:
    #     print(f"Distance Correlation: {dcor_score:.4f}")
    #     logger.debug(f"Distance Correlation : {dcor_score}")

    return dcor_score


def _distance_correlation(A, B):
    """
    Computes the Distance Correlation between two distance matrices A and B.
    A and B must be square matrices of the same size.
    """
    # 1. Double Centering
    # A_centered = A - row_mean - col_mean + grand_mean
    A_centered = _double_center(A)
    B_centered = _double_center(B)

    # 2. Compute Distance Covariance (dCov)
    # The mean of the product of the centered matrices
    dCov_sq = np.mean(A_centered * B_centered)

    # 3. Compute Distance Variances (dVar)
    dVar_A_sq = np.mean(A_centered * A_centered)
    dVar_B_sq = np.mean(B_centered * B_centered)

    # 4. Compute Distance Correlation
    # If variance is 0, correlation is 0
    if dVar_A_sq > 0 and dVar_B_sq > 0:
        dCor = np.sqrt(dCov_sq / (np.sqrt(dVar_A_sq) * np.sqrt(dVar_B_sq)))
    else:
        dCor = 0.0

    return dCor


def _double_center(D):
    """
    Performs double centering on a distance matrix D.
    Formula: A_ij = d_ij - mean(row_i) - mean(col_j) + mean(grand)
    """
    n = D.shape[0]
    
    # Mean of each row
    row_means = np.mean(D, axis=1, keepdims=True)
    # Mean of each column
    col_means = np.mean(D, axis=0, keepdims=True)
    # Grand mean of the whole matrix
    grand_mean = np.mean(D)
    
    # Broadcasting makes this efficient
    centered = D - row_means - col_means + grand_mean
    
    # There is a slight adjustment for unbiased estimators in literature,
    # but for large N (like scRNAseq), the standard centering is sufficient.
    return centered


def _compute_distances_with_progress(data, n_jobs=-1, batch_size=2000, desc="Computing distances", verbose=True):
    """Helper to compute distances with progress bar"""
    n_samples = data.shape[0]
    if n_samples <= batch_size or not verbose:
        return pairwise_distances(data, n_jobs=n_jobs)
    
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    n_batches = int(np.ceil(n_samples / batch_size))
    
    with tqdm(total=n_batches, desc=desc, disable=not verbose) as pbar:
        for i in range(0, n_samples, batch_size):
            end_i = min(i + batch_size, n_samples)
            distances[i:end_i, :] = pairwise_distances(data[i:end_i], data, n_jobs=n_jobs)
            pbar.update(1)
            
    np.fill_diagonal(distances, 0)
    return distances