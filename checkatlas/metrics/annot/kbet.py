import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.stats import chi2
from scipy.sparse import issparse


def run(adata, batch_label="batch", k=25, alpha=0.05, n_jobs=-1, verbose=True):
    """
    Calculate kBET rejection rate.
    
    kBET is a statistical test that determines whether the batch composition
    in local neighborhoods matches the global batch composition. Lower rejection
    rates indicate better batch mixing.
    
    `kBET readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/kbet/>`__

    :param adata: AnnData object
        Annotated data matrix
    :param batch_label: str, default='batch'
        Key in adata.obs containing batch labels
    :param k: int, default=25
        Number of nearest neighbors to consider
    :param alpha: float, default=0.05
        Significance level for the chi-squared test
    :param n_jobs: int, default=-1
        Number of parallel jobs for kNN computation. -1 uses all cores.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return: float
        kBET rejection rate. Range: [0, 1], where lower values indicate better mixing.
    """
    # Handle sparse matrices
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    
    # 'batch' col check in adata.obs
    if batch_label not in adata.obs:
        print(f"Batch label '{batch_label}' not found in adata.obs")
        return None
    
    batch_labels = np.asarray(adata.obs[batch_label])
    n_samples = X.shape[0]
    
    if verbose:
        print(f"Computing kBET ({n_samples:,} samples, k={k})...")
    
    # Get unique batches and their frequencies
    unique_batches, batch_counts = np.unique(batch_labels, return_counts=True)
    n_batches = len(unique_batches)
    
    if n_batches == 1:
        return 0.0
    
    # Expected frequencies (global batch distribution)
    expected_freqs = batch_counts / n_samples
    
    # Find k-nearest neighbors with parallelism
    k = min(k, n_samples - 1)
    if verbose:
        print(f"  Computing kNN (k={k}, n_jobs={n_jobs})...")
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm='auto', n_jobs=n_jobs)
    nbrs.fit(X)
    _, indices = nbrs.kneighbors(X)
    
    # Remove self-neighbor (first column)
    neighbor_indices = indices[:, 1:]  # shape: (n_samples, k)
    
    if verbose:
        print(f"  Computing chi-squared tests (vectorized)...")
    
    # ── Vectorized chi² computation ──────────────────────────────────
    # Map batch labels to integer codes for fast bincount
    batch_to_code = {b: i for i, b in enumerate(unique_batches)}
    batch_codes = np.array([batch_to_code[b] for b in batch_labels], dtype=np.intp)
    
    # Get batch codes for all neighbors: shape (n_samples, k)
    neighbor_batch_codes = batch_codes[neighbor_indices]
    
    # Count observed batch frequencies per cell using vectorized bincount
    # We need per-row bincount — use a flat approach with offsets
    offsets = np.arange(n_samples) * n_batches  # offset per cell
    flat_codes = (neighbor_batch_codes + offsets[:, np.newaxis]).ravel()
    observed_flat = np.bincount(flat_codes, minlength=n_samples * n_batches)
    observed_counts = observed_flat.reshape(n_samples, n_batches).astype(np.float64)
    
    # Expected counts per cell
    expected_counts = expected_freqs * k  # shape: (n_batches,)
    expected_counts = np.maximum(expected_counts, 1e-10)  # avoid division by zero
    
    # Vectorized chi-squared statistic for all cells at once
    chi2_stats = np.sum((observed_counts - expected_counts) ** 2 / expected_counts, axis=1)
    
    # Vectorized p-values
    df = n_batches - 1
    p_values = 1 - chi2.cdf(chi2_stats, df)
    
    # Count rejections
    rejections = np.sum(p_values < alpha)
    rejection_rate = rejections / n_samples
    
    if verbose:
        print(f"  kBET rejection rate = {rejection_rate:.4f} "
              f"({rejections:,}/{n_samples:,} cells rejected)")
    
    return rejection_rate
