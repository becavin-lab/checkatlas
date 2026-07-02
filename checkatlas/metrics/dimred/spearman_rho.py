import numpy as np
import scanpy as sc
import scipy.stats
from sklearn.metrics import pairwise_distances
from tqdm import tqdm

from .kruskal_stress import (
    _chunked_kruskal_stress,
    _compute_distances_chunked,
    _flat_to_ij,
    _gather_pairs,
    _row_block_upper,
)


def run(
    adata,
    low_dim_key="X_umap",
    high_dim_key="X",
    n_samples=None,
    seed=42,
    n_jobs=-1,
    verbose=True,
    batch_size=1000,
    sample_pairs=1000000,
    precomputed_high_dists=None,
    precomputed_low_dists=None,
):
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
        batch_size (int): Batch size for row-chunked reads of TriangularMatrix inputs.
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

    if precomputed_high_dists is not None and precomputed_low_dists is not None:
        if verbose:
            print("Using precomputed distance matrices...")
        high_dists = precomputed_high_dists
        low_dists = precomputed_low_dists
    else:
        if verbose:
            print("Precomputed distances not provided. Calculating locally...")

        if low_dim_key not in adata.obsm.keys():
            if verbose:
                print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        n_obs = adata.n_obs
        if n_samples is not None and n_samples < n_obs:
            if verbose:
                print(f"Subsampling {n_samples} cells...")
            np.random.seed(seed)
            sample_idx = np.random.choice(n_obs, n_samples, replace=False)
        else:
            sample_idx = np.arange(n_obs)

        if high_dim_key == "X":
            high_dim_data = adata.X[sample_idx]
            if hasattr(high_dim_data, "toarray"):
                high_dim_data = high_dim_data.toarray()
        else:
            if high_dim_key not in adata.obsm.keys():
                if verbose:
                    print(f"Calculating {high_dim_key}...")
                sc.tl.pca(adata, random_state=seed)
            high_dim_data = adata.obsm[high_dim_key][sample_idx]

        low_dim_data = adata.obsm[low_dim_key][sample_idx]

        high_dists = _compute_distances_chunked(
            high_dim_data, n_jobs, batch_size, "High-Dim Dists", verbose
        )
        low_dists = _compute_distances_chunked(
            low_dim_data, n_jobs, batch_size, "Low-Dim Dists", verbose
        )

    n = high_dists.shape[0]
    n_pairs_total = n * (n - 1) // 2

    if sample_pairs is not None and sample_pairs < n_pairs_total:
        if verbose:
            print(
                f"Using sampling-based Spearman "
                f"(sampling {sample_pairs:,} of {n_pairs_total:,} pairs)..."
            )
        rho = _sampled_spearman(high_dists, low_dists, sample_pairs, seed, verbose)
    else:
        if verbose:
            print(f"Computing full Spearman correlation ({n_pairs_total:,} pairs)...")
        rho = _chunked_spearman(high_dists, low_dists, batch_size, verbose)

    return rho


def _sampled_spearman(high_dists, low_dists, n_pairs, seed, verbose=True):
    """Compute Spearman correlation by sampling distance pairs (vectorised)."""
    n = high_dists.shape[0]
    n_pairs_total = n * (n - 1) // 2
    n_pairs = min(n_pairs, n_pairs_total)

    if verbose:
        print("Sampling random pair indices...")

    rng = np.random.default_rng(seed)
    flat_idx = rng.choice(n_pairs_total, size=n_pairs, replace=False)
    flat_idx.sort()
    rows, cols = _flat_to_ij(flat_idx, n)

    if verbose:
        print("Extracting sampled distance pairs...")

    H = _gather_pairs(high_dists, rows, cols).astype(np.float64, copy=False)
    L = _gather_pairs(low_dists, rows, cols).astype(np.float64, copy=False)

    if verbose:
        print("Computing Spearman correlation...")

    if H.std() == 0 or L.std() == 0:
        return 0.0
    rho = scipy.stats.spearmanr(H, L).statistic
    if rho is None or np.isnan(rho):
        return 0.0
    return float(rho)


def _chunked_spearman(high_dists, low_dists, batch_size, verbose=True):
    """Compute Spearman correlation on the full upper triangle via row-chunked reads.

    Mirrors :func:`_chunked_kruskal_stress` but keeps the per-block
    distance vectors in a list and concatenates at the end.  Uses the
    same row-chunked read path so the previous per-row overhead is
    avoided.
    """
    n = high_dists.shape[0]

    H_list = []
    L_list = []
    for i0 in tqdm(range(0, n, batch_size), desc="Extracting distances", disable=not verbose):
        i1 = min(n, i0 + batch_size)
        H_list.append(_row_block_upper(high_dists, i0, i1))
        L_list.append(_row_block_upper(low_dists, i0, i1))

    H_vec = np.concatenate(H_list).astype(np.float64, copy=False)
    L_vec = np.concatenate(L_list).astype(np.float64, copy=False)

    if verbose:
        print("Computing Spearman correlation...")

    if H_vec.std() == 0 or L_vec.std() == 0:
        return 0.0
    rho = scipy.stats.spearmanr(H_vec, L_vec).statistic
    if rho is None or np.isnan(rho):
        return 0.0
    return float(rho)
