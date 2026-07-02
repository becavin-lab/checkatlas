import numpy as np
import scanpy as sc
from sklearn.metrics import pairwise_distances
from tqdm import tqdm


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
        batch_size (int): Batch size for row-chunked reads of TriangularMatrix inputs.
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
                f"Using sampling-based stress calculation "
                f"(sampling {sample_pairs:,} of {n_pairs_total:,} pairs)..."
            )
        stress = _sampled_kruskal_stress(
            high_dists, low_dists, sample_pairs, seed, verbose
        )
    else:
        if verbose:
            print(f"Computing full Kruskal stress ({n_pairs_total:,} pairs)...")
        stress = _chunked_kruskal_stress(
            high_dists, low_dists, batch_size, verbose
        )

    return stress


def _upper_tri_row_starts(n: int) -> np.ndarray:
    """Return the start index of each row in a flat row-major upper-triangle buffer.

    ``starts[i]`` is the flat index where ``D[i, i+1]`` begins.  Used to
    convert a flat upper-triangle index into a (row, col) pair.
    """
    starts = np.empty(n + 1, dtype=np.int64)
    starts[0] = 0
    if n > 1:
        counts = np.arange(n - 1, 0, -1, dtype=np.int64)  # length n-1
        starts[1 : 1 + len(counts)] = np.cumsum(counts)
    starts[n] = n * (n - 1) // 2
    return starts


def _flat_to_ij(flat_idx: np.ndarray, n: int) -> tuple:
    """Convert a flat upper-triangle index into a (row, col) pair.

    Inverse of the row-major formula used in :class:`TriangularMatrix`:
    given ``flat_idx = row*n - row*(row+1)/2 + (col - row - 1)``,
    recover ``(row, col)``.
    """
    row_starts = _upper_tri_row_starts(n)
    rows = np.searchsorted(row_starts, flat_idx, side="right") - 1
    rows = np.clip(rows, 0, n - 2)
    cols = flat_idx - row_starts[rows] + rows + 1
    return rows.astype(np.int64), cols.astype(np.int64)


def _gather_pairs(dists, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """Element-wise distance lookup that works for ndarray, memmap, and TriangularMatrix."""
    if hasattr(dists, "_get_elements"):
        return dists._get_elements(rows, cols)
    return np.asarray(dists)[rows, cols]


def _sampled_kruskal_stress(high_dists, low_dists, n_pairs, seed, verbose=True):
    """Compute Kruskal stress by sampling distance pairs (vectorised)."""
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

    H = _gather_pairs(high_dists, rows, cols).astype(np.float32, copy=False)
    L = _gather_pairs(low_dists, rows, cols).astype(np.float32, copy=False)

    if verbose:
        print("Normalizing distances...")

    H_max = float(H.max())
    L_max = float(L.max())
    H_norm = H / H_max if H_max > 0 else H
    L_norm = L / L_max if L_max > 0 else L

    if verbose:
        print("Calculating stress score...")

    diff_sq = (H_norm - L_norm) ** 2
    numerator = float(diff_sq.sum())
    denominator = float((L_norm ** 2).sum())

    if denominator == 0:
        return 0.0
    return float(np.sqrt(numerator / denominator))


def _chunked_kruskal_stress(high_dists, low_dists, batch_size, verbose=True):
    """Compute Kruskal stress on the full upper triangle via row-chunked reads.

    For each row block we issue a single block-read on the distance
    matrices, so the per-block cost is one C-level call regardless of
    the block's row count.  The row-loop that the previous version used
    (one row at a time) is replaced by a single loop over row-blocks.
    """
    n = high_dists.shape[0]

    if verbose:
        print("Finding max distances for normalization (row-chunked)...")

    H_max = 0.0
    L_max = 0.0
    for i0 in tqdm(range(0, n, batch_size), desc="  H/L max", disable=not verbose):
        i1 = min(n, i0 + batch_size)
        H_block = _row_block_upper(high_dists, i0, i1)
        L_block = _row_block_upper(low_dists, i0, i1)
        H_max = max(H_max, float(H_block.max()))
        L_max = max(L_max, float(L_block.max()))

    if H_max <= 0:
        H_max = 1.0
    if L_max <= 0:
        L_max = 1.0

    if verbose:
        print("Computing stress values (row-chunked)...")

    numerator = 0.0
    denominator = 0.0
    for i0 in tqdm(range(0, n, batch_size), desc="  stress", disable=not verbose):
        i1 = min(n, i0 + batch_size)
        H_block = _row_block_upper(high_dists, i0, i1) / H_max
        L_block = _row_block_upper(low_dists, i0, i1) / L_max
        diff = H_block - L_block
        numerator += float((diff * diff).sum())
        denominator += float((L_block * L_block).sum())

    if denominator == 0:
        return 0.0
    return float(np.sqrt(numerator / denominator))


def _row_block_upper(dists, i0: int, i1: int) -> np.ndarray:
    """Return the strict upper-triangle of rows ``[i0, i1)`` as a 1-D array.

    The returned array is the row-major concatenation of
    ``D[i, i+1], D[i, i+2], ..., D[i, n-1]`` for ``i = i0, ..., i1-1``.
    For a :class:`TriangularMatrix` input this maps to a single
    contiguous read of the underlying memmap for those rows.  For a
    dense ndarray the strict upper triangle is sliced in-place.

    The 1-D layout matches :func:`np.triu_indices(n, k=1)`'s row-major
    order so the metric can sum / mean over the result directly.
    """
    n = dists.shape[0]
    if i0 >= n:
        return np.zeros(0, dtype=np.float32)
    if i1 > n:
        i1 = n

    block = _block(dists, i0, i1, i0, n)
    if block.size == 0:
        return np.zeros(0, dtype=np.float32)

    # block shape: (i1-i0, n-i0).  The strict upper triangle for rows
    # [i0, i1) is the elements block[i-i0, j-i0] for j > i.  Since the
    # columns j=i0..n-1 correspond to block columns 0..(n-1-i0), the
    # strict-upper mask is j > i ⇔ (j-i0) > (i-i0).
    n_cols = n - i0
    col_idx = np.arange(n_cols)[None, :]   # (1, n_cols)
    row_off = np.arange(i1 - i0)[:, None]  # (i1-i0, 1)
    mask = col_idx > row_off
    return block[mask].astype(np.float32, copy=False)


def _block(dists, i0: int, i1: int, j0: int, j1: int) -> np.ndarray:
    """Read a (rows, cols) block.  Works for ndarray, memmap, and TriangularMatrix."""
    if hasattr(dists, "_get_block"):
        rows = np.arange(i0, i1, dtype=np.int64)
        cols = np.arange(j0, j1, dtype=np.int64)
        return np.asarray(dists._get_block(rows, cols))
    if hasattr(dists, "_to_dense_f16") and not isinstance(dists, np.ndarray):
        return dists._to_dense_f16()[i0:i1, j0:j1].astype(np.float32, copy=False)
    return np.asarray(dists)[i0:i1, j0:j1]


def _compute_distances_chunked(
    data, n_jobs=-1, batch_size=1000, desc="Computing distances", verbose=True
):
    """Compute pairwise distances in chunks."""
    n_samples = data.shape[0]

    distances = np.zeros((n_samples, n_samples), dtype=np.float32)

    for i in tqdm(range(0, n_samples, batch_size), desc=desc, disable=not verbose):
        end = min(i + batch_size, n_samples)
        distances[i:end, :] = pairwise_distances(data[i:end], data, n_jobs=n_jobs)

    np.fill_diagonal(distances, 0)
    return distances
