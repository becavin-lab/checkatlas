import gc

import numpy as np
import scanpy as sc
from joblib import Parallel, delayed
from sklearn.metrics import pairwise_distances
from tqdm import tqdm, trange


def run(
    adata,
    low_dim_key="X_umap",
    high_dim_key="X",
    k_neighbors=30,
    n_samples=None,
    seed=42,
    verbose=True,
    n_jobs=-1,
    precomputed_high_dists=None,
    precomputed_low_dists=None,
    precomputed_high_knn=None,
    precomputed_low_knn=None,
    precomputed_high_knn_dists=None,
    precomputed_low_knn_dists=None,
    use_memmap=True,
):
    """
    Continuity (Memory-Optimized)
    Measures the preservation of original neighbors from the high-dimensional space in the low-dimensional embedding.
    It penalizes "broken trajectories" or missing neighbors (points that are neighbors in high-dim but not in low-dim).

    Optimized for memory usage by processing rows in parallel without storing full N*N rank matrices.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation. None = all.
        seed (int): Random seed for reproducibility.
        verbose (bool): Whether to print progress.
        n_jobs (int): Number of parallel jobs for computation.
        precomputed_high_dists (np.ndarray or memmap): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray or memmap): Precomputed low-dim distance matrix.
        use_memmap (bool): Hint that we are using memory mapped files.

    Returns:
        float: The Continuity score.

    Interpretation:
        Range 0 to 1.
        Higher is better (1 means perfect continuity).
    """

    # Determine number of workers
    if n_jobs == -1:
        import os

        n_workers = os.cpu_count() or 4
    else:
        n_workers = max(1, n_jobs)

    # 1. Use Precomputed Distances if available
    if precomputed_high_dists is not None and precomputed_low_dists is not None:
        if verbose:
            print("Using precomputed distance matrices...")
        high_dists = precomputed_high_dists
        low_dists = precomputed_low_dists
        n_cells = high_dists.shape[0]

    else:
        # Fallback to local computation
        if verbose:
            print("Precomputed distances not provided. Calculating locally...")

        # Check keys
        if low_dim_key not in adata.obsm.keys():
            if verbose:
                print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        # Prepare Data
        n_obs = adata.n_obs
        if n_samples is not None and n_samples < n_obs:
            if verbose:
                print(f"Subsampling {n_samples} cells...")
            np.random.seed(seed)
            indices = np.random.choice(n_obs, n_samples, replace=False)
        else:
            indices = np.arange(n_obs)

        # High Dim Data
        if high_dim_key == "X":
            high_dim_data = adata.X[indices]
            if hasattr(high_dim_data, "toarray"):
                high_dim_data = high_dim_data.toarray()
        else:
            if high_dim_key not in adata.obsm.keys():
                if verbose:
                    print(f"Calculating {high_dim_key}...")
                sc.tl.pca(adata, random_state=seed)
            high_dim_data = adata.obsm[high_dim_key][indices]

        low_dim_data = adata.obsm[low_dim_key][indices]
        n_cells = high_dim_data.shape[0]

        # Compute Distances
        if verbose and n_cells > 10000:
            print(
                "Warning: Running on large dataset without precomputed memmap distances."
            )

        high_dists = pairwise_distances(high_dim_data, n_jobs=n_workers)
        low_dists = pairwise_distances(low_dim_data, n_jobs=n_workers)

    # 2. Fast path: precomputed kNN available
    # Use precomputed high-dim neighbours + direct row access (no to_dense, no joblib)
    if (
        precomputed_high_dists is not None
        and precomputed_low_dists is not None
        and precomputed_high_knn is not None
    ):
        if verbose:
            print(f"Using precomputed kNN fast path (k={k_neighbors})...")

        # Continuity: high-dim neighbours, check ranks in low-dim
        high_neighbors = precomputed_high_knn[:, 1 : k_neighbors + 1]  # skip self

        _total_penalty = 0.0
        for i in trange(n_cells, desc="Continuity", disable=not verbose):
            l_row = low_dists[i]
            nbrs = high_neighbors[i].astype(np.int64, copy=False)
            nbr_dists = l_row[nbrs]
            for dist_j in nbr_dists:
                r = np.count_nonzero(l_row < dist_j)
                if r > k_neighbors:
                    _total_penalty += r - k_neighbors

        max_penalty = n_cells * k_neighbors * (2 * n_cells - 3 * k_neighbors - 1) / 2.0
        if max_penalty == 0:
            return 1.0
        score = 1.0 - (2.0 / max_penalty) * _total_penalty
        return max(0.0, min(1.0, score))

    # 3. Calculate Continuity Row-by-Row (original path)
    if verbose:
        print(f"Computing Continuity (k={k_neighbors}) row-by-row...")

    batch_size = 100
    batches = [
        range(i, min(i + batch_size, n_cells)) for i in range(0, n_cells, batch_size)
    ]

    def _process_batch(row_indices):
        batch_penalty = 0.0
        for i in row_indices:
            # 1. Get neighbors in HIGH dimension
            h_row = high_dists[i]
            # indices of k+1 smallest distances
            knn_indices = np.argpartition(h_row, k_neighbors + 1)[: k_neighbors + 1]

            # Sort to find self
            candidate_dists = h_row[knn_indices]
            sorted_args = np.argsort(candidate_dists)
            knn_indices_sorted = knn_indices[sorted_args]

            # Exclude self
            if knn_indices_sorted[0] == i:
                neighbors = knn_indices_sorted[1 : k_neighbors + 1]
            else:
                neighbors = knn_indices[knn_indices != i][:k_neighbors]

            # 2. Calculate rank of these neighbors in LOW dimension
            l_row = low_dists[i]

            neighbor_dists_low = l_row[neighbors]

            for dist_j in neighbor_dists_low:
                # Count points closer than j in low dim
                r = np.count_nonzero(l_row < dist_j)

                # If rank > k, penalty
                if r > k_neighbors:
                    batch_penalty += r - k_neighbors

        return batch_penalty

    # Execute parallel
    penalties = Parallel(n_jobs=n_workers)(
        delayed(_process_batch)(batch)
        for batch in tqdm(batches, desc="Computing penalties", disable=not verbose)
    )

    total_penalty = sum(penalties)

    # 3. Normalize Score
    if n_cells <= k_neighbors + 1:
        return 1.0

    max_penalty = n_cells * k_neighbors * (2 * n_cells - 3 * k_neighbors - 1) / 2.0

    if max_penalty == 0:
        return 1.0

    score = 1.0 - (2.0 / max_penalty) * total_penalty
    return max(0.0, min(1.0, score))
