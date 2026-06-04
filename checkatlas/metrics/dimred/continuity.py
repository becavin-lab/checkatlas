import numpy as np
import scanpy as sc

from ._rank_penalty import rank_penalty, self_excluded_knn


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
    rank_backend="auto",
):
    """
    Continuity (GPU-accelerated O(N^2 log N) inner loop).

    Symmetric twin of :func:`trustworthiness.run`.  Penalises "broken
    trajectories" — points that are neighbours in the original
    high-dim space but whose rank in the low-dim embedding exceeds
    ``k_neighbors``.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbours to consider.
        n_samples (int): Number of samples to subsample for calculation.
            ``None`` = all cells.
        seed (int): Random seed for reproducibility.
        verbose (bool): Whether to print progress.
        n_jobs (int): Used only for the standalone fallback path.
        precomputed_high_dists, precomputed_low_dists:
            Precomputed distance matrices.  Materialised once to a
            dense float32 buffer (cached on ``TriangularMatrix``).
        precomputed_high_knn, precomputed_high_knn_dists,
        precomputed_low_knn, precomputed_low_knn_dists:
            Precomputed k-NN index / distance arrays.  For continuity
            the only one that is read is ``precomputed_high_knn``;
            the others are accepted for API symmetry with
            :func:`trustworthiness.run`.
        use_memmap (bool): Hint, currently informational.
        rank_backend (str): ``"auto"``, ``"jax_single_shot"``,
            ``"jax_chunked"``, or ``"cpu"``.

    Returns:
        float: The Continuity score in [0, 1].
        Higher is better (1 means perfect continuity).
    """
    # 1. Resolve distance matrices
    if precomputed_high_dists is not None and precomputed_low_dists is not None:
        if verbose:
            print("Using precomputed distance matrices...")
        high_dists = precomputed_high_dists
        low_dists = precomputed_low_dists
        n_cells = high_dists.shape[0]
    else:
        # ── Standalone fallback: compute distances from the AnnData ──
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
            indices = np.random.choice(n_obs, n_samples, replace=False)
        else:
            indices = np.arange(n_obs)

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

        try:
            import psutil

            avail = psutil.virtual_memory().available
        except ImportError:
            avail = 16 * 1024**3

        needed = 2 * n_cells * n_cells * 4
        if needed > 0.8 * avail and verbose:
            print(
                f"Warning: standalone path will allocate ≈{needed / 1024**3:.1f} GB "
                f"of distance matrices (system has {avail / 1024**3:.1f} GB free). "
                "Consider running via cal_dimred to use the GPU + memmap path."
            )

        from sklearn.metrics import pairwise_distances

        high_dists = pairwise_distances(high_dim_data, n_jobs=n_jobs)
        low_dists = pairwise_distances(low_dim_data, n_jobs=n_jobs)

    # 2. Resolve HIGH-dim k-NN indices (precomputed or computed)
    if precomputed_high_knn is not None:
        high_neighbors = np.asarray(precomputed_high_knn)[:, 1 : k_neighbors + 1]
    else:
        if verbose:
            print("Computing high-dim k-NN...")
        high_neighbors = self_excluded_knn(np.asarray(high_dists), k=k_neighbors)

    # 3. Single O(N^2 log N) call.  Continuity is the *mirror* of
    #    trustworthiness: it uses the HIGH-dim neighbours and ranks
    #    them in the LOW-dim distance matrix.
    if verbose:
        print(f"Computing Continuity (k={k_neighbors}) via rank_penalty...")

    total_penalty = rank_penalty(
        low_dists,
        high_neighbors,
        k=k_neighbors,
        use_jax=rank_backend,
    )

    # 4. Normalise to [0, 1]
    if n_cells <= k_neighbors + 1:
        return 1.0
    max_penalty = n_cells * k_neighbors * (2 * n_cells - 3 * k_neighbors - 1) / 2.0
    if max_penalty == 0:
        return 1.0
    score = 1.0 - (2.0 / max_penalty) * total_penalty
    return max(0.0, min(1.0, score))
