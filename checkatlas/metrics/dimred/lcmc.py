import numpy as np
import scanpy as sc

from ._overlap import _overlap_count


def run(
    adata,
    low_dim_key="X_umap",
    high_dim_key="X",
    k_neighbors=30,
    n_samples=5000,
    seed=42,
    n_jobs=-1,
    verbose=True,
    precomputed_high_knn=None,
    precomputed_low_knn=None,
):
    """
    LCMC (Local Continuity Meta-Criterion)
    Measures the overlap of k-nearest neighbors between high-dimensional and low-dimensional spaces, adjusted for chance.
    It quantifies how well the local neighborhood is preserved relative to random chance.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs. (Unused: computation is vectorised.)
        verbose (bool): Whether to print progress.
        precomputed_high_knn (np.ndarray): Precomputed k-NN indices for high-dim data.
        precomputed_low_knn (np.ndarray): Precomputed k-NN indices for low-dim data.

    Returns:
        float: The LCMC score.

    Interpretation:
        Range 0 to 1.
        Higher is better (1 means perfect preservation of neighbors).
    """

    if precomputed_high_knn is not None and precomputed_low_knn is not None:
        if verbose:
            print("Using precomputed k-NN graphs...")
        indices_high = np.asarray(precomputed_high_knn)[:, 1 : k_neighbors + 1]
        indices_low = np.asarray(precomputed_low_knn)[:, 1 : k_neighbors + 1]
    else:
        if verbose:
            print("Precomputed k-NN not provided. Calculating locally...")

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

        from sklearn.neighbors import NearestNeighbors

        query_k = k_neighbors + 1

        if verbose:
            print(f"Finding neighbors (k={k_neighbors}) in High-Dim...")
        nbrs_high = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(
            high_dim_data
        )
        _, indices_high = nbrs_high.kneighbors(high_dim_data)
        indices_high = indices_high[:, 1 : k_neighbors + 1]

        low_dim_data = adata.obsm[low_dim_key][sample_idx]

        if verbose:
            print(f"Finding neighbors (k={k_neighbors}) in Low-Dim...")
        nbrs_low = NearestNeighbors(n_neighbors=query_k, n_jobs=n_jobs).fit(
            low_dim_data
        )
        _, indices_low = nbrs_low.kneighbors(low_dim_data)
        indices_low = indices_low[:, 1 : k_neighbors + 1]

    if verbose:
        print("Calculating LCMC (vectorised)...")

    n_cells = indices_high.shape[0]
    intersection = _overlap_count(indices_low, indices_high).astype(np.float64)
    mean_overlap = (intersection.sum() / n_cells) / k_neighbors

    baseline = k_neighbors / (n_cells - 1) if n_cells > 1 else 0.0
    denom = 1.0 - baseline
    if denom == 0:
        return 0.0
    return float((mean_overlap - baseline) / denom)
