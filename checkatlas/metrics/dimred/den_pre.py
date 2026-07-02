import numpy as np
from scipy.stats import pearsonr


def run(
    adata,
    low_dim_key="X_umap",
    high_dim_key="X",
    k_neighbors=30,
    n_samples=5000,
    seed=42,
    n_jobs=-1,
    verbose=True,
    log_transform=False,
    precomputed_high_knn_dists=None,
    precomputed_low_knn_dists=None,
):
    """
    Density Preservation
    Measures the preservation of local density by calculating the Pearson correlation between the radii of k-nearest neighbors in high-dimensional and low-dimensional spaces.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs. (Unused: computation is vectorised.)
        verbose (bool): Whether to print progress.
        log_transform (bool): Whether to log-transform the radii before correlation.
        precomputed_high_knn_dists (np.ndarray): Precomputed k-NN distances for high-dim data,
            shape (N, k_neighbors + 1) including the self-distance.
        precomputed_low_knn_dists (np.ndarray): Precomputed k-NN distances for low-dim data,
            shape (N, k_neighbors + 1) including the self-distance.

    Returns:
        float: The Density Preservation score (Pearson Correlation).

    Interpretation:
        Range -1 to 1.
        Higher is better (1 means perfect preservation of relative density).
    """

    if precomputed_high_knn_dists is None or precomputed_low_knn_dists is None:
        raise ValueError(
            "den_pre.run requires precomputed_high_knn_dists and "
            "precomputed_low_knn_dists. Run via cal_dimred or pass them "
            "explicitly. The metric does not recompute neighbours locally "
            "because the radii are a function of the (N, k) kNN distance "
            "matrix which cal_dimred already materialises."
        )

    if verbose:
        print("Using precomputed k-NN distances for density preservation...")

    dists_high = np.asarray(precomputed_high_knn_dists)
    dists_low = np.asarray(precomputed_low_knn_dists)

    if dists_high.shape != dists_low.shape:
        raise ValueError(
            f"kNN dists shape mismatch: high {dists_high.shape} vs low {dists_low.shape}"
        )
    if dists_high.ndim != 2:
        raise ValueError(f"expected 2-D kNN dists; got ndim={dists_high.ndim}")
    if dists_high.shape[1] < k_neighbors:
        raise ValueError(
            f"kNN dists have only {dists_high.shape[1]} columns; "
            f"need at least k_neighbors={k_neighbors}"
        )

    radii_high = dists_high[:, k_neighbors - 1].astype(np.float64, copy=False)
    radii_low = dists_low[:, k_neighbors - 1].astype(np.float64, copy=False)

    if log_transform:
        radii_high = np.log1p(radii_high)
        radii_low = np.log1p(radii_low)

    if verbose:
        print("Calculating correlation of radii...")

    if radii_high.std() == 0 or radii_low.std() == 0:
        return 0.0
    correlation, _ = pearsonr(radii_high, radii_low)
    return float(correlation)
