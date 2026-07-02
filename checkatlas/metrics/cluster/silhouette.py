import numpy as np
from sklearn.metrics import pairwise_distances, silhouette_score
from tqdm import tqdm


_TRIANGULAR_DENSE_THRESHOLD = 30_000


def run(
    count_repr,
    annotations,
    sample_size=None,
    n_jobs=-1,
    random_state=42,
    verbose=True,
    batch_size=1000,
    precomputed_dists=None,
):
    """
    Calculate the Average Silhouette Width (Silhouette Score).

    Optimized with parallel distance computation and optional precomputed distances.

    :param count_repr: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points
    :param annotations: array-like of shape (n_samples,)
        Cluster labels for each sample
    :param sample_size: int, optional (default=None)
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data. If None, no sampling is used.
    :param n_jobs: int, optional (default=-1)
        Number of parallel jobs for distance computation.
        -1 means using all processors.
    :param random_state: int, optional (default=42)
        Random state for reproducibility when sampling.
    :param verbose: bool, optional (default=True)
        Whether to show progress during computation.
    :param batch_size: int, optional (default=1000)
        Batch size for chunked distance computation.
    :param precomputed_dists: np.ndarray or TriangularMatrix, optional
        Precomputed pairwise distance matrix. If provided, skips distance computation.
    :return: float
        The mean Silhouette Coefficient over all samples.
        Range: [-1, 1], where higher values indicate better clustering.
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)
    n_samples = count_repr.shape[0]

    # Path 1: Use precomputed distances (fastest, no recomputation)
    if precomputed_dists is not None:
        if hasattr(precomputed_dists, "to_dense"):
            if precomputed_dists.n <= _TRIANGULAR_DENSE_THRESHOLD:
                if verbose:
                    print("TriangularMatrix — densifying for Silhouette (sklearn)...")
                dense = precomputed_dists.to_dense()
                return float(silhouette_score(dense, annotations, metric="precomputed"))
            else:
                if verbose:
                    print("TriangularMatrix — computing Silhouette from precomputed distances...")
                from ..annot.average_silhouette_width import _silhouette_from_triangular

                return _silhouette_from_triangular(precomputed_dists, annotations, verbose=verbose)
        if verbose:
            print("Using precomputed distances for Silhouette...")
        return silhouette_score(precomputed_dists, annotations, metric="precomputed")

    # Path 2: Sampling (fast approximate)
    if sample_size is not None:
        return silhouette_score(
            count_repr,
            annotations,
            sample_size=sample_size,
            random_state=random_state,
        )

    # Path 3: Full computation with parallel batched distances + progress
    if verbose:
        n_batches = (n_samples + batch_size - 1) // batch_size

        # Compute full distance matrix in batches with progress bar
        distances = np.zeros((n_samples, n_samples), dtype=np.float32)

        with tqdm(
            total=n_batches,
            desc="Computing pairwise distances",
            disable=not verbose,
        ) as pbar:
            for i in range(0, n_samples, batch_size):
                end_i = min(i + batch_size, n_samples)
                distances[i:end_i, :] = pairwise_distances(
                    count_repr[i:end_i], count_repr, n_jobs=n_jobs
                )
                pbar.update(1)

        np.fill_diagonal(distances, 0)
        if verbose:
            print("Computing silhouette scores...")
        return silhouette_score(distances, annotations, metric="precomputed")

    else:
        # No progress bar, but still parallel
        distances = pairwise_distances(count_repr, n_jobs=n_jobs)
        np.fill_diagonal(distances, 0)
        return silhouette_score(distances, annotations, metric="precomputed")
