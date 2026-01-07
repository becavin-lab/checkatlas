from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.metrics import pairwise_distances
import numpy as np
from tqdm import tqdm


def run(count_repr, annotations, sample_size=None, n_jobs=None, random_state=42, 
        varbose=True, batch_size=1000):
    """
    Calculate the Average Silhouette Width (Silhouette Score).
    
    :param count_repr: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points
    :param annotations: array-like of shape (n_samples,)
        Cluster labels for each sample
    :param sample_size: int, optional (default=None)
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data. If None, no sampling is used.
    :param n_jobs: int, optional (default=None)
        Number of parallel jobs for distance computation.
        -1 means using all processors. Only used when sample_size is None.
    :param random_state: int, optional (default=42)
        Random state for reproducibility when sampling.
    :param varbose: bool, optional (default=True)
        Whether to show tqdm progress bar during computation.
    :param batch_size: int, optional (default=1000)
        Batch size for progress tracking during distance computation.
    :return: float
        The mean Silhouette Coefficient over all samples.
        Range: [-1, 1], where higher values indicate better clustering.
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)
    n_samples = count_repr.shape[0]
    
    if sample_size is not None:
        # Use sklearn's built-in sampling (faster, approximate)
        return silhouette_score(
            count_repr, 
            annotations, 
            sample_size=sample_size,
            random_state=random_state
        )
        
    elif n_jobs is not None and varbose:
        # Batch-wise distance computation with progress tracking
        n_batches = (n_samples + batch_size - 1) // batch_size
        
        # Compute full distance matrix in batches with progress bar
        distances = np.zeros((n_samples, n_samples), dtype=np.float32)
        
        with tqdm(total=n_batches, desc="Computing pairwise distances", 
                  disable=not varbose) as pbar:
            for i in range(0, n_samples, batch_size):
                end_i = min(i + batch_size, n_samples)
                # Compute distances from batch to all points
                distances[i:end_i, :] = pairwise_distances(
                    count_repr[i:end_i], 
                    count_repr, 
                    n_jobs=n_jobs
                )
                pbar.update(1)
        
        # Fix floating-point precision issues on diagonal
        np.fill_diagonal(distances, 0)
        
        # Now compute silhouette with precomputed distances
        if varbose:
            print("Computing silhouette scores...")
        return silhouette_score(distances, annotations, metric="precomputed")
    
    elif n_jobs is not None:
        # Precompute distances with parallel processing (no progress bar)
        distances = pairwise_distances(count_repr, n_jobs=n_jobs)
        np.fill_diagonal(distances, 0)  # Fix floating-point precision issues
        return silhouette_score(distances, annotations, metric="precomputed")
    
    else:
        # Default: no parallelization, no sampling
        return silhouette_score(count_repr, annotations)