import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.stats import chi2


def run(X, batch_labels, k=25, alpha=0.05):
    """
    Calculate k-nearest neighbor Batch Effect Test (kBET) rejection rate.
    
    kBET is a statistical test that determines whether the batch composition
    in local neighborhoods matches the global batch composition. Lower rejection
    rates indicate better batch mixing.
    
    `kBET readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/kbet/>`__

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix (e.g., PCA embedding, integrated space)
    :param batch_labels: array-like of shape (n_samples,)
        Batch labels for each sample
    :param k: int, default=25
        Number of nearest neighbors to consider
    :param alpha: float, default=0.05
        Significance level for the chi-squared test
    :return: float
        kBET rejection rate. Range: [0, 1], where lower values indicate better mixing.
        
    """
    X = np.asarray(X)
    batch_labels = np.asarray(batch_labels)
    
    n_samples = X.shape[0]
    
    # Get unique batches and their frequencies
    unique_batches, batch_counts = np.unique(batch_labels, return_counts=True)
    n_batches = len(unique_batches)
    
    if n_batches == 1:
        # Only one batch, no batch effect
        return 0.0
    
    # Expected frequencies (global batch distribution)
    expected_freqs = batch_counts / n_samples
    
    # Find k-nearest neighbors
    k = min(k, n_samples - 1)
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm='auto')  # +1 to exclude self
    nbrs.fit(X)
    distances, indices = nbrs.kneighbors(X)
    
    # Perform chi-squared test for each cell
    rejections = 0
    
    for i in range(n_samples):
        # Get neighbors (excluding self)
        neighbor_indices = indices[i, 1:]  # Skip first neighbor (itself)
        neighbor_batches = batch_labels[neighbor_indices]
        
        # Observed frequencies in neighborhood
        observed_counts = np.zeros(n_batches)
        for j, batch in enumerate(unique_batches):
            observed_counts[j] = np.sum(neighbor_batches == batch)
        
        # Expected frequencies for k neighbors
        expected_counts = expected_freqs * k
        
        # Avoid division by zero
        expected_counts = np.maximum(expected_counts, 1e-10)
        
        # Chi-squared statistic
        chi2_stat = np.sum((observed_counts - expected_counts) ** 2 / expected_counts)
        
        # Degrees of freedom
        df = n_batches - 1
        
        # P-value
        p_value = 1 - chi2.cdf(chi2_stat, df)
        
        # Reject if p-value < alpha
        if p_value < alpha:
            rejections += 1
    
    # Calculate rejection rate
    rejection_rate = rejections / n_samples
    
    return rejection_rate
