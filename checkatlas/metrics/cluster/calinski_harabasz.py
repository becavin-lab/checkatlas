from sklearn.metrics import calinski_harabasz_score
import numpy as np


def run(count_repr, annotations, n_jobs=-1, verbose=True):
    """
    Calculate the Calinski-Harabasz Index (Variance Ratio Criterion).
    
    Measures the ratio of between-cluster dispersion to within-cluster dispersion.
    Higher values indicate better-defined clusters.

    `Calinski-Harabasz readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/calinski_harabasz/>`__

    :param count_repr: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points.
    :param annotations: array-like of shape (n_samples,)
        Cluster labels for each sample.
    :param n_jobs: int, optional (default=-1)
        Not used directly (O(N) computation), kept for API consistency.
    :param verbose: bool, optional (default=True)
        Whether to print progress.
    :return: float
        The Calinski-Harabasz Index. Higher is better.
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)
    
    if verbose: print("Computing Calinski-Harabasz Index...")
    return calinski_harabasz_score(count_repr, annotations)
