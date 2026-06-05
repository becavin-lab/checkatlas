import numpy as np
from sklearn.metrics import davies_bouldin_score


def run(count_repr, annotations, n_jobs=-1, verbose=True):
    """
    Calculate the Davies-Bouldin Index.

    Measures the average similarity ratio of each cluster with the cluster
    that is most similar to it. Lower values indicate better clustering.

    `Davies-Bouldin readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/dbi/>`__

    :param count_repr: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points.
    :param annotations: array-like of shape (n_samples,)
        Cluster labels for each sample.
    :param n_jobs: int, optional (default=-1)
        Not used directly (sklearn DBI doesn't parallelize), kept for API consistency.
    :param verbose: bool, optional (default=True)
        Whether to print progress.
    :return: float
        The Davies-Bouldin Index. Lower is better (0 is optimal).
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)

    if verbose:
        print("Computing Davies-Bouldin Index...")
    return davies_bouldin_score(count_repr, annotations)
