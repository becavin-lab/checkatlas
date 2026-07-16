import numpy as np
from sklearn.utils.parallel import Parallel, delayed
from scipy.sparse import issparse
from scipy.stats import ks_2samp
from tqdm import tqdm


def run(X, labels, n_jobs=-1, verbose=True):
    """
    Calculate the Kolmogorov-Smirnov (KS) Statistic for clustering separation.

    Evaluates the degree of separation between clusters by comparing
    the distributions of features across different clusters. For each pair of clusters,
    it computes the KS statistic for each feature in parallel and averages them.

    A higher KS score indicates that the clusters have distinct feature distributions.

    Range: [0, 1], where 1 indicates completely distinct distributions (perfect separation).

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix.
    :param labels: array-like of shape (n_samples,)
        Cluster labels.
    :param n_jobs: int, optional (default=-1)
        Number of parallel jobs for feature-wise KS computation. -1 uses all processors.
    :param verbose: bool, optional (default=True)
        Whether to show progress.
    :return: float
        Average KS statistic over all cluster pairs.
    """
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    n_features = X.shape[1]

    # Build cluster pairs
    pairs = []
    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            pairs.append((unique_labels[i], unique_labels[j]))

    if verbose:
        print(
            f"Computing KS statistic for {len(pairs)} cluster pairs × {n_features} features..."
        )

    def _ks_for_pair(label_i, label_j):
        """Compute average KS statistic across all features for one cluster pair."""
        X_i = X[labels == label_i]
        X_j = X[labels == label_j]

        # Vectorized: compute KS for all features at once
        ks_stats = np.zeros(n_features)
        for k in range(n_features):
            stat, _ = ks_2samp(X_i[:, k], X_j[:, k])
            ks_stats[k] = stat

        return np.mean(ks_stats)

    # Parallelize across cluster pairs
    ks_scores = Parallel(n_jobs=n_jobs)(
        delayed(_ks_for_pair)(li, lj)
        for li, lj in tqdm(pairs, desc="KS test (cluster pairs)", disable=not verbose)
    )

    if not ks_scores:
        return 0.0

    final_score = np.mean(ks_scores)

    if verbose:
        print(f"KS Score: {final_score:.4f}")
    return final_score
