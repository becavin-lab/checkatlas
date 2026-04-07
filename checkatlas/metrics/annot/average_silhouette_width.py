from sklearn.metrics import silhouette_score
from scipy.sparse import issparse
import numpy as np


def run(X, labels, n_jobs=-1, verbose=True):
    """
    Calculate the Average Silhouette Width for clustering quality evaluation.
    
    The Average Silhouette Width (Silhouette Score) measures how similar an object is
    to its own cluster compared to other clusters. It combines both cohesion and 
    separation aspects of clustering.
    
    `Average Silhouette Width readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/silhouette/>`__

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points
    :param labels: array-like of shape (n_samples,)
        Cluster labels for each sample
    :param n_jobs: int, default=-1
        Number of parallel jobs for distance computation. -1 uses all cores.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return: float
        The mean Silhouette Coefficient over all samples.
        Range: [-1, 1], where higher values indicate better clustering.
        
    Notes
    -----
    The Silhouette Coefficient for a single sample is defined as:
    
    .. math::
        s = \\frac{b - a}{\\max(a, b)}
    
    where:
    - a: mean intra-cluster distance (cohesion)
    - b: mean nearest-cluster distance (separation)
    
    The average over all samples gives the Average Silhouette Width.
    
    Interpretation:
    - Values close to +1: Sample is well matched to its own cluster
    - Values close to 0: Sample is on or very close to decision boundary
    - Values close to -1: Sample might be assigned to wrong cluster
    """
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    if verbose:
        n_clusters = len(np.unique(labels))
        print(f"Computing Average Silhouette Width ({len(X):,} samples, "
              f"{n_clusters} clusters, n_jobs={n_jobs})...")
    
    score = silhouette_score(X, labels, metric='euclidean', n_jobs=n_jobs)
    
    if verbose:
        print(f"  ASW = {score:.6f}")
    
    return score
