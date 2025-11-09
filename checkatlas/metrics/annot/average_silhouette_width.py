from sklearn.metrics import silhouette_score


def run(X, labels):
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
    return silhouette_score(X, labels, metric='euclidean')
