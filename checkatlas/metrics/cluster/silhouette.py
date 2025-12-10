from sklearn.metrics import silhouette_score


def run(count_repr, annotations):
    """
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
   
    """
    return silhouette_score(count_repr, annotations)
