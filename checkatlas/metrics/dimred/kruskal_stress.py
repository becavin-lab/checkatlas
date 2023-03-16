import numpy as np
from sklearn.metrics.pairwise import euclidean_distances


def run(high_dim_counts, low_dim_counts):
    """
    Computes the kruskal stress for the low dimension representation
    stored at adata.obsm.key_repr

    :param adata:
    :param key_repr:
    :return:
    """
    low_dim_distances = euclidean_distances(low_dim_counts, low_dim_counts)
    high_dim_distances = euclidean_distances(high_dim_counts, high_dim_counts)
    stress = np.sqrt(
        np.square(high_dim_distances - low_dim_distances).sum()
        / np.square(low_dim_distances).sum()
    )
    return stress
