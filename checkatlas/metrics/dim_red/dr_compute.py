import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import euclidean_distances


def kruskal_stress(adata, key_repr):
    """
    Computes the kruskal stress for the low dimension representation
    stored at adata.obsm.key_repr

    :param adata:
    :param key_repr:
    :return:
    """
    high_dim_counts = adata.X
    low_dim_counts = adata.obsm[key_repr]
    low_dim_distances = euclidean_distances(low_dim_counts, low_dim_counts)
    high_dim_distances = euclidean_distances(high_dim_counts, high_dim_counts)
    stress = np.sqrt(
        np.square(high_dim_distances - low_dim_distances).sum()
        / np.square(low_dim_distances).sum()
    )
    return stress


def spearmans_rho(adata, key_repr):
    return
