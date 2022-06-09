import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix


def euc(a, b):
    return np.sum(np.square(a - b))


# euclidian(np.array([1,2,3]),np.array([2,3,4]))


def kruskal_stress(adata, key_repr):
    """
    Computes the kruskal stress for the low dimension representation
    stored at adata.obsm.key_repr

    Parameters
    ----------
    adata
    key_repr

    Returns
    -------

    """
    for key in key_repr:
        print(key)
        high_dim_counts = adata.X
        low_dim_counts = adata.obsm[key]
        high_dim_distances = distance_matrix(high_dim_counts, high_dim_counts)
        low_dim_distances = distance_matrix(low_dim_counts, low_dim_counts)
        stress = np.sqrt(
            np.square(high_dim_distances - low_dim_distances).sum()
            / np.square(low_dim_distances).sum()
        )
        print(
            f"dist_diff = "
            f"{np.square(high_dim_distances - low_dim_distances).sum()}"
        )
        print(f"lowdimsum = {np.square(low_dim_distances).sum()}")
        df = pd.DataFrame(low_dim_distances)
        print(df.describe())
        print(stress)
    return stress


# kruskal_stress(test_set, key_repr=['X_pca','X_umap','X_tsne'])


def spearmans_rho(adata, key_repr):
    return
