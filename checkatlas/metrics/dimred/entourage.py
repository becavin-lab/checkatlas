import numpy as np
from sklearn.neighbors import NearestNeighbors
import logging

logger = logging.getLogger("checkatlas")


def run(high_dim_counts, low_dim_counts):
    """
    `Entourage readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/dimred/entourage/>`__

    :param high_dim_counts:
    :param low_dim_counts:
    :return:
    """
    # for i in range(4):
    k_neighbors = 4

    # high_dim_counts = high_dim_counts[0:100,]
    # low_dim_counts = low_dim_counts[0:100,]

    X = high_dim_counts
    neigh = NearestNeighbors(n_neighbors=k_neighbors)
    neigh.fit(X)
    A = neigh.kneighbors_graph(X, mode="distance")

    X_low = low_dim_counts
    neigh_low = NearestNeighbors(n_neighbors=k_neighbors)
    neigh_low.fit(X_low)
    B = neigh_low.kneighbors_graph(X_low, mode="distance")

    n = X.shape[0]
    k = k_neighbors

    total_inter = 0
    for i in range(n):
        if i % 100 == 0:
            print(f"i{i} n{n}")
        indice_high = A[i].indices.tolist()
        indice_low = B[i].indices.tolist()

        new_ind_high = np.array(indice_high[1:])
        new_ind_low = np.array(indice_low[1:])

        inter = np.intersect1d(new_ind_high, new_ind_low)
        card = len(inter)
        total_inter = total_inter + card

    entourage_score = total_inter / (n * k)
    logger.debug(f"Entourage : {entourage_score}")

    return entourage_score
