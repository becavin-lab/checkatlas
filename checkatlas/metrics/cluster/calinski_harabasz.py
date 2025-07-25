from sklearn.metrics import calinski_harabasz_score


def run(count_repr, annotations):
    """
    `Calinski-Harabasz readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/calinski_harabasz/>`__

    :param count_repr:
    :param annotations:
    :return:
    """
    return calinski_harabasz_score(count_repr, annotations)
