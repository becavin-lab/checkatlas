from sklearn.metrics import davies_bouldin_score


def run(count_repr, annotations):
    """

    `Davies-Bouldin readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/dbi/>`__

    :param count_repr:
    :param annotations:
    :return:
    """
    # print(count_repr, annotations)
    return davies_bouldin_score(count_repr, annotations)
