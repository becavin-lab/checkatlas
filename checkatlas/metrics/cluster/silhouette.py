from sklearn.metrics import silhouette_score


def run(count_repr, annotations):
    """

    :param count_repr:
    :param annotations:
    :return:
    """
    return silhouette_score(count_repr, annotations)
