from sklearn.metrics import calinski_harabasz_score


def run(count_repr, annotations):
    """

    :param count_repr:
    :param annotations:
    :return:
    """
    return calinski_harabasz_score(count_repr, annotations)
