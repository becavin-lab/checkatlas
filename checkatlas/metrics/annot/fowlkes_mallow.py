from sklearn.metrics import fowlkes_mallows_score


def run(annotation, ref_annotation):
    """



    :param annotation:
    :param ref_annotation:
    :return:
    """
    return fowlkes_mallows_score(annotation, ref_annotation)
