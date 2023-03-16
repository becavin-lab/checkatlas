from sklearn.metrics import adjusted_rand_score


def run(annotation, ref_annotation):
    """



    :param annotation:
    :param ref_annotation:
    :return:
    """
    return adjusted_rand_score(annotation, ref_annotation)
