from sklearn.metrics import v_measure_score


def run(annotation, ref_annotation):
    """



    :param annotation:
    :param ref_annotation:
    :return:
    """
    return v_measure_score(annotation, ref_annotation)
