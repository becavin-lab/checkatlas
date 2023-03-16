from sklearn.metrics import adjusted_mutual_info_score


def run(annotation, ref_annotation):
    """



    :param annotation:
    :param ref_annotation:
    :return:
    """
    return adjusted_mutual_info_score(annotation, ref_annotation)
