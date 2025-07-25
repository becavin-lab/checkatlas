from sklearn.metrics import f1_score


def run(annotation, ref_annotation):
    """

    `Isolated f1 score readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/f1_score/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """

    return f1_score(annotation, ref_annotation, average="weighted")
