from sklearn.metrics import rand_score


def run(annotation, ref_annotation):
    """

    `Rand Index readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/rand/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """
    return rand_score(annotation, ref_annotation)
