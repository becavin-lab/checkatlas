from sklearn.metrics import fowlkes_mallows_score


def run(annotation, ref_annotation):
    """

    `Fowlkes-Mallows readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/fowlkes_mallows/>`__

    :param annotation:
    :param ref_annotation:
    :return:
    """
    return fowlkes_mallows_score(annotation, ref_annotation)
