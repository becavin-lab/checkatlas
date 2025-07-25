from sklearn.metrics.cluster import adjusted_rand_score


def run(annotation, ref_annotation):
    """

    `ARI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/adjusted_rand_index/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """
    return adjusted_rand_score(annotation, ref_annotation)
