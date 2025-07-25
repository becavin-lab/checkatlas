from sklearn.metrics import normalized_mutual_info_score


def run(annotation, ref_annotation):
    """

    `NMI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/nmi/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """

    return normalized_mutual_info_score(annotation, ref_annotation)
