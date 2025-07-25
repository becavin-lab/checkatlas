from sklearn.metrics import mutual_info_score


def run(annotation, ref_annotation):
    """

    `Mutual Information readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/mutual_information/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """

    return mutual_info_score(annotation, ref_annotation)
