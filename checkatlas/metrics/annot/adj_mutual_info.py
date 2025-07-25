from sklearn.metrics import adjusted_mutual_info_score


def run(annotation, ref_annotation):
    """

    `AMI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/adj_mutual_info/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """
    return adjusted_mutual_info_score(annotation, ref_annotation)
