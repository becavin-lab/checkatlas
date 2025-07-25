from sklearn.metrics import v_measure_score


def run(annotation, ref_annotation):
    """

    `V-measure readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/v_measure/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """
    return v_measure_score(annotation, ref_annotation)
