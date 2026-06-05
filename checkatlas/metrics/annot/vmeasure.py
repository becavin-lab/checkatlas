from sklearn.metrics import v_measure_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `V-measure readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/v_measure/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). V-measure is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing V-measure...")
    return v_measure_score(annotation, ref_annotation)
