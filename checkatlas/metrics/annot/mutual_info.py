from sklearn.metrics import mutual_info_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `Mutual Information readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/mutual_information/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). MI is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Mutual Information...")
    return mutual_info_score(annotation, ref_annotation)
