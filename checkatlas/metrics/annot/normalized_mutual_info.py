from sklearn.metrics import normalized_mutual_info_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `NMI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/nmi/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). NMI is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Normalized Mutual Information...")
    return normalized_mutual_info_score(annotation, ref_annotation)
