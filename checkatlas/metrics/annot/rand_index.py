from sklearn.metrics import rand_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `Rand Index readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/rand/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). Rand Index is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Rand Index...")
    return rand_score(annotation, ref_annotation)
