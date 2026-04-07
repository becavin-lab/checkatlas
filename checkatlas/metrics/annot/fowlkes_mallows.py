from sklearn.metrics import fowlkes_mallows_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `Fowlkes-Mallows readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/fowlkes_mallows/>`__

    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). FM is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Fowlkes-Mallows Score...")
    return fowlkes_mallows_score(annotation, ref_annotation)
