from sklearn.metrics.cluster import adjusted_rand_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `ARI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/adjusted_rand_index/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). ARI is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Adjusted Rand Index...")
    return adjusted_rand_score(annotation, ref_annotation)
