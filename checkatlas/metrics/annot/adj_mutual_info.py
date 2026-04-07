from sklearn.metrics import adjusted_mutual_info_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `AMI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/adj_mutual_info/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). AMI is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Adjusted Mutual Information...")
    return adjusted_mutual_info_score(annotation, ref_annotation)
