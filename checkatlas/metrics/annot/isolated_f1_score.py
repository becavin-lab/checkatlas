from sklearn.metrics import f1_score


def run(annotation, ref_annotation, n_jobs=-1, verbose=True):
    """

    `Isolated f1 score readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/cellannotation/f1_score/>`__


    :param annotation:
    :param ref_annotation:
    :param n_jobs: int, default=-1
        Not used (API consistency). F1 is O(N) and not parallelizable.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return:
    """
    if verbose:
        print("Computing Isolated F1 Score...")
    return f1_score(annotation, ref_annotation, average="weighted")
