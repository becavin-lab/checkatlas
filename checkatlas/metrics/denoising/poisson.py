from sklearn.metrics import mean_poisson_deviance


def run(annotation, ref_annotation):
    """
    DENOISING

    `Poisson readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/denoising/poisson/>`__

    :param annotation:
    :param ref_annotation:
    :return:
    """
    return mean_poisson_deviance(ref_annotation, annotation)
