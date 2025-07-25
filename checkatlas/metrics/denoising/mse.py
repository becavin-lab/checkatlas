from sklearn.metrics import mean_squared_error


def run(annotation, ref_annotation):
    """
    DENOISING

    `Mean Squared Error score readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/denoising/mse/>`__


    :param annotation:
    :param ref_annotation:
    :return:
    """
    return mean_squared_error(annotation, ref_annotation)
