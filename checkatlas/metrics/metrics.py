import logging
from sklearn.metrics import (
    adjusted_rand_score,
    davies_bouldin_score,
    silhouette_score,
)
from sklearn.preprocessing import LabelEncoder

from .dim_red import dr_compute

METRICS_CLUST = ["silhouette", "davies_bouldin"]

METRICS_ANNOT = ["rand_index"]

METRICS_DIMRED = ["kruskal_stress"]

logger = logging.getLogger("checkatlas")


def calc_metric_cluster(metric, count_representation, annotation):
    if metric == "silhouette":
        return silhouette_score(count_representation, annotation)
    elif metric == "davies_bouldin":
        return davies_bouldin_score(count_representation, annotation)
    else:
        logger.warning(f"{metric} is not a recognized " f"cluster metric.")
        return -1


def calc_metric_annot(metric, annotation, ref_annotation):
    if metric == "rand_index":
        return adjusted_rand_score(
            *annotation_to_num(annotation, ref_annotation)
        )
    else:
        logger.warning(f"{metric} is not a recognized " f"annotation metric.")
        return -1


def calc_metric_dimred(metric, high_dim_counts, low_dim_counts):
    if metric == "kruskal_stress":
        return dr_compute.kruskal_stress(high_dim_counts, low_dim_counts)
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"dimensionality reduction metric."
        )
        return -1


def annotation_to_num(annotation, ref_annotation):
    """
    Transforms the annotations from categorical to numerical

    Parameters
    ----------
    adata
    partition_key
    reference

    Returns
    -------

    """
    annotation = annotation.to_numpy()
    ref_annotation = ref_annotation.to_numpy()
    le = LabelEncoder()
    le.fit(annotation)
    annotation = le.transform(annotation)
    le2 = LabelEncoder()
    le2.fit(ref_annotation)
    ref_annotation = le2.transform(ref_annotation)
    return annotation, ref_annotation
