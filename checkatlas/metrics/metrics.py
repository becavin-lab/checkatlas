import logging

import rpy2.robjects as ro
import rpy2.robjects as robjects
from sklearn.preprocessing import LabelEncoder

from . import annot, cluster, dimred

METRICS_CLUST = cluster.__all__
METRICS_ANNOT = annot.__all__
METRICS_DIMRED = dimred.__all__

logger = logging.getLogger("checkatlas")

R_ANNOT = robjects.r(
    "type <- function(seurat, obs_key){ "
    "return(seurat[[obs_key]][[obs_key]])}"
)
R_REDUCTION = robjects.r(
    "reduc <- function(seurat, obsm_key){"
    " return(Embeddings(object = seurat, "
    "reduction = obsm_key))}"
)


def calc_metric_cluster_scanpy(
    metric, adata, obs_key, obsm_key_representation
):
    if metric in METRICS_CLUST:
        metric_module = getattr(cluster, metric)
        annotations = adata.obs[obs_key]
        if obsm_key_representation != "":
            count_repr = adata.obsm[obsm_key_representation]
            return metric_module.run(count_repr, annotations)
        else:
            original_count = adata.X.toarray()
            return metric_module.run(original_count, annotations)
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_cluster_seurat(
    metric, seurat, obs_key, obsm_key_representation
):
    if metric in METRICS_CLUST:
        metric_module = getattr(cluster, metric)
        annotations = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        count_repr = ro.conversion.rpy2py(
            R_REDUCTION(seurat, obsm_key_representation)
        )
        return metric_module.run(count_repr, annotations)
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_annot_scanpy(metric, adata, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        metric_module = getattr(annot, metric)
        annotation = adata.obs[obs_key]
        ref_annotation = adata.obs[ref_obs]
        annotation, ref_annotation = annotation_to_num(
            annotation, ref_annotation
        )
        return metric_module.run(annotation, ref_annotation)
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def calc_metric_annot_seurat(metric, seurat, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        metric_module = getattr(annot, metric)
        annotation = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        ref_annotation = ro.conversion.rpy2py(R_ANNOT(seurat, ref_obs))
        # annotation, ref_annotation = annotation_to_num(
        #    annotation, ref_annotation
        # )
        return metric_module.run(annotation, ref_annotation)
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def calc_metric_dimred(metric, adata, obsm_key):
    if metric in METRICS_DIMRED:
        metric_module = getattr(dimred, metric)
        high_dim_counts = adata.X
        low_dim_counts = adata.obsm[obsm_key]
        return metric_module.run(high_dim_counts, low_dim_counts)
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"dimensionality reduction metric."
            f"\nList of dim. red. metrics: {METRICS_DIMRED}"
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
