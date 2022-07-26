import logging

from .cluster import clust_compute
from .dim_red import dr_compute

METRICS_CLUST = ["silhouette", "davies_bouldin"]

METRICS_ANNOT = ["rand_index"]

METRICS_DIMRED = ["kruskal_stress"]

logger = logging.getLogger("checkatlas")


def calc_metric_cluster(metric, adata, obs_key, obsm_key):
    if metric == "silhouette":
        return clust_compute.silhouette(adata, obs_key, obsm_key)
    elif metric == "davies_bouldin":
        return clust_compute.davies_bouldin(adata, obs_key, obsm_key)
    else:
        logger.warning(f"{metric} is not a recognized " f"cluster metric.")
        return -1


def calc_metric_annot(metric, adata, obs_key, ref_obs):
    if metric == "rand_index":
        return clust_compute.rand(adata, obs_key, ref_obs)
    else:
        logger.warning(f"{metric} is not a recognized " f"annotation metric.")
        return -1


def calc_metric_dimred(metric, adata, obsm_key):
    if metric == "kruskal_stress":
        return dr_compute.kruskal_stress(adata, obsm_key)
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"dimensionality reduction metric."
        )
        return -1
