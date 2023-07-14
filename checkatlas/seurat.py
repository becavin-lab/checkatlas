import argparse
import logging
import os
import re
import warnings

import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.rinterface_lib.sexp import NULLType
from rpy2.robjects import pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FactorVector, StrVector

from checkatlas import atlas, checkatlas
from checkatlas.metrics import metrics
from checkatlas.utils import folders

"""
Module for management of Atlas n Seurat format

- Adata and summary tables are created from seuratà
- QC are created using Seurat
- UMAP, t-SNE are created using Seurat
- For creating metrics tables, obs tables are
  first converted from R to python and then
  scikitlearn is used.
"""

logger = logging.getLogger("checkatlas")
warnings.filterwarnings(action="ignore", message="R[write to console]")

SEURAT_TYPE = "Seurat"
SEURAT_EXTENSION = ".rds"


SCANPY_TO_SEURAT_OBS = {
    "total_counts": "nCount_RNA",
    "n_genes_by_counts": "nFeature_RNA",
    "pct_counts_mt": "percent.mito",
    "pct_counts_ribo": "percent.ribo",
    "violin_plot": "violin_plot",
}

SEURAT_TO_SCANPY_OBS = {
    "nCount_RNA": "total_counts",
    "nFeature_RNA": "n_genes_by_counts",
    "percent.mito": "pct_counts_mt",
    "percent.ribo": "pct_counts_ribo",
}


def detect_seurat(atlas_path: str) -> dict:
    if atlas_path.endswith(SEURAT_EXTENSION):
        atlas_info = dict()
        atlas_info[checkatlas.ATLAS_NAME_KEY] = os.path.splitext(
            os.path.basename(atlas_path)
        )[0]
        atlas_info[checkatlas.ATLAS_TYPE_KEY] = SEURAT_TYPE
        atlas_info[checkatlas.ATLAS_EXTENSION_KEY] = SEURAT_EXTENSION
        atlas_info[checkatlas.ATLAS_PATH_KEY] = atlas_path
        return atlas_info
    else:
        return dict()


def check_seurat_install() -> None:
    """Check if Seurat is installed, run installation if not"""
    # import R's utility package
    utils = rpackages.importr("utils")
    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
    # R package names
    packnames = ("Seurat", "SeuratObject")
    # Selectively install what needs to be install.
    # We are fancy, just because we can.
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        # create personal library
        rcode = """dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)"""
        robjects.r(rcode)
        # add to the path
        rcode = """.libPaths(Sys.getenv("R_LIBS_USER"))"""
        robjects.r(rcode)
        logger.debug(f"Set Rlibpaths: {robjects.r(rcode)}")
        utils.install_packages(StrVector(names_to_install))


def read_atlas(atlas_info: dict) -> RS4:
    """Read Seurat object in python using rpy2

    Args:
        atlas_path (str): _description_

    Returns:
        RS4: _description_
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    atlas_path = atlas_info[checkatlas.ATLAS_PATH_KEY]
    logger.info(f"Load {atlas_name} in " f"{atlas_path}")
    rcode = f'readRDS("{atlas_path}")'
    seurat = robjects.r(rcode)
    rclass = robjects.r["class"]
    if rclass(seurat)[0] == "Seurat":
        importr("Seurat")
        return seurat
    else:
        logger.info(f"{atlas_name} is not a Seurat object")
        return None


def get_viable_obs_qc(seurat: RS4, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_QC values
    Extract sorted obs_keys in same order then OBS_QC

    Args:
        seurat (RS4): _description_
        args (argparse.Namespace): _description_

    Returns:
        list: _description_
    """
    r_obs = robjects.r(
        "obs <- function(seurat){ return(colnames(seurat@meta.data))}"
    )
    obs_keys = list()
    for obs_qc in args.qc_display:
        obs_qc = SCANPY_TO_SEURAT_OBS[obs_qc]
        if obs_qc in r_obs(seurat):
            obs_keys.append(obs_qc)
    return obs_keys


def get_viable_obs_annot(seurat: RS4, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_CLUSTERS values
    ! Remove obs_key with only one category !
    Extract sorted obs_keys in same order then OBS_CLUSTERS

    Args:
        seurat (RS4): _description_
        args (argparse.Namespace): _description_

    Returns:
        list: _description_
    """
    obs_keys = list()
    r_obs = robjects.r(
        "obs <- function(seurat){ return(colnames(seurat@meta.data))}"
    )
    obs_key_seurat = r_obs(seurat)
    r_annot = robjects.r(
        "type <- function(seurat, obs_key){ "
        "return(seurat[[obs_key]][[obs_key]])}"
    )
    # Get keys from OBS_CLUSTERS
    for obs_key in obs_key_seurat:
        for obs_key_celltype in args.obs_cluster:
            if obs_key_celltype in obs_key:
                if isinstance(r_annot(seurat, obs_key), FactorVector):
                    obs_keys.append(obs_key)
    # Remove keys with only one category
    obs_keys_final = list()
    for obs_key in obs_keys:
        annotations = r_annot(seurat, obs_key)
        if len(annotations.levels) != 1:
            logger.debug(
                f"Add obs_key {obs_key} with cat {annotations.levels}"
            )
            obs_keys_final.append(obs_key)
    return sorted(obs_keys_final)


def get_viable_obsm(seurat, args):
    """
    Search viable obsm for dimensionality reduction metric
    calc.
    ! No filter on osbm is appled for now !
    :param seurat:
    :param args:
    :return:
    """
    obsm_keys = list()
    # for obsm_key in adata.obsm_keys():
    #   if obsm_key in args.obsm_dimred:
    r_obsm = robjects.r(
        "f<-function(seurat){return(names(seurat@reductions))}"
    )
    obsm_keys_r = r_obsm(seurat)
    obsm_keys = list()
    for obsm_key in obsm_keys_r:
        print(obsm_key)
        obsm_keys.append(obsm_key)
    logger.debug(f"Add obsm {obsm_keys}")
    return obsm_keys


def create_summary_table(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Create a table with all interesting variables
    :param seurat:
    :param atlas_name:
    :param csv_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    logger.debug(f"Create Summary table for {atlas_name}")
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.SUMMARY),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    # Create summary table
    header = [
        "AtlasFileType",
        "NbCells",
        "NbGenes",
        "AnnData.raw",
        "AnnData.X",
        "File_extension",
        "File_path",
    ]
    r_nrow = robjects.r["nrow"]
    r_ncol = robjects.r["ncol"]
    ncells = r_ncol(seurat)[0]
    ngenes = r_nrow(seurat)[0]
    x_raw = False
    x_norm = True
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["AtlasFileType"][atlas_name] = atlas_info[
        checkatlas.ATLAS_TYPE_KEY
    ]
    df_summary["NbCells"][atlas_name] = ncells
    df_summary["NbGenes"][atlas_name] = ngenes
    df_summary["AnnData.raw"][atlas_name] = x_raw
    df_summary["AnnData.X"][atlas_name] = x_norm
    df_summary["File_extension"][atlas_name] = atlas_info[
        checkatlas.ATLAS_EXTENSION_KEY
    ]
    df_summary["File_path"][atlas_name] = atlas_info[checkatlas.ATLAS_PATH_KEY]
    df_summary.to_csv(csv_path, index=False, sep="\t")


def create_anndata_table(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Create a table with all AnnData-like arguments in Seurat object
    :param seurat:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    logger.debug(f"Create Adata table for {atlas_name}")
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.ANNDATA),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    # Create AnnData table
    header = ["atlas_obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)

    # Create r_functions
    r_obs = robjects.r(
        "obs <- function(seurat){ return(colnames(seurat@meta.data))}"
    )
    r_obsm = robjects.r(
        "f<-function(seurat){return(names(seurat@reductions))}"
    )
    r_uns = robjects.r(
        "uns <- function(seurat){ return(colnames(seurat@misc))}"
    )

    obs_list = r_obs(seurat)
    obsm_list = r_obsm(seurat)
    var_list = [""]
    varm_list = [""]
    uns_list = [""]
    if not isinstance(r_uns(seurat), NULLType):
        uns_list = r_uns(seurat)

    df_summary["atlas_obs"][atlas_name] = (
        "<code>" + "</code><br><code>".join(obs_list) + "</code>"
    )
    df_summary["obsm"][atlas_name] = (
        "<code>" + "</code><br><code>".join(obsm_list) + "</code>"
    )
    df_summary["var"][atlas_name] = (
        "<code>" + "</code><br><code>".join(var_list) + "</code>"
    )
    df_summary["varm"][atlas_name] = (
        "<code>" + "</code><br><code>".join(varm_list) + "</code>"
    )
    df_summary["uns"][atlas_name] = (
        "<code>" + "</code><br><code>".join(uns_list) + "</code>"
    )
    df_summary.to_csv(csv_path, index=False, quoting=False, sep="\t")


def create_qc_tables(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Display the atlas QC of seurat
    Search for the metadata variable which correspond
    to the total_RNA, total_UMI, MT_ratio, RT_ratio
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    qc_path = os.path.join(
        folders.get_folder(args.path, folders.QC),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    logger.debug(f"Create QC tables for {atlas_name}")
    obs_keys = get_viable_obs_qc(seurat, args)
    r_meta = robjects.r("obs <- function(seurat){ return(seurat@meta.data)}")
    r_metadata = r_meta(seurat)
    with (ro.default_converter + pandas2ri.converter).context():
        df_metadata = ro.conversion.get_conversion().rpy2py(r_metadata)
        df_annot = df_metadata[obs_keys]
        # rename columns with scanpy names
        new_columns = list()
        for column in df_annot.columns:
            new_columns.append(SEURAT_TO_SCANPY_OBS[column])
        df_annot.columns = new_columns

        # Rank cell by qc metric
        for header in df_annot.columns:
            if header != atlas.CELLINDEX_HEADER:
                new_header = f"cellrank_{header}"
                df_annot = df_annot.sort_values(header, ascending=False)
                df_annot.loc[:, [new_header]] = range(1, len(df_annot) + 1)

        # Sample QC table when more cells than args.plot_celllimit are present
        df_annot = atlas.atlas_sampling(df_annot, "QC", args)
        df_annot.loc[:, [atlas.CELLINDEX_HEADER]] = range(1, len(df_annot) + 1)
        df_annot.to_csv(qc_path, index=False, quoting=False, sep="\t")


def create_qc_plots(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Display the atlas QC
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    qc_path = os.path.join(
        folders.get_folder(args.path, folders.QC_FIG),
        atlas_name + checkatlas.QC_FIG_EXTENSION,
    )
    logger.debug(f"Create QC violin plot for {atlas_name}")
    importr("ggplot2")
    r_cmd = (
        "vln_plot <- function(seurat, obs, qc_path){"
        "vln <- VlnPlot(seurat, features = obs, ncol = length(obs));"
        "ggsave(qc_path, vln, width = 10, "
        "height = 4, dpi = 150)}"
    )
    r_violin = robjects.r(r_cmd)
    obs_keys = list(SEURAT_TO_SCANPY_OBS.keys())
    r_obs = robjects.StrVector(obs_keys)
    r_violin(seurat, r_obs, qc_path)


def create_umap_fig(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    # Search if tsne reduction exists
    r = re.compile(".*umap*.")
    r_names = robjects.r["names"]
    obsm_list = r_names(seurat)
    importr("ggplot2")
    if len(list(filter(r.match, obsm_list))) > 0:
        logger.debug(f"Create UMAP figure for {atlas_name}")
        # Setting up figures directory
        umap_path = os.path.join(
            folders.get_folder(args.path, folders.UMAP),
            atlas_name + checkatlas.UMAP_EXTENSION,
        )
        # Exporting umap
        obs_keys = get_viable_obs_annot(seurat, args)
        r_cmd = (
            "umap <- function(seurat, obs_key, umap_path){"
            "umap_plot <- DimPlot(seurat, group.by = obs_key, "
            'reduction = "umap");'
            "ggsave(umap_path, umap_plot, width = 10, "
            "height = 6, dpi = 76)}"
        )
        r_umap = robjects.r(r_cmd)
        r_umap(seurat, obs_keys[0], umap_path)


def create_tsne_fig(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Display the TSNE of celltypes
    Search for the OBS variable which correspond to the celltype annotation
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    # Search if tsne reduction exists
    r = re.compile(".*tsne*.")
    r_names = robjects.r["names"]
    obsm_list = r_names(seurat)
    importr("ggplot2")
    if len(list(filter(r.match, obsm_list))) > 0:
        logger.debug(f"Create t-SNE figure for {atlas_name}")
        # Setting up figures directory
        tsne_path = os.path.join(
            folders.get_folder(args.path, folders.TSNE),
            atlas_name + checkatlas.TSNE_EXTENSION,
        )
        # Exporting tsne
        obs_keys = get_viable_obs_annot(seurat, args)
        r_cmd = (
            "tsne <- function(seurat, obs_key, tsne_path){"
            "tsne_plot <- DimPlot(seurat, group.by = obs_key, "
            'reduction = "tsne");'
            "ggsave(tsne_path, tsne_plot, width = 10, "
            "height = 6, dpi = 76)}"
        )
        r_tsne = robjects.r(r_cmd)
        r_tsne(seurat, obs_keys[0], tsne_path)


def create_metric_cluster(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Calc clustering metrics
    :param seurat:
    :param atlas_path:
    :param atlas_info:
    :param args:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.CLUSTER),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    header = ["Clust_Sample", "obs"] + args.metric_cluster
    df_cluster = pd.DataFrame(columns=header)
    obs_keys = get_viable_obs_annot(seurat, args)
    obsm_key_representation = "umap"
    if len(obs_keys) > 0:
        logger.debug(f"Calc clustering metrics for {atlas_name}")
        for obs_key in obs_keys:
            dict_line = {
                "Clust_Sample": [atlas_name + "_" + obs_key],
                "obs": [obs_key],
            }
            for metric in args.metric_cluster:
                logger.debug(
                    f"Calc {metric} for {atlas_name} "
                    f"with obs {obs_key} and obsm {obsm_key_representation}"
                )
                metric_value = metrics.calc_metric_cluster_seurat(
                    metric, seurat, obs_key, obsm_key_representation
                )
                dict_line[metric] = metric_value
            df_line = pd.DataFrame(dict_line)
            df_cluster = pd.concat(
                [df_cluster, df_line], ignore_index=True, axis=0
            )
        df_cluster.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obs_key was found for {atlas_name}")


def create_metric_annot(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Calc annotation metrics
    :param adata:
    :param atlas_path:
    :param atlas_info:
    :param args:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.ANNOTATION),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    header = ["Annot_Sample", "Reference", "obs"] + args.metric_annot
    df_annot = pd.DataFrame(columns=header)
    obs_keys = get_viable_obs_annot(seurat, args)
    if len(obs_keys) > 1:
        logger.debug(f"Calc annotation metrics for {atlas_name}")
        if len(obs_keys) != 0:
            ref_obs = obs_keys[0]
            for i in range(1, len(obs_keys)):
                obs_key = obs_keys[i]
                dict_line = {
                    "Annot_Sample": [atlas_name + "_" + obs_key],
                    "Reference": [ref_obs],
                    "obs": [obs_key],
                }
                for metric in args.metric_annot:
                    logger.debug(
                        f"Calc {metric} for {atlas_name} "
                        f"with obs {obs_key} vs ref_obs {ref_obs}"
                    )
                    metric_value = metrics.calc_metric_annot_seurat(
                        metric, seurat, obs_key, ref_obs
                    )
                    dict_line[metric] = metric_value
                df_line = pd.DataFrame(dict_line)
                df_annot = pd.concat(
                    [df_annot, df_line], ignore_index=True, axis=0
                )
            df_annot.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obs_key was found for {atlas_name}")


def create_metric_dimred(
    seurat: RS4, atlas_info: dict, args=argparse.Namespace
) -> None:
    """
    Calc dimensionality reduction metrics
    :param adata:
    :param atlas_path:
    :param atlas_info:
    :param args:
    :return:
    """
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.DIMRED),
        atlas_name + checkatlas.TSV_EXTENSION,
    )
    header = ["Dimred_Sample", "obsm"] + args.metric_dimred
    df_dimred = pd.DataFrame(columns=header)
    # r_reduction = robjects.r(
    #     "reduc <- function(seurat, obsm_key){"
    #     " return(Embeddings(object = seurat, reduction = obsm_key))}"
    # )
    obsm_keys = get_viable_obsm(seurat, args)
    if len(obsm_keys) > 0:
        logger.debug(f"Calc dim red metrics for {atlas_name}")
        for obsm_key in obsm_keys:
            dict_line = {
                "Dimred_Sample": [atlas_name + "_" + obsm_key],
                "obsm": [obsm_key],
            }
            for metric in args.metric_dimred:
                logger.debug(
                    f"Calc {metric} for {atlas_name} with obsm {obsm_key}"
                )
                # r_countmatrix = robjects.r(
                #     "mat <- function(seurat)
                #     { return(seurat@assays$RNA@counts)}"
                # )
                # high_dim_counts = ro.conversion.rpy2py(r_countmatrix(seurat))
                # low_dim_counts = ro.conversion.rpy2py(
                #    r_reduction(seurat, obsm_key)
                # )
                # metric_value = metrics.calc_metric_dimred(
                # metric, high_dim_counts, low_dim_counts)
                logger.warning(
                    "!!! Dim reduction metrics not available for Seurat"
                    " at the moment !!!"
                )
                # metric_value = -1
                # dict_line[metric] = str(metric_value)
            df_line = pd.DataFrame(dict_line)
            df_dimred = pd.concat(
                [df_dimred, df_line], ignore_index=True, axis=0
            )
        df_dimred.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obsm_key was found for {atlas_name}")
