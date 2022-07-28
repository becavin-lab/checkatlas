import logging
import os
import re
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.sexp import NULLType
from rpy2.robjects.vectors import FactorVector

from . import checkatlas, folders


"""
Module for management of Atlas n Seurat format

- Adata and sumarry tables are created from seuratà
- QC are created using Seurat
- UMAP, t-SNE are created using Seurat
- For creating metrics tables, obs tables are first converted from R to python and then scikitlearn is used.


"""

logger = logging.getLogger("checkatlas")

SCANPY_TO_SEURAT_OBS = {"total_counts": "nCount_RNA",
                            "n_genes_by_counts": "nFeature_RNA",
                            "pct_counts_mt": "percent.mito",
                            "pct_counts_ribo": "percent.ribo",
                            "violin_plot": "violin_plot"}

SEURAT_TO_SCANPY_OBS = {"nCount_RNA": "total_counts",
                        "nFeature_RNA": "n_genes_by_counts",
                        "percent.mito": "pct_counts_mt",
                        "percent.ribo": "pct_counts_ribo"}


def read_atlas(atlas_path, atlas_info):
    """
    Read Seurat object in python using rpy2
    :param atlas_path:
    :param atlas_info:
    :return:
    """
    importr("Seurat")
    importr("SeuratObject")
    logger.info(f"Load {atlas_info[0]} in {atlas_info[-1]}")
    rcode = f'readRDS("{atlas_path}")'
    seurat = robjects.r(rcode)
    rclass = robjects.r['class']
    if rclass(seurat)[0] == "Seurat":
        r_seurat = importr('Seurat')
        return seurat
    else:
        logger.info(f"{atlas_info[0]} is not a Seurat object")
        return None


def get_viable_obs_qc(seurat, args):
    """
    Search in obs_keys a match to OBS_QC values
    Extract sorted obs_keys in same order then OBS_QC
    :param adata:
    :return:
    """
    r_obs = robjects.r('obs <- function(seurat){ return(colnames(seurat@meta.data))}')
    obs_keys = list()
    for obs_qc in args.qc_display:
        obs_qc = SCANPY_TO_SEURAT_OBS[obs_qc]
        if obs_qc in r_obs(seurat):
            obs_keys.append(obs_qc)
    return obs_keys


def get_viable_obs_annot(seurat, args):
    """
    Search in obs_keys a match to OBS_CLUSTERS values
    ! Remove obs_key with only one category !
    Extract sorted obs_keys in same order then OBS_CLUSTERS
    :param adata:
    :return:
    """
    obs_keys = list()
    r_obs = robjects.r('obs <- function(seurat){ return(colnames(seurat@meta.data))}')
    obs_key_seurat = r_obs(seurat)
    r_annot = robjects.r('type <- function(seurat, obs_key){ return(seurat[[obs_key]][[obs_key]])}')
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
            logger.debug(f"Add obs_key {obs_key} with cat {annotations.levels}")
            obs_keys_final.append(obs_key)
    return sorted(obs_keys_final)


def create_summary_table(seurat, atlas_path, atlas_info, args) -> None:
    """
    Create a table with all interesting variables
    :param seurat:
    :param atlas_name:
    :param csv_path:
    :return:
    """
    atlas_name = atlas_info[0]
    logger.debug(f"Create Summary table for {atlas_name}")
    atlas_file_type = atlas_info[1]
    atlas_extension = atlas_info[2]
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.SUMMARY),
        atlas_name + checkatlas.SUMMARY_EXTENSION,
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
    r_nrow = robjects.r['nrow']
    r_ncol = robjects.r['ncol']
    ncells = r_ncol(seurat)[0]
    ngenes = r_nrow(seurat)[0]
    x_raw = False
    x_norm = True
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["AtlasFileType"][atlas_name] = atlas_file_type
    df_summary["NbCells"][atlas_name] = ncells
    df_summary["NbGenes"][atlas_name] = ngenes
    df_summary["AnnData.raw"][atlas_name] = x_raw
    df_summary["AnnData.X"][atlas_name] = x_norm
    df_summary["File_extension"][atlas_name] = atlas_extension
    df_summary["File_path"][atlas_name] = atlas_path.replace(args.path, "")
    df_summary.to_csv(csv_path, index=False, sep="\t")


def create_anndata_table(seurat, atlas_path, atlas_info, args) -> None:
    """
    Create a table with all AnnData-like arguments in Seurat object
    :param seurat:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[0]
    logger.debug(f"Create Adata table for {atlas_name}")
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.ANNDATA),
        atlas_name + checkatlas.ADATA_EXTENSION,
        )
    # Create AnnData table
    header = ["obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)

    # Create r_functions
    r_obs = robjects.r('obs <- function(seurat){ return(colnames(seurat@meta.data))}')
    r_names = robjects.r['names']
    r_uns = robjects.r('uns <- function(seurat){ return(colnames(seurat@misc))}')
    obs_list = r_obs(seurat)
    obsm_list = r_names(seurat)
    var_list = []
    varm_list = []
    uns_list = []
    if not isinstance(r_uns(seurat), NULLType):
        uns_list = r_uns(seurat)

    df_summary["obs"][atlas_name] = (
            "<code>"
            + "</code><br><code>".join(obs_list)
            + "</code>"
    )
    df_summary["obsm"][atlas_name] = (
            "<code>"
            + "</code><br><code>".join(obsm_list)
            + "</code>"
    )
    df_summary["var"][atlas_name] = (
            "<code>" + "</code><br><code>".join(var_list) + "</code>"
    )
    df_summary["varm"][atlas_name] = (
            "<code>"
            + "</code><br><code>".join(varm_list)
            + "</code>"
    )
    df_summary["uns"][atlas_name] = (
            "<code>" + "</code><br><code>".join(uns_list) + "</code>"
    )
    df_summary.to_csv(csv_path, index=False, quoting=False, sep="\t")


def create_qc_tables(seurat, atlas_path, atlas_info, args) -> None:
    """
    Display the atlas QC of seurat
    Search for the metadata variable which correspond to the total_RNA, total_UMI,
     MT_ratio, RT_ratio
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[0]
    qc_path = os.path.join(
        folders.get_folder(args.path, folders.QC),
        atlas_name + checkatlas.QC_EXTENSION,
        )
    logger.debug(f"Create QC tables for {atlas_name}")
    obs_keys = get_viable_obs_qc(seurat, args)
    r_meta = robjects.r('obs <- function(seurat){ return(seurat@meta.data)}')
    r_metadata = r_meta(seurat)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_metadata = ro.conversion.rpy2py(r_metadata)
        df_annot = df_metadata[obs_keys]
        # rename columns with scanpy names
        new_columns = list()
        for column in df_annot.columns:
            new_columns.append(SEURAT_TO_SCANPY_OBS[column])
        df_annot.columns = new_columns
        df_annot.to_csv(qc_path, index=False, quoting=False, sep="\t")


def create_qc_plots(seurat, atlas_path, atlas_info, args) -> None:
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
    atlas_name = atlas_info[0]
    qc_path = os.path.join(
        folders.get_folder(args.path, folders.QC_FIG),
        atlas_name + checkatlas.QC_FIG_EXTENSION,
        )
    logger.debug(f"Create QC violin plot for {atlas_name}")
    importr('ggplot2')
    r_cmd = 'vln_plot <- function(seurat, obs, qc_path){' \
            'vln <- VlnPlot(seurat, features = obs, ncol = length(obs));' \
            'ggsave(qc_path, vln, width = 10, ' \
            'height = 4, dpi = 150)}'
    r_violin = robjects.r(r_cmd)
    obs_keys = list(SEURAT_TO_SCANPY_OBS.keys())
    r_obs = robjects.StrVector(obs_keys)
    r_violin(seurat, r_obs, qc_path)


def create_umap_fig(seurat, atlas_path, atlas_info, args) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[0]
    # Search if tsne reduction exists
    r = re.compile(".*umap*.")
    r_names = robjects.r['names']
    obsm_list = r_names(seurat)
    importr('ggplot2')
    if len(list(filter(r.match, obsm_list))) > 0:
        logger.debug(f"Create UMAP figure for {atlas_name}")
        # Setting up figures directory
        umap_path = os.path.join(
            folders.get_folder(args.path, folders.UMAP),
            atlas_name + checkatlas.UMAP_EXTENSION,
            )
        # Exporting umap
        obs_keys = get_viable_obs_annot(seurat, args)
        r_cmd = 'umap <- function(seurat, obs_key, umap_path){' \
                'umap_plot <- DimPlot(seurat, group.by = obs_key, ' \
                'reduction = "umap");' \
                'ggsave(umap_path, umap_plot, width = 10, ' \
                'height = 6, dpi = 76)}'
        r_umap = robjects.r(r_cmd)
        r_umap(seurat, obs_keys[0], umap_path)


def create_tsne_fig(seurat, atlas_path, atlas_info, args) -> None:
    """
    Display the TSNE of celltypes
    Search for the OBS variable which correspond to the celltype annotation
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    atlas_name = atlas_info[0]
    # Search if tsne reduction exists
    r = re.compile(".*tsne*.")
    r_names = robjects.r['names']
    obsm_list = r_names(seurat)
    importr('ggplot2')
    if len(list(filter(r.match, obsm_list))) > 0:
        logger.debug(f"Create t-SNE figure for {atlas_name}")
        # Setting up figures directory
        tsne_path = os.path.join(
            folders.get_folder(args.path, folders.TSNE),
            atlas_name + checkatlas.TSNE_EXTENSION,
            )
        # Exporting tsne
        obs_keys = get_viable_obs_annot(seurat, args)
        r_cmd = 'tsne <- function(seurat, obs_key, tsne_path){' \
                'tsne_plot <- DimPlot(seurat, group.by = obs_key, ' \
                'reduction = "tsne");' \
                'ggsave(tsne_path, tsne_plot, width = 10, ' \
                'height = 6, dpi = 76)}'
        r_tsne = robjects.r(r_cmd)
        r_tsne(seurat, obs_keys[0], tsne_path)


def metric_cluster():
    print('yo')
def metric_annot():
    print('yo')
def metric_dimred():
    print('yo')

# def install_library():
#     # R package names
#     packnames = ["Seurat"]
#     utils = rpackages.importr("utils")
#
#     # select a mirror for R packages
#     utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
#
#     # Selectively install what needs to be install.
#     # We are fancy, just because we can.
#     names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
#     if len(names_to_install) > 0:
#         utils.install_packages(StrVector(names_to_install))
