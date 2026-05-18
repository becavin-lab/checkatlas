import argparse
import logging
import os
import re
import warnings
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata import _io as _io
from sklearn.utils.fixes import _object_dtype_isnan

try:
    from . import cellranger, check
    from .metrics import metrics
    from .utils import files, folders
    from .utils.col_detector import CheckAtlasColumnDetector
except ImportError:
    from checkatlas import cellranger, check
    from checkatlas.metrics import metrics
    from checkatlas.utils import files, folders
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

"""
Atlas module for AnnData (Scanpy)
All the function to screen the atlases
"""

ANNDATA_TYPE = "AnnData"
ANNDATA_EXTENSION = ".h5ad"

CELLINDEX_HEADER = "cell_index"

# Default obs columns to search for clustering/annotation keys
OBS_CLUSTERS = [
    "leiden",
    "louvain",
    "celltype",
    "cell_type",
    "cell_annotation",
    "cluster",
    "seurat_clusters",
    "annotation",
]

# Default obs columns for QC display
OBS_QC = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "n_genes",
    "n_counts",
]

logger = logging.getLogger("checkatlas")

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
sc.settings.verbosity = 0



def detect_scanpy(atlas_path: str) -> dict:
    if atlas_path.endswith(ANNDATA_EXTENSION):
        atlas_info = dict()
        atlas_info[check.ATLAS_NAME_KEY] = os.path.splitext(
            os.path.basename(atlas_path)
        )[0]
        atlas_info[check.ATLAS_TYPE_KEY] = ANNDATA_TYPE
        atlas_info[check.ATLAS_EXTENSION_KEY] = ANNDATA_EXTENSION
        atlas_info[check.ATLAS_PATH_KEY] = atlas_path
        return atlas_info
    else:
        return dict()


def read_atlas(atlas_info: dict) -> AnnData:
    """
    Read Scanpy or Cellranger data : .h5ad or .h5

    Args:
        atlas_info (dict): info dict about the atlas

    Returns:
        AnnData: scanpy object from .h5ad
    """
    logger.info(
        f"Load {atlas_info[check.ATLAS_NAME_KEY]} "
        f"in {atlas_info[check.ATLAS_PATH_KEY]}"
    )
    try:
        if (
            atlas_info[check.ATLAS_TYPE_KEY]
            == cellranger.CELLRANGER_TYPE_CURRENT
        ):
            logger.debug(
                "Read Cellranger >= v3 results "
                f"{atlas_info[check.ATLAS_PATH_KEY]}"
            )
            adata = cellranger.read_cellranger_current(atlas_info)
        elif (
            atlas_info[check.ATLAS_TYPE_KEY]
            == cellranger.CELLRANGER_TYPE_OBSOLETE
        ):
            logger.debug(
                "Read Cellranger < v3 results "
                f"{atlas_info[check.ATLAS_PATH_KEY]}"
            )
            adata = cellranger.read_cellranger_obsolete(atlas_info)
        else:
            logger.debug(
                f"Read Scanpy file {atlas_info[check.ATLAS_PATH_KEY]}"
            )
            adata = sc.read_h5ad(atlas_info[check.ATLAS_PATH_KEY])
        return adata
    except _io.utils.AnnDataReadError:
        logger.warning(
            "AnnDataReadError, cannot read: "
            f"{atlas_info[check.ATLAS_PATH_KEY]}"
        )
        return dict()


def clean_scanpy_atlas(adata: AnnData, atlas_info: dict) -> AnnData:
    """
    Clean the Scanpy object to be sure to get all information out of it

    - Make var names unique
    - Make var unique for Raw matrix
    - If OBS_CLUSTERS are present and in int32 -> be sure to
    transform them in categorical

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas

    Returns:
        AnnData: cleaned atlas
    """
    logger.debug(f"Clean scanpy: {atlas_info[check.ATLAS_NAME_KEY]}")
    # Make var names unique
    list_var = adata.var_names
    if len(set(list_var)) == len(list_var):
        logger.debug("Var names unique")
    else:
        logger.debug(
            "Var names not unique, ran : adata.var_names_make_unique()"
        )
        adata.var_names_make_unique()
        # Test a second time if it is unique (sometimes it helps)
        list_var = adata.var_names
        if len(set(list_var)) == len(list_var):
            logger.debug("Var names unique")
        else:
            logger.debug(
                "Var names not unique, ran : adata.var_names_make_unique()"
            )
            adata.var_names_make_unique()
            # If it is still not unique, create unique var_names "by hand"
            list_var = adata.var_names
            if len(set(list_var)) == len(list_var):
                logger.debug("Var names unique")
            else:
                logger.debug(
                    "Var names not unique, ran : adata.var_names_make_unique()"
                )
                adata.var.index = [
                    x + "_" + str(i)
                    for i, x in zip(range(len(adata.var)), adata.var_names)
                ]
                list_var = adata.var_names
                if len(set(list_var)) == len(list_var):
                    logger.debug("Var names unique")
    # Make var unique for Raw matrix
    if adata.raw is not None:
        list_var = adata.raw.var_names
        if len(set(list_var)) == len(list_var):
            logger.debug("Var names for Raw unique, transform ")
        else:
            logger.debug("Var names for Raw not unique")
            adata.raw.var.index = [
                x + "_" + str(i)
                for i, x in zip(range(len(adata.raw.var)), adata.raw.var_names)
            ]
            list_var = adata.raw.var_names
            if len(set(list_var)) == len(list_var):
                logger.debug("Var names for Raw unique")

    # If OBS_CLUSTERS are present and in int32 -> be sure to
    # transform them in categorical
    for obs_key in adata.obs_keys():
        for obs_key_celltype in OBS_CLUSTERS:
            if obs_key_celltype in obs_key:
                if (
                    adata.obs[obs_key].dtype == np.int32
                    or adata.obs[obs_key].dtype == np.int64
                ):
                    adata.obs[obs_key] = pd.Categorical(adata.obs[obs_key])
    return adata


def get_viable_obs_qc(adata: AnnData, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_QC values
    Extract sorted obs_keys in same order then OBS_QC

    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obs_keys
    """
    obs_keys = list()
    for obs_key in adata.obs_keys():
        if obs_key in args.qc_display:
            obs_keys.append(obs_key)
    return obs_keys


def get_viable_obs_annot(adata: AnnData, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_CLUSTERS values
    ! Remove obs_key with only one category !
    Extract sorted obs_keys in same order then OBS_CLUSTERS

    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obs_keys
    """
    obs_keys = list()
    # Get keys from OBS_CLUSTERS
    for obs_key in adata.obs_keys():
        for obs_key_celltype in args.obs_cluster:
            if obs_key_celltype in obs_key:
                if isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype):
                    obs_keys.append(obs_key)
    # Remove keys with only one category and no NaN in the array
    obs_keys_final = list()
    for obs_key in obs_keys:
        annotations = adata.obs[obs_key]
        if not _object_dtype_isnan(annotations).any():
            categories_temp = annotations.cat.categories
            # remove nan if found
            categories = categories_temp.dropna()
            if True in categories.isin(["nan"]):
                index = categories.get_loc("nan")
                categories = categories.delete(index)
            # Add obs_key with more than one category (with Nan removed)
            if len(categories) != 1:
                logger.debug(
                    f"Add obs_key {obs_key} with cat {categories_temp}"
                )
                obs_keys_final.append(obs_key)
    return sorted(obs_keys_final)


def get_viable_obsm(adata: AnnData, args: argparse.Namespace) -> list:
    """
    TO DO
    Search viable obsm for dimensionality reduction metric
    calc.
    ! No filter on osbm is appled for now !
    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obsm_keys
    """
    obsm_keys = list()
    # for obsm_key in adata.obsm_keys():
    #   if obsm_key in args.obsm_dimred:
    obsm_keys = adata.obsm_keys()
    logger.debug(f"Add obsm {obsm_keys}")
    return obsm_keys


def create_summary_table(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Create a table with all summarizing variables

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (str): info dict of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    atlas_type = atlas_info[check.ATLAS_TYPE_KEY]
    atlas_path = atlas_info[check.ATLAS_PATH_KEY]
    logger.debug(f"Create Summary table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.SUMMARY, check.TSV_EXTENSION, args.path
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
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["AtlasFileType"][atlas_name] = atlas_type
    df_summary["NbCells"][atlas_name] = adata.n_obs
    df_summary["NbGenes"][atlas_name] = adata.n_vars
    df_summary["AnnData.raw"][atlas_name] = adata.raw is not None
    df_summary["AnnData.X"][atlas_name] = adata.X is not None
    df_summary["File_extension"][atlas_name] = atlas_name
    df_summary["File_path"][atlas_name] = atlas_path
    df_summary.to_csv(csv_path, index=False, sep="\t")


def create_anndata_table(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Create an html table with all AnnData arguments
    The html code will make all elements of the table visible in MultiQC
    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info dict on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]

    logger.debug(f"Create Adata table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.ANNDATA, check.TSV_EXTENSION, args.path
    )
    # Create AnnData table
    header = ["atlas_obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    # html_element = "<span class=\"label label-primary\">"
    # new_line = ''
    # for value in list(adata.obs.columns):
    #     new_line += html_element + value + "</span><br>"
    #     print(new_line)
    df_summary["atlas_obs"][atlas_name] = (
        "<code>"
        + "</code><br><code>".join(list(adata.obs.columns))
        + "</code>"
    )
    df_summary["obsm"][atlas_name] = (
        "<code>"
        + "</code><br><code>".join(list(adata.obsm_keys()))
        + "</code>"
    )
    df_summary["var"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.var_keys())) + "</code>"
    )
    df_summary["varm"][atlas_name] = (
        "<code>"
        + "</code><br><code>".join(list(adata.varm_keys()))
        + "</code>"
    )
    df_summary["uns"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.uns_keys())) + "</code>"
    )
    df_summary.to_csv(csv_path, index=False, quoting=False, sep="\t")


def create_qc_tables(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Display the atlas QC table
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    qc_path = files.get_file_path(
        atlas_name, folders.QC, check.TSV_EXTENSION, args.path
    )
    logger.debug(f"Create QC tables for {atlas_name}")
    qc_genes = []
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("mt")
        logger.debug(f"Mitochondrial genes in {atlas_name} for QC")
    else:
        logger.debug(f"No mitochondrial genes in {atlas_name} for QC")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("ribo")
        logger.debug(f"Ribosomal genes in {atlas_name} for QC")
    else:
        logger.debug(f"No ribosomal genes in {atlas_name} for QC")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_genes,
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    df_annot = adata.obs[get_viable_obs_qc(adata, args)]
    # Rank cell by qc metric
    for header in df_annot.columns:
        if header != CELLINDEX_HEADER:
            new_header = f"cellrank_{header}"
            df_annot = df_annot.sort_values(header, ascending=False)
            df_annot.loc[:, [new_header]] = range(1, adata.n_obs + 1)

    # Sample QC table when more cells than args.plot_celllimit are present
    df_annot = atlas_sampling(df_annot, "QC", args)
    df_annot.loc[:, [CELLINDEX_HEADER]] = range(1, len(df_annot) + 1)
    df_annot.to_csv(qc_path, index=False, quoting=False, sep="\t")


def create_qc_plots(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Display the atlas QC plot
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.settings.figdir = folders.get_workingdir(args.path)
    sc.set_figure_params(dpi_save=80)
    qc_path = os.sep + atlas_name + check.QC_FIG_EXTENSION
    logger.debug(f"Create QC violin plot for {atlas_name}")
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    sc.pl.violin(
        adata,
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
        ],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=qc_path,
    )


def create_umap_fig(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.set_figure_params(dpi_save=150)
    # Search if umap reduction exists
    obsm_keys = get_viable_obsm(adata, args)
    r = re.compile(".*umap*.")
    obsm_umap_keys = list(filter(r.match, obsm_keys))
    if len(obsm_umap_keys) > 0:
        obsm_umap = obsm_umap_keys[0]
        logger.debug(
            f"Create UMAP figure for {atlas_name} with obsm={obsm_umap}"
        )
        # Set the umap to display
        if isinstance(adata.obsm[obsm_umap], pd.DataFrame):
            # Transform to numpy if it is a pandas dataframe
            adata.obsm["X_umap"] = adata.obsm[obsm_umap].to_numpy()
        else:
            adata.obsm["X_umap"] = adata.obsm[obsm_umap]
        # Setting up figures directory
        sc.settings.figdir = folders.get_workingdir(args.path)
        umap_path = os.sep + atlas_name + check.UMAP_EXTENSION
        # Exporting umap
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.umap(adata, color=obs_keys[0], show=False, save=umap_path)
        else:
            sc.pl.umap(adata, show=False, save=umap_path)


def create_tsne_fig(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Display the TSNE of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.set_figure_params(dpi_save=150)
    # Search if tsne reduction exists
    obsm_keys = get_viable_obsm(adata, args)
    r = re.compile(".*tsne*.")
    obsm_tsne_keys = list(filter(r.match, obsm_keys))
    if len(obsm_tsne_keys) > 0:
        obsm_tsne = obsm_tsne_keys[0]
        logger.debug(
            f"Create t-SNE figure for {atlas_name} with obsm={obsm_tsne}"
        )
        # Set the t-sne to display
        if isinstance(adata.obsm[obsm_tsne], pd.DataFrame):
            # Transform to numpy if it is a pandas dataframe
            adata.obsm["X_tsne"] = adata.obsm[obsm_tsne].to_numpy()
        else:
            adata.obsm["X_tsne"] = adata.obsm[obsm_tsne]
        # Setting up figures directory
        sc.settings.figdir = sc.settings.figdir = folders.get_workingdir(
            args.path
        )
        tsne_path = os.sep + atlas_name + check.TSNE_EXTENSION
        # Exporting tsne
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.tsne(adata, color=obs_keys[0], show=False, save=tsne_path)
        else:
            sc.pl.tsne(adata, show=False, save=tsne_path)


def create_metric_cluster(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calc clustering metrics

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name,
        folders.CLUSTER,
        check.TSV_EXTENSION,
        args.path,
    )
    header = ["Clust_Sample", "obs"] + args.metric_cluster
    df_cluster = pd.DataFrame(columns=header)
    obs_keys = get_viable_obs_annot(adata, args)
    obsm_keys = get_viable_obsm(adata, args)
    r = re.compile(".*umap*.")
    obsm_umap_keys = list(filter(r.match, obsm_keys))
    r = re.compile(".*tsne*.")
    obsm_tsne_keys = list(filter(r.match, obsm_keys))
    obsm_key_representation = ""
    if len(obsm_umap_keys) > 0:
        obsm_key_representation = obsm_umap_keys[0]
        print("reach", obsm_key_representation)
    elif len(obsm_tsne_keys) > 0:
        obsm_key_representation = obsm_tsne_keys[0]
        print("reach", obsm_key_representation)

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
                metric_value, running_time = (
                    metrics.calc_metric_cluster_scanpy(
                        metric, adata, obs_key, obsm_key_representation
                    )
                )
                dict_line[metric] = metric_value
                dict_line[f"{metric}_running_time"] = running_time
            df_line = pd.DataFrame(dict_line)
            df_cluster = pd.concat(
                [df_cluster, df_line], ignore_index=True, axis=0
            )
        df_cluster.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obs_key was found for {atlas_name}")


def create_metric_annot(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calculate all annotation metrics via the comprehensive ``cal_annot``
    pipeline (ref-vs-pred, embedding-based, batch/integration, graph).

    The pipeline auto-detects reference/predicted columns, embedding keys,
    and batch labels; runs every metric listed in ``METRICS_ANNOT``; and
    writes results as tab-separated files in the annotation folder.

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    annotation_dir = folders.get_folder(args.path, folders.ANNOTATION)

    logger.info("Running full annotation pipeline for %s", atlas_name)

    # cal_annot auto-detects columns, runs specified metrics (ref-vs-pred
    # plus embedding/batch/graph-dependent ones), and saves its own CSV.
    df = metrics.cal_annot(
        adata,
        atlas_name=atlas_name,
        metric_list=args.metric_annot,
        all=True,               # fallback: run every metric in METRICS_ANNOT
        file_dir=annotation_dir,
        n_jobs=-1,
        verbose=True,
    )

    # Also write a MultiQC-compatible TSV in the annotation folder.
    if not df.empty:
        csv_path = files.get_file_path(
            atlas_name,
            folders.ANNOTATION,
            check.TSV_EXTENSION,
            args.path,
        )
        df.to_csv(csv_path, index=False, sep="\t")
        logger.info("Annotation metrics saved to %s", csv_path)
    else:
        logger.warning("No annotation metrics calculated for %s", atlas_name)


def create_metric_dimred(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calc dimensionality reduction metrics

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name,
        folders.DIMRED,
        check.TSV_EXTENSION,
        args.path,
    )
    header = ["Dimred_Sample", "obsm"] + args.metric_dimred
    df_dimred = pd.DataFrame(columns=header)
    obsm_keys = get_viable_obsm(adata, args)
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
                metric_value, running_time = metrics.calc_metric_dimred(
                    metric, adata, obsm_key
                )
                dict_line[metric] = metric_value
                dict_line[f"{metric}_running_time"] = running_time
            df_line = pd.DataFrame(dict_line)
            df_dimred = pd.concat(
                [df_dimred, df_line], ignore_index=True, axis=0
            )
        df_dimred.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obsm_key was found for {atlas_name}")


def atlas_sampling(
    df_annot: pd.DataFrame, type_df: str, args: argparse.Namespace
) -> pd.DataFrame:
    """
    If args.plot_celllimit != 0 and args.plot_celllimit < len(df_annot)
    The atlas qC table will be sampled for MultiQC

    Args:
        df_annot (pd.DataFrame): Table to sample
        type_df (str): type of table
        args (argparse.Namespace): arguments of checkatlas workflow

    Returns:
        pd.DataFrame: Sampled QC table
    """
    if args.plot_celllimit != 0 and args.plot_celllimit < len(df_annot):
        logger.debug(f"Sample {type_df} table with {len(df_annot)} cells")
        df_annot = df_annot.sample(args.plot_celllimit)
        logger.debug(f"{type_df} table sampled to {len(df_annot)} cells")
    return df_annot


# Public API functions for column detection

def col_annotation_ref(adata: AnnData, 
                       min_score: float = 0.5,
                       return_all: bool = False) -> Optional[str]:
    """
    Detect reference (ground truth) annotation column in AnnData object.
    
    This function uses intelligent semantic and statistical analysis to identify
    the most likely reference/ground truth cell type annotation column.
    
    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return list of all candidates with scores. Default: False
        
    Returns:
        str or List[Tuple[str, float]] or None: 
            - If return_all=False: Best reference column name, or None if none found
            - If return_all=True: List of (column_name, score) tuples sorted by score
            
    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> ref_col = atlas.col_annotation_ref(adata)
        >>> print(f"Reference column: {ref_col}")
        >>> 
        >>> # Get all candidates with scores
        >>> all_refs = atlas.col_annotation_ref(adata, return_all=True)
        >>> for col, score in all_refs:
        ...     print(f"{col}: {score:.3f}")
    """
    
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(
        min_reference_score=min_score,
        min_predicted_score=0.3
    )
    
    ref_candidates = results['annotation']['reference']
    
    if return_all:
        return ref_candidates
    else:
        return ref_candidates[0][0] if ref_candidates else None


def col_annotation_pred(adata: AnnData,
                        min_score: float = 0.5,
                        return_all: bool = False,
                        max_results: int = 5) -> Optional[List[str]]:
    """
    Detect predicted/cluster annotation columns in AnnData object.
    
    This function identifies columns containing cluster labels or automated
    cell type predictions (e.g., leiden, louvain, seurat_clusters, celltypist).
    
    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return with scores. Default: False
        max_results (int): Maximum number of columns to return. Default: 5
        
    Returns:
        List[str] or List[Tuple[str, float]] or None:
            - If return_all=False: List of column names sorted by confidence
            - If return_all=True: List of (column_name, score) tuples
            - None if no columns found
            
    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> pred_cols = atlas.col_annotation_pred(adata)
        >>> print(f"Predicted columns: {pred_cols}")
        >>> 
        >>> # Get with scores
        >>> pred_with_scores = atlas.col_annotation_pred(adata, return_all=True)
        >>> for col, score in pred_with_scores:
        ...     print(f"{col}: {score:.3f}")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(
        min_reference_score=0.3,
        min_predicted_score=min_score
    )
    
    pred_candidates = results['annotation']['predicted'][:max_results]
    
    if not pred_candidates:
        return None
    
    if return_all:
        return pred_candidates
    else:
        return [col for col, score in pred_candidates]


def col_cluster(adata: AnnData,
                min_score: float = 0.5,
                return_all: bool = False,
                max_results: int = 5) -> Optional[List[str]]:
    """
    Detect cluster label columns in AnnData object.

    Uses the dedicated cluster-label detector (Leiden, Louvain, k‑means,
    Seurat clusters, PhenoGraph, etc.) which applies semantic 70 % +
    statistical 30 % scoring tuned for algorithmic clustering outputs.

    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return with scores. Default: False
        max_results (int): Maximum number of columns to return. Default: 5

    Returns:
        List[str] or List[Tuple[str, float]] or None:
            - If return_all=False: List of column names sorted by confidence
            - If return_all=True: List of (column_name, score) tuples
            - None if no columns found

    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> cluster_cols = atlas.col_cluster(adata)
        >>> print(f"Cluster columns: {cluster_cols}")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(min_cluster_score=min_score)

    clust_candidates = results["clustering"]["cluster_labels"][:max_results]

    if not clust_candidates:
        return None

    if return_all:
        return clust_candidates
    else:
        return [col for col, _score in clust_candidates]


def col_dimred(adata: AnnData,
               return_all: bool = False,
               max_results: int = 10) -> Optional[List[str]]:
    """
    Detect dimensionality reduction representations in AnnData.obsm.
    
    This function identifies embedding keys like X_pca, X_umap, X_tsne, etc.
    
    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        return_all (bool): If True, return with metadata. Default: False
        max_results (int): Maximum number of representations to return. Default: 10
        
    Returns:
        List[str] or List[Dict] or None:
            - If return_all=False: List of obsm keys (e.g., ['X_umap', 'X_pca'])
            - If return_all=True: List of dicts with 'key', 'shape', 'n_components'
            - None if no representations found
            
    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> dimred_keys = atlas.col_dimred(adata)
        >>> print(f"Dimensionality reductions: {dimred_keys}")
        >>> 
        >>> # Get with metadata
        >>> dimred_detailed = atlas.col_dimred(adata, return_all=True)
        >>> for emb in dimred_detailed:
        ...     print(f"{emb['key']}: {emb['n_components']} components")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters()
    
    embeddings = results['clustering']['embeddings'][:max_results]
    
    if not embeddings:
        return None
    
    if return_all:
        return [
            {
                'key': key,
                'shape': meta['shape'],
                'n_components': meta['n_components']
            }
            for key, meta in embeddings
        ]
    else:
        return [key for key, meta in embeddings]
