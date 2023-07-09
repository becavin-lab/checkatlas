import argparse
import logging
import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata import _io as _io

from checkatlas import checkatlas

from .metrics import metrics
from .utils import files, folders

"""
Atlas module
All the function to screen the atlases
"""

OBS_CLUSTERS = [
    "cell_type",
    "CellType",
    "celltype",
    "ann_finest_level",
    "cellranger_graphclust",
    "seurat_clusters",
    "RNA_snn_res.0.5",
    "louvain",
    "leiden",
    "orig.ident",
]

OBSM_DIMRED = [
    "X_umap",
    "X_pca",
    "X-tsne",
]

OBS_QC = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "percent.mito",
    "percent.ribo",
    "dropouts",
    "nCount_SCT",
    "nFeature_SCT",
    "nCount_RNA",
    "nFeature_RNA",
    "nCount_HTO",
    "nFeature_HTO",
    "n_genes",
    "n_genes_by_counts",
    "total_counts",
    "total_counts_mt",
    "pct_counts_mt",
    "total_counts_ribo",
    "pct_counts_ribo",
]

CELLINDEX_HEADER = "cell_index"

logger = logging.getLogger("checkatlas")


def read_atlas(atlas_path: str) -> AnnData:
    """
    Read Scanpy or Cellranger data : .h5ad or .h5

    Args:
        atlas_path (str): path of the .h5ad atlas

    Returns:
        AnnData: scanpy object from .h5ad
    """
    logger.info(
        f"Load {checkatlas.get_atlas_name(atlas_path)} "
        f"in {checkatlas.get_atlas_directory(atlas_path)}"
    )
    try:
        if atlas_path.endswith(".h5"):
            logger.debug(f"Read Cellranger results {atlas_path}")
            adata = read_cellranger(atlas_path)
        else:
            logger.debug(f"Read Scanpy file {atlas_path}")
            adata = sc.read_h5ad(atlas_path)
        return adata
    except _io.utils.AnnDataReadError:
        logger.warning(
            f"AnnDataReadError, cannot read: "
            f"{checkatlas.get_atlas_name(atlas_path)}"
        )
        return None


def read_cellranger(atlas_path: str) -> AnnData:
    """
    Read cellranger files.

    Load first /outs/filtered_feature_bc_matrix.h5
    Then add (if found):
    - Clustering
    - PCA-
    - UMAP
    - TSNE
    Args:
        atlas_path (str): path of the atlas

    Returns:
        AnnData: scanpy object from cellranger
    """
    cellranger_path = atlas_path.replace(checkatlas.CELLRANGER_FILE, "")
    cellranger_path = os.path.join(cellranger_path, "outs")
    clust_path = os.path.join(
        cellranger_path, "analysis", "clustering", "graphclust", "clusters.csv"
    )
    rna_umap = os.path.join(
        cellranger_path, "analysis", "umap", "2_components", "projection.csv"
    )
    rna_tsne = os.path.join(
        cellranger_path, "analysis", "tsne", "2_components", "projection.csv"
    )
    rna_pca = os.path.join(
        cellranger_path, "analysis", "pca", "10_components", "projection.csv"
    )
    adata = sc.read_10x_h5(atlas_path)
    adata.var_names_make_unique()
    # Add cluster
    if os.path.exists(clust_path):
        df_cluster = pd.read_csv(clust_path, index_col=0)
        adata.obs["cellranger_graphclust"] = df_cluster["Cluster"]
    # Add reduction
    if os.path.exists(rna_umap):
        df_umap = pd.read_csv(rna_umap, index_col=0)
        adata.obsm["X_umap"] = df_umap
    if os.path.exists(rna_tsne):
        df_tsne = pd.read_csv(rna_tsne, index_col=0)
        adata.obsm["X_tsne"] = df_tsne
    if os.path.exists(rna_pca):
        df_pca = pd.read_csv(rna_pca, index_col=0)
        adata.obsm["X_pca"] = df_pca
    return adata


def clean_scanpy_atlas(adata: AnnData, atlas_path: str) -> AnnData:
    """
    Clean the Scanpy object to be sure to get all information out of it

    - Make var names unique
    - Make var unique for Raw matrix
    - If OBS_CLUSTERS are present and in int32 -> be sure to
    transform them in categorical

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path to the atlas

    Returns:
        AnnData: cleaned atlas
    """
    logger.debug(f"Clean scanpy: {checkatlas.get_atlas_name(atlas_path)}")
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
                if type(adata.obs[obs_key].dtype) == pd.CategoricalDtype:
                    obs_keys.append(obs_key)
    # Remove keys with only one category
    obs_keys_final = list()
    for obs_key in obs_keys:
        annotations = adata.obs[obs_key]
        categories_temp = annotations.cat.categories
        # remove nan if found
        categories = categories_temp.dropna()
        if True in categories.isin(["nan"]):
            index = categories.get_loc("nan")
            categories = categories.delete(index)
        # Add obs_key with more than one category (with Nan removed)
        if len(categories) != 1:
            logger.debug(f"Add obs_key {obs_key} with cat {categories_temp}")
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
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Create a table with all summarizing variables

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    logger.debug(f"Create Summary table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.SUMMARY, checkatlas.SUMMARY_EXTENSION, args.path
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
    df_summary["AtlasFileType"][atlas_name] = checkatlas.get_atlas_type(
        atlas_path
    )
    df_summary["NbCells"][atlas_name] = adata.n_obs
    df_summary["NbGenes"][atlas_name] = adata.n_vars
    df_summary["AnnData.raw"][atlas_name] = adata.raw is not None
    df_summary["AnnData.X"][atlas_name] = adata.X is not None
    df_summary["File_extension"][atlas_name] = checkatlas.get_atlas_extension(
        atlas_path
    )
    df_summary["File_path"][atlas_name] = atlas_path.replace(args.path, "")
    df_summary.to_csv(csv_path, index=False, sep="\t")


def create_anndata_table(
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Create an html table with all AnnData arguments
    The html code will make all elements of the table visible in MultiQC
    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    logger.debug(f"Create Adata table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.ANNDATA, checkatlas.ADATA_EXTENSION, args.path
    )
    # Create AnnData table
    header = ["obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    # html_element = "<span class=\"label label-primary\">"
    # new_line = ''
    # for value in list(adata.obs.columns):
    #     new_line += html_element + value + "</span><br>"
    #     print(new_line)
    df_summary["obs"][atlas_name] = (
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
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Display the atlas QC table
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    qc_path = files.get_file_path(
        atlas_name, folders.QC, checkatlas.QC_EXTENSION, args.path
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
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Display the atlas QC plot
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    sc.settings.figdir = folders.get_workingdir(args.path)
    sc.set_figure_params(dpi_save=80)
    qc_path = os.sep + atlas_name + checkatlas.QC_FIG_EXTENSION
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
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
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
        umap_path = os.sep + atlas_name + checkatlas.UMAP_EXTENSION
        # Exporting umap
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.umap(adata, color=obs_keys[0], show=False, save=umap_path)
        else:
            sc.pl.umap(adata, show=False, save=umap_path)


def create_tsne_fig(
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Display the TSNE of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
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
        tsne_path = os.sep + atlas_name + checkatlas.TSNE_EXTENSION
        # Exporting tsne
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.tsne(adata, color=obs_keys[0], show=False, save=tsne_path)
        else:
            sc.pl.tsne(adata, show=False, save=tsne_path)


def metric_cluster(
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Calc clustering metrics

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name,
        folders.CLUSTER,
        checkatlas.METRIC_CLUSTER_EXTENSION,
        args.path,
    )
    header = ["Sample", "obs"] + args.metric_cluster
    df_cluster = pd.DataFrame(columns=header)
    obs_keys = get_viable_obs_annot(adata, args)
    obsm_key_representation = "X_umap"

    if len(obs_keys) > 0:
        logger.debug(f"Calc clustering metrics for {atlas_name}")
        for obs_key in obs_keys:
            dict_line = {
                "Sample": [atlas_name + "_" + obs_key],
                "obs": [obs_key],
            }
            for metric in args.metric_cluster:
                logger.debug(
                    f"Calc {metric} for {atlas_name} "
                    f"with obs {obs_key} and obsm {obsm_key_representation}"
                )
                metric_value = metrics.calc_metric_cluster_scanpy(
                    metric, adata, obs_key, obsm_key_representation
                )
                dict_line[metric] = metric_value
            df_line = pd.DataFrame(dict_line)
            df_cluster = pd.concat(
                [df_cluster, df_line], ignore_index=True, axis=0
            )
        df_cluster.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obs_key was found for {atlas_name}")


def metric_annot(
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Calc annotation metrics

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name,
        folders.ANNOTATION,
        checkatlas.METRIC_ANNOTATION_EXTENSION,
        args.path,
    )
    header = ["Sample", "Reference", "obs"] + args.metric_annot
    df_annot = pd.DataFrame(columns=header)
    obs_keys = get_viable_obs_annot(adata, args)
    if len(obs_keys) > 1:
        logger.debug(f"Calc annotation metrics for {atlas_name}")
        ref_obs = obs_keys[0]
        for i in range(1, len(obs_keys)):
            obs_key = obs_keys[i]
            dict_line = {
                "Sample": [atlas_name + "_" + obs_key],
                "Reference": [ref_obs],
                "obs": [obs_key],
            }
            for metric in args.metric_annot:
                logger.debug(
                    f"Calc {metric} for {atlas_name} "
                    f"with obs {obs_key} vs ref_obs {ref_obs}"
                )
                metric_value = metrics.calc_metric_annot_scanpy(
                    metric, adata, obs_key, ref_obs
                )
                dict_line[metric] = metric_value
            df_line = pd.DataFrame(dict_line)
            df_annot = pd.concat(
                [df_annot, df_line], ignore_index=True, axis=0
            )
        df_annot.to_csv(csv_path, index=False, sep="\t")
    else:
        logger.debug(f"No viable obs_key was found for {atlas_name}")


def metric_dimred(
    adata: AnnData, atlas_path: str, args: argparse.Namespace
) -> None:
    """
    Calc dimensionality reduction metrics

    Args:
        adata (AnnData): atlas to analyse
        atlas_path (str): path of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name,
        folders.DIMRED,
        checkatlas.METRIC_DIMRED_EXTENSION,
        args.path,
    )
    header = ["Sample", "obsm"] + args.metric_dimred
    df_dimred = pd.DataFrame(columns=header)
    obsm_keys = get_viable_obsm(adata, args)
    if len(obsm_keys) > 0:
        logger.debug(f"Calc dim red metrics for {atlas_name}")
        for obsm_key in obsm_keys:
            dict_line = {
                "Sample": [atlas_name + "_" + obsm_key],
                "obsm": [obsm_key],
            }
            for metric in args.metric_dimred:
                logger.debug(
                    f"Calc {metric} for {atlas_name} with obsm {obsm_key}"
                )
                metric_value = metrics.calc_metric_dimred(
                    metric, adata, obsm_key
                )
                dict_line[metric] = metric_value
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
