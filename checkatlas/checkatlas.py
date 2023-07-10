import logging
import os

import pandas as pd

from . import atlas, cellranger, seurat

"""
checkatlas base module.
This is the principal module of the checkatlas project.

"""

EXTENSIONS = [".rds", ".h5ad", ".h5", ".mtx"]
SCANPY_EXTENSION = ".h5ad"
CELLRANGER_EXTENSION = ".h5"
SEURAT_EXTENSION = ".rds"

CELLRANGER_FILE = "filtered_feature_bc_matrix.h5"
CELLRANGER_MATRIX_FILE = "matrix.mtx"

SUMMARY_EXTENSION = "_checkatlas_summ.tsv"
ADATA_EXTENSION = "_checkatlas_adata.tsv"
QC_FIG_EXTENSION = "_checkatlas_qc.png"
QC_EXTENSION = "_checkatlas_qc.tsv"
UMAP_EXTENSION = "_checkatlas_umap.png"
TSNE_EXTENSION = "_checkatlas_tsne.png"
METRIC_CLUSTER_EXTENSION = "_checkatlas_mcluster.tsv"
METRIC_ANNOTATION_EXTENSION = "_checkatlas_mannot.tsv"
METRIC_DIMRED_EXTENSION = "_checkatlas_mdimred.tsv"
METRIC_SPECIFICITY_EXTENSION = "_checkatlas_mspecificity.tsv"

ATLAS_NAME_KEY = "Atlas_name"
ATLAS_TYPE_KEY = "Atlas_type"
ATLAS_EXTENSION_KEY = "Atlas_extension"
ATLAS_PATH_KEY = "Atlas_path"


logger = logging.getLogger("checkatlas")


def list_atlases(path: str) -> list:
    """
    List all atlases files in the path
    Detect .rds, .h5, .h5ad

    Args:
        path: Path for searching single-cell atlases.

    Returns:
        list: List of file atlases to check.
    """
    atlas_list = list()
    for root, dirs, files in os.walk(path):
        for file in files:
            for extension in EXTENSIONS:
                if file.endswith(extension):
                    atlas_list.append(os.path.join(root, file))
    print(atlas_list)
    return atlas_list


def get_atlas_name(atlas_path: str) -> str:
    """
    From atlas_path extract the atlas_name
    Args:
        atlas_path:
    Returns:
        str: The atlas_name
    """
    atlas_type = get_atlas_type(atlas_path)
    if atlas_type == "CellRanger":
        # If cellranger take the name of the first folder
        return atlas_path.split(os.sep)[-3]
    else:
        return os.path.splitext(os.path.basename(atlas_path))[0]


def get_atlas_type(atlas_path: str) -> str:
    """
    Return the type of atlas using its extension
    TO DO : Need to be more smart then just using extension.
    It should use function type() or class() in R
    Args:
        atlas_path (str): path of the atlas

    Returns:
        str: Type of atlas among : Scanpy, Cellranger, Seurat
    """
    atlas_extension = get_atlas_extension(atlas_path)
    if atlas_extension == SCANPY_EXTENSION:
        return "Scanpy"
    elif atlas_extension == CELLRANGER_EXTENSION:
        return "CellRanger"
    else:
        return "Seurat"


def get_atlas_extension(atlas_path: str) -> str:
    """
    From atlas_path extract the atlas file extension
    Args:
        atlas_path:
    Returns:
        None
    """
    return os.path.splitext(os.path.basename(atlas_path))[1]


def get_atlas_directory(atlas_path: str) -> str:
    """
    From atlas_path extract the atlas directory
    Args:
        atlas_path:
    Returns:
        None
    """
    return os.path.dirname(atlas_path)


def clean_list_atlases(atlas_list: list, checkatlas_path: str) -> tuple:
    """
    Go through all files and detect Seurat, CellRanger or Scanpy Atlas
    The "cleaning means that we test if the atlas is valid or not.
    Args:
        atlas_list: list of atlases found with proper extension
        checkatlas_path: the path where checkatlas files are saved
    Returns:
         tuple: clean_atlas_scanpy, clean_atlas_seurat, clean_atlas_cellranger
    """

    clean_atlas_scanpy = list()
    clean_atlas_seurat = list()
    clean_atlas_cellranger = list()
    for atlas_path in atlas_list:
        atlas_info = atlas.detect_scanpy(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_atlas_scanpy.append(atlas_info)
        atlas_info = cellranger.detect_cellranger(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_atlas_cellranger.append(atlas_info)
        atlas_info = seurat.detect_seurat(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_atlas_seurat.append(atlas_info)
    # Save the list of atlas taken into account
    df_atlas = pd.DataFrame.from_dict(atlas_info)
    print(df_atlas)
    return clean_atlas_scanpy, clean_atlas_seurat, clean_atlas_cellranger


def get_pipeline_functions(module, args) -> list:
    """
    Using arguments of checkatlas program -> build
    the list of functions to run on each adata
    and seurat object
    Args:
        module: Module to use either atlas or atlas_seurat
        args: List of args for checkatlas program
    Returns:
         list: list of functions to run
    """
    checkatlas_functions = list()

    if not args.NOADATA:
        checkatlas_functions.append(module.create_anndata_table)
    if not args.NOQC:
        if "violin_plot" in args.qc_display:
            checkatlas_functions.append(module.create_qc_plots)
        if (
            "total-counts" in args.qc_display
            or "n_genes_by_counts" in args.qc_display
            or "pct_counts_mt" in args.qc_display
        ):
            checkatlas_functions.append(module.create_qc_tables)
    if not args.NOREDUCTION:
        checkatlas_functions.append(module.create_umap_fig)
        checkatlas_functions.append(module.create_tsne_fig)
    if not args.NOMETRIC:
        if len(args.metric_cluster) > 0:
            checkatlas_functions.append(module.metric_cluster)
        else:
            logger.debug(
                "No clustering metric was specified in --metric_cluster"
            )
        if len(args.metric_annot) > 0:
            checkatlas_functions.append(module.metric_annot)
        else:
            logger.debug(
                "No annotation metric was specified in --metric_annot"
            )
        if len(args.metric_dimred) > 0:
            checkatlas_functions.append(module.metric_dimred)
        else:
            logger.debug("No dim red metric was specified in --metric_dimred")
    # Create summary by default, it is ran at last so it marks
    # the end of the pipeline
    # This table is then used by the resume option
    checkatlas_functions.append(module.create_summary_table)
    return checkatlas_functions


if __name__ == "__main__":
    path = "/data/analysis/data_becavin/checkatlas_test/tuto"
    path = "/Users/christophebecavin/Documents/testatlas/"
    # atlas_path = "/Users/christophebecavin/Documents/testatlas/"
    atlas_info = ["test_version", "Scanpy", ".h5ad", "data/test_version.h5ad"]
    # folders.checkatlas_folders(path)
    # atlas_list = list_atlases(path)
