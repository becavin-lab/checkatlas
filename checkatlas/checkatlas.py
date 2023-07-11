import logging
import os

import pandas as pd

from . import atlas, cellranger, seurat
from .utils import files as chk_files

"""
checkatlas base module.
This is the principal module of the checkatlas project.

"""

PROCESS_TYPE = [
    "summary",
    "qc",
    "metric_cluster",
    "metric_annot",
    "metric_dimred",
]


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
ATLAS_TABLE_HEADER = [
    ATLAS_NAME_KEY,
    ATLAS_TYPE_KEY,
    ATLAS_EXTENSION_KEY,
    ATLAS_PATH_KEY,
]


logger = logging.getLogger("checkatlas")


def list_all_atlases(checkatlas_path: str) -> None:
    # Get all files with matching extension
    EXTENSIONS = [
        atlas.ANNDATA_EXTENSION,
        cellranger.CELLRANGER_FILE,
        cellranger.CELLRANGER_MATRIX_FILE,
        seurat.SEURAT_EXTENSION,
    ]
    atlas_list = list()
    for root, dirs, files in os.walk(checkatlas_path):
        for file in files:
            for extension in EXTENSIONS:
                if file.endswith(extension):
                    atlas_list.append(os.path.join(root, file))

    # Filter the lists keepng only atlases
    clean_scanpy_list = list()
    clean_cellranger_list = list()
    clean_seurat_list = list()
    for atlas_path in atlas_list:
        atlas_info = atlas.detect_scanpy(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_scanpy_list.append(atlas_info)
        atlas_info = cellranger.detect_cellranger(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_cellranger_list.append(atlas_info)
        atlas_info = seurat.detect_seurat(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a seurat output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_seurat_list.append(atlas_info)

    # Save the list of atlas taken into account
    chk_files.save_list_scanpy(clean_scanpy_list, checkatlas_path)
    chk_files.save_list_cellranger(clean_cellranger_list, checkatlas_path)
    chk_files.save_list_seurat(clean_seurat_list, checkatlas_path)


def read_list_atlases(checkatlas_path: str) -> tuple:
    clean_scanpy_list = pd.read_csv(
        chk_files.get_table_scanpy_path(checkatlas_path)
    )
    clean_scanpy_list.index = clean_scanpy_list[ATLAS_NAME_KEY]
    clean_cellranger_list = pd.read_csv(
        chk_files.get_table_cellranger_path(checkatlas_path)
    )
    clean_cellranger_list.index = clean_cellranger_list[ATLAS_NAME_KEY]
    clean_seurat_list = pd.read_csv(
        chk_files.get_table_seurat_path(checkatlas_path)
    )
    clean_seurat_list.index = clean_seurat_list[ATLAS_NAME_KEY]
    return clean_scanpy_list, clean_cellranger_list, clean_seurat_list


if __name__ == "__main__":
    path = "/data/analysis/data_becavin/checkatlas_test/tuto"
    path = "/Users/christophebecavin/Documents/testatlas/"
    # atlas_path = "/Users/christophebecavin/Documents/testatlas/"
    atlas_info = ["test_version", "Scanpy", ".h5ad", "data/test_version.h5ad"]
    # folders.checkatlas_folders(path)
    # atlas_list = list_atlases(path)
