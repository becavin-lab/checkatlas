import logging
import os
import warnings

import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import checkatlas

EXTENSIONS_CELLRANGER = [".h5", ".mtx"]
CELLRANGER_FILE = "filtered_feature_bc_matrix.h5"
CELLRANGER_MATRIX_FILE = "matrix.mtx"
CELLRANGER_TYPE_OBSOLETE = "Cellranger < v3"
CELLRANGER_TYPE_CURRENT = "Cellranger >= v3"

logger = logging.getLogger("checkatlas")
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
sc.settings.verbosity = 0


def detect_cellranger(atlas_path: str) -> dict:
    if atlas_path.endswith(CELLRANGER_FILE):
        atlas_info = dict()
        atlas_info[checkatlas.ATLAS_NAME_KEY] = atlas_path.split(os.sep)[-3]
        atlas_info[checkatlas.ATLAS_TYPE_KEY] = CELLRANGER_TYPE_CURRENT
        atlas_info[checkatlas.ATLAS_EXTENSION_KEY] = EXTENSIONS_CELLRANGER[0]
        atlas_info[checkatlas.ATLAS_PATH_KEY] = atlas_path
        return atlas_info
    elif atlas_path.endswith(CELLRANGER_MATRIX_FILE):
        atlas_info = dict()
        atlas_info[checkatlas.ATLAS_NAME_KEY] = atlas_path.split(os.sep)[-5]
        atlas_info[checkatlas.ATLAS_TYPE_KEY] = CELLRANGER_TYPE_OBSOLETE
        atlas_info[checkatlas.ATLAS_EXTENSION_KEY] = EXTENSIONS_CELLRANGER[1]
        atlas_info[checkatlas.ATLAS_PATH_KEY] = atlas_path
        return atlas_info
    else:
        return dict()


def read_cellranger_current(atlas_info: dict) -> AnnData:
    """
    Read cellranger files.

    Load first /outs/filtered_feature_bc_matrix.h5
    Then add (if found):
    - Clustering
    - PCA-
    - UMAP
    - TSNE
    Args:
        atlas_path (dict): info on the atlas

    Returns:
        AnnData: scanpy object from cellranger
    """
    cellranger_path = atlas_info[checkatlas.ATLAS_PATH_KEY].replace(
        CELLRANGER_FILE, ""
    )
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
    adata = sc.read_10x_h5(atlas_info[checkatlas.ATLAS_PATH_KEY])
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


def read_cellranger_obsolete(atlas_info: dict) -> AnnData:
    """
    Read cellranger files.

    Load first /outs/filtered_feature_bc_matrix.h5
    Then add (if found):
    - Clustering
    - PCA-
    - UMAP
    - TSNE
    Args:
        atlas_path (dict): info on the atlas

    Returns:
        AnnData: scanpy object from cellranger
    """
    cellranger_path = atlas_info[checkatlas.ATLAS_PATH_KEY].replace(
        CELLRANGER_MATRIX_FILE, ""
    )
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
    adata = sc.read_10x_h5(atlas_info[checkatlas.ATLAS_PATH_KEY])
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
