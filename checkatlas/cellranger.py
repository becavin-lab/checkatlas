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
    cellranger_out_path = os.path.dirname(
        atlas_info[checkatlas.ATLAS_PATH_KEY]
    )
    cellranger_analysis_path = os.path.join(cellranger_out_path, "analysis")
    cellranger_clust_path = os.path.join(
        cellranger_analysis_path, "clustering"
    )
    cellranger_umap_path = os.path.join(cellranger_analysis_path, "umap")
    cellranger_tsne_path = os.path.join(cellranger_analysis_path, "tsne")
    cellranger_pca_path = os.path.join(cellranger_analysis_path, "pca")

    # Search graphclust
    graphclust_path = ""
    for root, dirs, files in os.walk(cellranger_clust_path):
        for dir in dirs:
            if dir.endswith("graphclust"):
                cluster_path = os.path.join(root, dir, "clusters.csv")
                if os.path.exists(cluster_path) and not root.endswith("atac"):
                    graphclust_path = cluster_path
                    break
    # Search kmeans
    kmeans_path = ""
    k_value = 0
    found_kmeans = False
    for root, dirs, files in os.walk(cellranger_clust_path):
        for dir in dirs:
            # Search the highest kmeans = 10
            dir_prefix = "kmeans_10"
            if dir_prefix in dir and not found_kmeans:
                cluster_path = os.path.join(root, dir, "clusters.csv")
                if os.path.exists(cluster_path):
                    kmeans_path = cluster_path
                    k_value = 10
                    found_kmeans = True
                    break
            # Or search the highest kmeans = 5 (for multiome atlas)
            dir_prefix = os.path.join("gex", "kmeans_5")
            if dir_prefix in os.path.join(root, dir) and not found_kmeans:
                cluster_path = os.path.join(root, dir, "clusters.csv")
                if os.path.exists(cluster_path):
                    kmeans_path = cluster_path
                    k_value = 5
                    found_kmeans = True
                    break

    # Search umap
    rna_umap = ""
    for root, dirs, files in os.walk(cellranger_umap_path):
        for file in files:
            if file.endswith("projection.csv") and not root.endswith("atac"):
                rna_umap = os.path.join(root, file)
                break

    # Search t-SNE
    rna_tsne = ""
    for root, dirs, files in os.walk(cellranger_tsne_path):
        for file in files:
            if file.endswith("projection.csv") and not root.endswith("atac"):
                rna_tsne = os.path.join(root, file)
                break

    rna_pca = ""
    for root, dirs, files in os.walk(cellranger_pca_path):
        for file in files:
            if file.endswith("projection.csv") and not root.endswith("atac"):
                rna_pca = os.path.join(root, file)
                break

    # Manage multiome cellranger files
    dim_red_path = os.path.join(
        cellranger_analysis_path, "dimensionality_reduction"
    )
    if os.path.exists(dim_red_path):
        gex_path = os.path.join(dim_red_path, "gex")
        if os.path.exists(gex_path):
            rna_umap = os.path.join(gex_path, "umap_projection.csv")
            rna_tsne = os.path.join(gex_path, "tsne_projection.csv")
            rna_pca = os.path.join(gex_path, "pca_projection.csv")

    # Read 10x h5 file
    adata = sc.read_10x_h5(atlas_info[checkatlas.ATLAS_PATH_KEY])
    adata.var_names_make_unique()

    # Add cluster
    if os.path.exists(graphclust_path):
        df_cluster = pd.read_csv(graphclust_path, index_col=0)
        adata.obs["cellranger_graphclust"] = df_cluster["Cluster"]
    if os.path.exists(kmeans_path):
        df_cluster = pd.read_csv(kmeans_path, index_col=0)
        adata.obs["cellranger_kmeans_" + str(k_value)] = df_cluster["Cluster"]

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
    cellranger_out_path = os.path.join(cellranger_path, os.pardir, os.pardir)
    cellranger_analysis_path = os.path.join(
        cellranger_out_path, "analysis_csv"
    )

    cellranger_umap_path = os.path.join(cellranger_analysis_path, "umap")
    cellranger_tsne_path = os.path.join(cellranger_analysis_path, "tsne")
    cellranger_pca_path = os.path.join(cellranger_analysis_path, "pca")
    print(cellranger_out_path)
    print(cellranger_analysis_path)
    print(cellranger_umap_path)

    # Search graphclust
    graphclust_path = ""
    for root, dirs, files in os.walk(cellranger_out_path):
        for dir in dirs:
            if dir.endswith("graphclust"):
                cluster_path = os.path.join(root, dir, "clusters.csv")
                if os.path.exists(cluster_path):
                    graphclust_path = cluster_path
                    break
    # Search kmeans
    kmeans_path = ""
    k_value = 0
    for root, dirs, files in os.walk(cellranger_out_path):
        for dir in dirs:
            if dir.endswith("kmeans"):
                # Search the highest kmeans from 15 to 3
                for k in reversed(range(3, 16)):
                    cluster_path = os.path.join(
                        root, dir, str(k) + "_clusters", "clusters.csv"
                    )
                    if os.path.exists(cluster_path):
                        kmeans_path = cluster_path
                        k_value = k
                        break

    rna_umap = os.path.join(cellranger_umap_path, "projection.csv")
    rna_tsne = os.path.join(cellranger_tsne_path, "projection.csv")
    rna_pca = os.path.join(cellranger_pca_path, "projection.csv")

    # get matrix folder
    matrix_folder = os.path.dirname(atlas_info[checkatlas.ATLAS_PATH_KEY])
    adata = sc.read_10x_mtx(matrix_folder)
    adata.var_names_make_unique()

    # Add cluster
    if os.path.exists(graphclust_path):
        df_cluster = pd.read_csv(graphclust_path, index_col=0)
        adata.obs["cellranger_graphclust"] = df_cluster["Cluster"]
    if os.path.exists(kmeans_path):
        df_cluster = pd.read_csv(kmeans_path, index_col=0)
        adata.obs["cellranger_kmeans_" + str(k_value)] = df_cluster["Cluster"]

    # Add reduction
    if os.path.exists(rna_umap):
        df_umap = pd.read_csv(rna_umap, index_col=0)
        adata.obsm["X_umap"] = df_umap
    if os.path.exists(rna_tsne):
        df_tsne = pd.read_csv(rna_tsne, index_col=0)
        if len(df_tsne) == len(adata):
            adata.obsm["X_tsne"] = df_tsne
    if os.path.exists(rna_pca):
        df_pca = pd.read_csv(rna_pca, index_col=0)
        if len(df_pca) == len(adata):
            adata.obsm["X_pca"] = df_pca
    return adata
