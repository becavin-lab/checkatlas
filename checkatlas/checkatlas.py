import os
import pandas as pd
import scanpy
from . import atlas
from dask.distributed import Client
from dask.distributed import wait
import webbrowser

"""
checkatlas base module.
This is the principal module of the checkatlas project.

"""

EXTENSIONS = [".rds", ".h5ad"]
RSCRIPT = "checkatlas/convertSeurat.R"
SUMMARY_PATH = "_checkatlas.tsv"
ADATA_PATH = "_checkatlas_adata.tsv"
UMAP_PATH = "_checkatlas_umap.png"


def list_atlases(path) -> list:
    """
    List all atlases files in the path
    Detect .rds, .h5, .h5ad
    :param path:
    :return: List of files
    """
    atlas_list = list()
    for root, dirs, files in os.walk(path):
        for file in files:
            for extension in EXTENSIONS:
                if file.endswith(extension):
                    atlas_list.append(os.path.join(root, file))
    return atlas_list


def convert_seurat_atlases(path, atlas_list) -> None:
    """
    Convert all Seurat atlas to Scanpy
    :param path:
    :param atlas_list:
    :return:
    """
    for atlas_path in atlas_list:
        if atlas_path.endswith(".rds"):
            atlas_h5 = atlas_path.replace(".rds", ".h5ad")
            if os.path.exists(atlas_h5):
                print("Seurat file already converted to Scanpy:", atlas_h5)
            else:
                atlas.convert_atlas(atlas_path)




def parse_atlases(atlas_path, path, atlas_name) -> None:
    """
    Main function of checkatlas
    For every atlas create summary tables with all attributes of the atlas
    Calc UMAP, tSNE, andd all metrics
    :param atlas_path:
    :return:
    """
    adata = scanpy.read_h5ad(atlas_path)

    # create summary table
    atlas.create_summary_table(adata, atlas_name, atlas_path)

    # Create AnnData table
    atlas.create_anndata_table(adata, atlas_name, atlas_path)

    # Create umap
    #atlas.create_umap_fig(path, adata, atlas_name, atlas_path)


def calc_checkatlas(path, atlas_list) -> None:
    """
    Run parsing of atlas in distributed threads
    :param path:
    :param atlas_list:
    :return:
    """
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = list()
        for atlas_path in atlas_list:
            if atlas_path.endswith(".h5ad"):
                atlas_name = (
                    os.path.basename(atlas_path).strip().replace(".h5ad", "")
                )
                print('Running thread for: ',atlas_path)
                future_temp = client.submit(parse_atlases, atlas_path, path, atlas_name, key=atlas_name)
                futures.append(future_temp)
        print('Running parsing of adata')
        print('Opening Dask browser: http://localhost:8787/status')
        webbrowser.open('http://localhost:8787/status')
        wait(futures)