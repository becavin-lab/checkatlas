import os
import pandas as pd
import scanpy
from . import checkatlas

"""
Atlas moduel
All the function to screen the atlases
"""

OBS_CLUSTERS = [
    "CellType",
    "celltype",
    "seurat_clusters",
    "orig.ident",
]


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
            for extension in checkatlas.EXTENSIONS:
                if file.endswith(extension):
                    atlas_list.append(os.path.join(root, file))
    return atlas_list


def convert_atlas(atlas_path) -> None:
    """
    Convert a atlas to Scanpy
    :param atlas_path:
    :return:
    """
    print("Convert Seurat object to Scanpy: ", atlas_path)
    atlas_name = os.path.basename(atlas_path).replace(".rds", "")
    rscript_cmd = (
            "Rscript "
            + checkatlas.RSCRIPT
            + " "
            + os.path.dirname(atlas_path)
            + " "
            + atlas_name
    )
    print(rscript_cmd)
    os.system(rscript_cmd)


def create_summary_table(adata, atlas_name, atlas_path) -> None:
    """
    Create a table with all interesting variables
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    # Create summary table
    header = ["Layers", "NbCells", "NbGenes", "RawData", "Normalized"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["Layers"][atlas_name] = len(adata.layers)
    df_summary["NbCells"][atlas_name] = adata.n_obs
    df_summary["NbGenes"][atlas_name] = adata.n_vars
    df_summary["RawData"][atlas_name] = adata.raw is None
    df_summary["Normalized"][atlas_name] = adata.raw is None
    df_summary.to_csv(
        atlas_path.replace(".h5ad", checkatlas.SUMMARY_PATH), index=False
    )


def create_anndata_table(adata, atlas_name, atlas_path) -> None:
    """
    Create a table with all AnnData arguments
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    # Create AnnData table
    header = ["obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["obs"][atlas_name] = ";".join(list(adata.obs.columns))
    df_summary["obsm"][atlas_name] = ";".join(list(adata.obsm_keys()))
    df_summary["var"][atlas_name] = ";".join(list(adata.var_keys()))
    df_summary["varm"][atlas_name] = ";".join(list(adata.varm_keys()))
    df_summary["uns"][atlas_name] = ";".join(list(adata.uns_keys()))
    df_summary.to_csv(
        atlas_path.replace(".h5ad", checkatlas.ADATA_PATH), index=False
    )


def create_umap_fig(path, adata, atlas_name, atlas_path) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation
    :param path:
    :param adata:
    :param atlas_name:
    :param atlas_path:
    :return:
    """
    # Create umap
    found = False
    i = 0
    # Setting up figures directory
    scanpy.settings.figdir = path
    print(scanpy.settings.figdir)
    if not os.path.exists(path + "/umap"):
        os.mkdir(path + "/umap")
    # Exporting umap
    while not found:
        if OBS_CLUSTERS[i] in adata.obs_keys():
            scanpy.pl.umap(
                adata,
                color=OBS_CLUSTERS[i],
                show=False,
                save="/" + atlas_name + checkatlas.UMAP_PATH,
            )
            found = True
        i += 1
    if not found:
        scanpy.pl.umap(
            adata, show=False, save="/" + atlas_name + checkatlas.UMAP_PATH
        )