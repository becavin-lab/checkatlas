import diopy
import os

import pandas as pd
import scanpy
import matplotlib.pyplot as plt

"""
checkatlas base module.

This is the principal module of the checkatlas project.
here you put your main classes and objects.

Be creative! do whatever you want!

If you want to replace this with a Flask application run:

    $ make init

and then choose `flask` as template.
"""

EXTENSIONS = ['.rds', '.h5ad']
RSCRIPT = 'checkatlas/convertSeurat.R'
SUMMARY_PATH = '_checkatlas.tsv'
ADATA_PATH = '_checkatlas_adata.tsv'
UMAP_PATH = '_checkatlas_umap.png'


class Atlas:
    def read_atlas(self) -> str:
        """
        Base method.
        """
        return "hello from BaseClass"

    def __call__(self) -> str:
        return ""


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


def convert_seurat_atlases(path, atlas_list) -> str:
    for atlas_path in atlas_list:
        if atlas_path.endswith('.rds'):
            atlas_h5 = atlas_path.replace('.rds', '.h5ad')
            if os.path.exists(atlas_h5):
                print("Seurat file already converted to Scanpy:", atlas_h5)
            else:
                print('Convert Seurat object to Scanpy: ', atlas_path)
                atlas_name = os.path.basename(atlas_path).replace('.rds', '')
                rscript_cmd = 'Rscript ' + RSCRIPT + ' ' + os.path.dirname(atlas_path) + ' ' + atlas_name
                print(rscript_cmd)
                os.system(rscript_cmd)


def parse_atlases(path, atlas_list) -> str:
    list_atlas_name = list()
    for atlas_path in atlas_list:
        if atlas_path.endswith('.h5ad'):
            list_atlas_name.append(os.path.basename(atlas_path).strip().replace('.h5ad', ''))
    for atlas_path in atlas_list:
        if atlas_path.endswith('.h5ad'):
            atlas_name = os.path.basename(atlas_path).strip().replace('.h5ad', '')
            print(atlas_path)
            # Create summary table
            header = ['Layers', 'NbCells', 'NbGenes', 'RawData', 'Normalized']
            df_summary = pd.DataFrame(index=[atlas_name], columns=header)
            adata = scanpy.read_h5ad(atlas_path)
            print(adata)
            df_summary['Layers'][atlas_name] = len(adata.layers)
            df_summary['NbCells'][atlas_name] = adata.n_obs
            df_summary['NbGenes'][atlas_name] = adata.n_vars
            df_summary['RawData'][atlas_name] = adata.raw is None
            df_summary['Normalized'][atlas_name] = adata.raw is None
            df_summary.to_csv(atlas_path.replace('.h5ad', SUMMARY_PATH), index=False)

            # Create AnnData table
            header = ['obs', 'obsm', 'var', 'varm', 'uns']
            df_summary = pd.DataFrame(index=[atlas_name], columns=header)
            adata = scanpy.read_h5ad(atlas_path)
            print(adata)
            df_summary['obs'][atlas_name] = ';'.join(list(adata.obs.columns))
            df_summary['obsm'][atlas_name] = ';'.join(list(adata.obsm_keys()))
            df_summary['var'][atlas_name] = ';'.join(list(adata.var_keys()))
            df_summary['varm'][atlas_name] = ';'.join(list(adata.varm_keys()))
            df_summary['uns'][atlas_name] = ';'.join(list(adata.uns_keys()))
            df_summary.to_csv(atlas_path.replace('.h5ad', ADATA_PATH), index=False)

            # Create umap
            # search cluster obs
            obs_clusters = ['CellType', 'celltype', 'seurat_clusters', 'orig.ident']
            found = False
            i = 0
            # Setting up figures directory
            scanpy.settings.figdir = path
            print(scanpy.settings.figdir)
            if not os.path.exists(path + '/umap'):
                os.mkdir(path + '/umap')
            # Exporting umap
            while not found:
                if obs_clusters[i] in adata.obs_keys():
                    fig = scanpy.pl.umap(adata, color=obs_clusters[i], show=False, save='/' + atlas_name + UMAP_PATH)
                    found = True
                i += 1
            if not found:
                fig = scanpy.pl.umap(adata, show=False, save='/' + atlas_name + UMAP_PATH)