import csv
import inspect
import os
import webbrowser

import anndata
import matplotlib
import scanpy as sc
from dask.distributed import Client, LocalCluster, wait

from . import atlas, folders, multiqc

# try:
#     from . import atlas, folders, multiqc
# except ImportError:
#     import atlas
#     import folders
#     import multiqc


"""
checkatlas base module.
This is the principal module of the checkatlas project.

"""

EXTENSIONS = [".h5ad", ".h5"]
# EXTENSIONS = [".rds", ".h5ad", ".h5"]
CELLRANGER_FILE = "/outs/filtered_feature_bc_matrix.h5"
RSCRIPT = inspect.getfile(atlas).replace("atlas.py", "convertSeurat.R")
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


def get_atlas_name(atlas_path):
    """
    From atlas_path extract the atlas_name
    :param atlas_path:
    :return:
    """
    return os.path.splitext(os.path.basename(atlas_path))[0]


def get_atlas_extension(atlas_path):
    """
    From atlas_path extract the atlas file extension
    :param atlas_path:
    :return:
    """
    return os.path.splitext(os.path.basename(atlas_path))[1]


def convert_seurat_atlas(atlas_path, atlas_name) -> bool:
    """
    Convert all Seurat atlas to Scanpy
    :param path:
    :param atlas_list:
    :return:
    """
    atlas.convert_atlas(atlas_path, atlas_name)
    return True


def clean_list_atlases(atlas_list, path) -> dict:
    """
    Go through all files and detect Seurat or Scanpy Atlas
    Then:
    - Convert Seurat files to Scanpy
    - Clean Scanpy files
    :param atlas_list:
    :return: clean_list_atlas will only cleaned Scanpy atlas. A dict with
    these info ['Atlas_Name','Type','Atlas_file_extension',
    'Checkatlas_folder']
    """
    # Create dict with these info ['Atlas_Name','Type','Converted',
    # 'Checkatlas_path']
    clean_atlas_dict = dict()
    for atlas_path in atlas_list:
        atlas_name = get_atlas_name(atlas_path)
        if atlas_path.endswith(".rds"):
            atlas_h5 = atlas_path.replace(".rds", ".h5ad")
            if os.path.exists(atlas_h5):
                print("Seurat file already converted to Scanpy:", atlas_h5)
            else:
                convert_seurat_atlas(atlas_path, atlas_name)
                if os.path.exists(atlas_h5):
                    print(
                        "Include Atlas: " + atlas_name + " from " + atlas_path
                    )
                    info = [
                        atlas_name,
                        "Seurat",
                        ".rds",
                        os.path.dirname(atlas_path) + "/",
                    ]
                    clean_atlas_dict[atlas_h5] = info
        elif atlas_path.endswith(".h5"):
            # detect if its a cellranger output
            if atlas_path.endswith(CELLRANGER_FILE):
                atlas_h5 = atlas_path.replace(CELLRANGER_FILE, "")
                atlas_name = get_atlas_name(atlas_h5)
                print("Include Atlas: " + atlas_name + " from " + atlas_path)
                info = [
                    atlas_name,
                    "Cellranger",
                    ".h5",
                    os.path.dirname(atlas_h5) + "/",
                ]
                clean_atlas_dict[atlas_path] = info
        elif atlas_path.endswith(".h5ad"):
            atlas_rds = atlas_path.replace(".h5ad", ".rds")
            if os.path.exists(atlas_rds):
                print("Include Atlas: " + atlas_name + " from " + atlas_path)
                info = [
                    atlas_name,
                    "Seurat",
                    ".rds",
                    os.path.dirname(atlas_path) + "/",
                ]
                clean_atlas_dict[atlas_path] = info
            else:
                print("Include Atlas: " + atlas_name + " from " + atlas_path)
                info = [
                    atlas_name,
                    "Scanpy",
                    ".h5ad",
                    os.path.dirname(atlas_path) + "/",
                ]
                clean_atlas_dict[atlas_path] = info
    # open file for writing, "w" is writing
    dict_file = open(path + "/checkatlas_files/list_atlases.csv", "w")
    w = csv.writer(dict_file)
    # loop over dictionary keys and values
    for key, val in clean_atlas_dict.items():
        # write every key and value to file
        w.writerow([key, ",".join(val)])
    dict_file.close()
    return clean_atlas_dict


def read_atlas(atlas_path, atlas_info):
    print("--- Load " + atlas_info[0] + " in " + atlas_info[-1])
    try:
        if atlas_path.endswith(".h5"):
            print(atlas_path)
            adata = sc.read_10x_h5(atlas_path)
            adata.var_names_make_unique()
        else:
            adata = sc.read_h5ad(atlas_path)
        return adata
    except anndata._io.utils.AnnDataReadError:
        print("AnnDataReadError, cannot read:", atlas_info[0])
        return None


def summary_table(adata, atlas_path, atlas_name, resume) -> None:
    """
    Extract summary tables from the adata
    :param adata:
    :param atlas_path:
    :param atlas_name:
    :param resume: True if you want to run the function only if the tables
     does not exist already
    :return:
    """
    # create summary table
    csv_path = atlas_path.replace(".h5ad", SUMMARY_EXTENSION)
    print(csv_path)
    if resume:
        if not os.path.isfile(csv_path):
            atlas.create_summary_table(adata, csv_path, atlas_name, path)
    else:
        atlas.create_summary_table(adata, csv_path, atlas_name, path)

    # Create AnnData table
    csv_path = atlas_path.replace(".h5ad", ADATA_EXTENSION)
    if resume:
        if not os.path.isfile(csv_path):
            atlas.create_anndata_table(adata, csv_path, atlas_name, path)
    else:
        atlas.create_anndata_table(adata, csv_path, atlas_name, path)


def start_multithread_client():
    # setup the cluster
    cluster = LocalCluster()
    client = Client(cluster)
    print("Opening Dask browser: http://localhost:8787/status")
    webbrowser.open("http://localhost:8787/status")
    return client


def run(path, atlas_list, multithread, n_cpus):
    """
    Main function of checkatlas
    Run all functions for all atlases:
    - Clean files list by getting only Scanpy atlas (converted from Seurat
    if necessary)
    - Extract summary tables
    - Create UMAP and T-sne figures
    - Calculate every metrics
    Dask powered -- everything is done in parallel

    :param path:
    :param atlas_list:
    :param n_cpus:
    :return:
    """

    folders.checkatlas_folders(path)

    if multithread:
        client = start_multithread_client()
        futures = list()

    # First clean atlas list and keep only the h5ad files
    clean_atlas_dict = clean_list_atlases(atlas_list, path)

    # Create summary files
    for atlas_path, atlas_info in clean_atlas_dict.items():
        if multithread:
            matplotlib.pyplot.switch_backend("Agg")
            # #### Multi-threaded
            # read adata
            atlas_name = atlas_info[0]
            future_name = "Read_" + atlas_name
            adata = client.submit(
                read_atlas, atlas_path, atlas_info, key=future_name
            )
            future_name = "Clean_" + atlas_name
            adata = client.submit(
                atlas.clean_scanpy_atlas, adata, atlas_info, key=future_name
            )
            # Create summary files
            future_name = "Summary_" + atlas_name
            future_sum = client.submit(
                atlas.create_summary_table,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_sum)
            # Adata explorer
            future_name = "Adata_" + atlas_name
            future_sum = client.submit(
                atlas.create_anndata_table,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_sum)
            # Create QC plots
            future_name = "QC_" + atlas_name
            future_sum = client.submit(
                atlas.create_qc_plots,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_sum)
            # Create Umap and TSNE figures
            future_name = "UMAP_" + get_atlas_name(atlas_path)
            future_fig_umap = client.submit(
                atlas.create_umap_fig,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_fig_umap)
            future_name = "TSNE_" + get_atlas_name(atlas_path)
            future_fig_tsne = client.submit(
                atlas.create_tsne_fig,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_fig_tsne)
            # Calculate metrics
            future_name = "Metric_Cluster_" + get_atlas_name(atlas_path)
            future_met = client.submit(
                atlas.metric_cluster,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name,
            )
            futures.append(future_met)
            future_name2 = "Metric_Annot_" + get_atlas_name(atlas_path)
            future_met2 = client.submit(
                atlas.metric_annot,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name2,
            )
            futures.append(future_met2)
            future_name3 = "Metric_DimRed_" + get_atlas_name(atlas_path)
            future_met3 = client.submit(
                atlas.metric_dimred,
                adata,
                atlas_path,
                atlas_info,
                path,
                key=future_name3,
            )
            futures.append(future_met3)
            # Wait for all thread to end
            wait(futures)
        else:
            # ##### Sequential
            atlas_name = atlas_info[0]
            csv_summary_path = os.path.join(
                folders.get_folder(path, folders.SUMMARY),
                atlas_name + SUMMARY_EXTENSION,
            )
            print("Search", csv_summary_path)
            # if not os.path.exists(csv_summary_path):
            adata = read_atlas(atlas_path, atlas_info)
            if adata is not None:
                adata = atlas.clean_scanpy_atlas(adata, atlas_info)
                atlas.create_summary_table(adata, atlas_path, atlas_info, path)
                atlas.create_anndata_table(adata, atlas_path, atlas_info, path)
                atlas.create_qc_plots(adata, atlas_path, atlas_info, path)
                atlas.create_qc_tables(adata, atlas_path, atlas_info, path)
                atlas.create_umap_fig(adata, atlas_path, atlas_info, path)
                atlas.create_tsne_fig(adata, atlas_path, atlas_info, path)
                atlas.metric_cluster(adata, atlas_path, atlas_info, path)
                atlas.metric_annot(adata, atlas_path, atlas_info, path)
                atlas.metric_dimred(adata, atlas_path, atlas_info, path)
            # else:
            #     print("Checkatlas already ran for:", atlas_name)

    print("Run MultiQC")
    multiqc.run_multiqc(path)


if __name__ == "__main__":
    path = "/Users/christophebecavin/Documents/testatlas/"
    # atlas_path = "/Users/christophebecavin/Documents/testatlas/"
    atlas_info = ["test_version", "Scanpy", ".h5ad", "data/test_version.h5ad"]
    # folders.checkatlas_folders(path)
    # atlas_list = list_atlases(path)
    # clean_atlas_dict = clean_list_atlases(atlas_list, path)
    #
    # for atlas_path, atlas_info in clean_atlas_dict.items():
    #     print(atlas_path, atlas_info)
    #     adata = read_atlas(atlas_path, atlas_info)
    #     print(adata is not None)
    #     if adata is not None:
    #         adata = atlas.clean_scanpy_atlas(adata, atlas_info)
    #         atlas.create_summary_table(adata, atlas_path, atlas_info, path)
    #         atlas.create_anndata_table(adata, atlas_path, atlas_info, path)
    #         # atlas.create_qc_plots(adata, atlas_path, atlas_info, path)
    #         atlas.create_umap_fig(adata, atlas_path, atlas_info, path)
    #         # atlas.create_tsne_fig(adata, atlas_path, atlas_info, path)
    #         atlas.metric_cluster(adata, atlas_path, atlas_info, path)
    #         atlas.metric_annot(adata, atlas_path, atlas_info, path)
    #         atlas.metric_dimred(adata, atlas_path, atlas_info, path)
# atlas_path = '/Users/christophebecavin/Documents/testatlas/hca/
# HCA_Barbry_Grch38_Raw_filter_Norm.h5ad'
# atlas_path = '/Users/christophebecavin/Documents/testatlas/Endothelial.h5ad'
# atlas_name = get_atlas_name(atlas_path)
# print('Read atlas')
#     adata = read_atlas(atlas_path, atlas_info)
#     sc.pp.subsample(adata, fraction=0.05)
#     adata.write('/Users/christophebecavin/Documents/testatlas/Endothelial_lite.h5ad')
#     print(adata)
# print('Calc QC')
# figure_path = '/Users/christophebecavin/Documents/testatlas/'
# atlas.create_qc_plots(adata, atlas_path, atlas_name, figure_path)
