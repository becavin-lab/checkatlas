import csv
import inspect
import logging
import os
import webbrowser
import matplotlib
from dask.distributed import Client, LocalCluster, wait

from . import atlas, folders, multiqc, atlas_seurat

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

EXTENSIONS = [".rds", ".h5ad", ".h5"]
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

logger = logging.getLogger("checkatlas")


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


def clean_list_atlases(atlas_list, path) -> tuple:
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
    clean_atlas_scanpy = dict()
    clean_atlas_seurat = dict()
    clean_atlas_cellranger = dict()
    for atlas_path in atlas_list:
        atlas_name = get_atlas_name(atlas_path)
        if atlas_path.endswith(".rds"):
            logger.debug(f"Include Atlas: {atlas_name} from {atlas_path}")
            info = [
                atlas_name,
                "Seurat",
                ".rds",
                os.path.dirname(atlas_path) + "/",
            ]
            clean_atlas_seurat[atlas_path] = info
        elif atlas_path.endswith(".h5"):
            # detect if its a cellranger output
            if atlas_path.endswith(CELLRANGER_FILE):
                atlas_h5 = atlas_path.replace(CELLRANGER_FILE, "")
                atlas_name = get_atlas_name(atlas_h5)
                logger.debug(f"Include Atlas: {atlas_name} from {atlas_path}")
                info = [
                    atlas_name,
                    "Cellranger",
                    ".h5",
                    os.path.dirname(atlas_h5) + "/",
                ]
                clean_atlas_cellranger[atlas_path] = info
        elif atlas_path.endswith(".h5ad"):
            logger.debug(f"Include Atlas: {atlas_name} from {atlas_path}")
            info = [
                atlas_name,
                "Scanpy",
                ".h5ad",
                os.path.dirname(atlas_path) + "/",
            ]
            clean_atlas_scanpy[atlas_path] = info
    # open file for writing, "w" is writing
    dict_file = open(folders.get_workingdir(path) + "list_atlases.csv", "w")
    w = csv.writer(dict_file)
    # loop over dictionary keys and values
    for key, val in clean_atlas_scanpy.items():
        w.writerow([key, ",".join(val)])
    for key, val in clean_atlas_seurat.items():
        w.writerow([key, ",".join(val)])
    for key, val in clean_atlas_cellranger.items():
        w.writerow([key, ",".join(val)])
    dict_file.close()
    return clean_atlas_scanpy, clean_atlas_seurat, clean_atlas_cellranger


def start_multithread_client():
    """
    Open a dask localcluster and a status browser
    :return:
    """
    # setup the cluster
    cluster = LocalCluster()
    client = Client(cluster)
    print("Opening Dask browser: http://localhost:8787/status")
    webbrowser.open("http://localhost:8787/status")
    return client


def get_pipeline_functions(module, args):
    """
    Using arguments of checkatlas program -> build
    the list of functions to run on each adata
    and seurat object
    :param args:
    :return: list of functions to run
    """
    checkatlas_functions = list()
    # Create summary by default, this table is used by the resume option
    checkatlas_functions.append(module.create_summary_table)
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
    return checkatlas_functions


def run(args):
    """
    Main function of checkatlas
    Run all functions for all atlases:
    - Clean files list by getting only Scanpy atlas (converted from Seurat
    if necessary)
    - Extract summary tables
    - Create UMAP and T-sne figures
    - Calculate every metrics
    Dask powered -- everything is done in parallel

    :param args:
    :param logger:
    :return:
    """
    logger.debug(f"Check checkatlas folders in:{args.path}")
    folders.checkatlas_folders(args.path)

    logger.info("Searching Seurat, Cellranger and Scanpy files")
    atlas_list = list_atlases(args.path)
    # First clean atlas list and keep only the h5ad files
    (
        clean_atlas_scanpy,
        clean_atlas_seurat,
        clean_atlas_cellranger,
    ) = clean_list_atlases(atlas_list, args.path)
    logger.info(
        f"Found {len(clean_atlas_scanpy)} potential "
        f"scanpy files with .h5ad extension"
    )
    logger.info(
        f"Found {len(clean_atlas_seurat)} potential "
        f"seurat files with .rds extension"
    )
    logger.info(
        f"Found {len(clean_atlas_cellranger)} cellranger "
        f"file with .h5 extension"
    )

    # Run all checkatlas analysis
    clean_atlas_adata = dict(clean_atlas_scanpy)
    clean_atlas_adata.update(clean_atlas_cellranger)
    run_scanpy(clean_atlas_adata, args)
    run_seurat(clean_atlas_seurat, args)

    if not args.NOMULTIQC:
        logger.info("Run MultiQC")
        multiqc.run_multiqc(args)


def run_scanpy(clean_atlas_scanpy, args):
    """
    Run Checkatlas pipeline for all Scanpy and Cellranger objects
    :param clean_atlas_scanpy:
    :param args:
    :return:
    """
    pipeline_functions = get_pipeline_functions(atlas, args)
    logger.debug(
        f"List of functions which will be ran "
        f"for each atlas: {pipeline_functions}"
    )

    if args.thread > 1:
        client = start_multithread_client()
        futures = list()
        matplotlib.pyplot.switch_backend("Agg")

    # Create summary files
    for atlas_path, atlas_info in clean_atlas_scanpy.items():
        atlas_name = atlas_info[0]

        # Load adata only if resume is not selected
        # and if csv_summary_path do not exist
        csv_summary_path = os.path.join(
            folders.get_folder(args.path, folders.SUMMARY),
            atlas_name + SUMMARY_EXTENSION,
        )
        logger.debug(f"Search {csv_summary_path}")
        if args.resume and os.path.exists(csv_summary_path):
            adata = None
        else:
            adata = atlas.read_atlas(atlas_path, atlas_info)
        if adata is not None:
            # Clean adata
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            logger.info(
                f"Run checkatlas pipeline for {atlas_name} Scanpy atlas"
            )
            # Run pipeline functions
            for function in pipeline_functions:
                if args.thread == 1:
                    function(adata, atlas_path, atlas_info, args)
                else:
                    future_name = function.__name__ + "_" + atlas_name
                    future = client.submit(
                        function,
                        adata,
                        atlas_path,
                        atlas_info,
                        args,
                        key=future_name,
                    )
                    futures.append(future)

            if args.thread > 1:
                # Wait for all thread to end
                wait(futures)


def run_seurat(clean_atlas_seurat, args):
    """
    Run Checkatlas pipeline for all Scanpy objects
    :param clean_atlas_scanpy:
    :param args:
    :return:
    """
    pipeline_functions = get_pipeline_functions(atlas_seurat, args)
    logger.debug(
        f"List of functions which will be ran "
        f"for each Seurat atlas: {pipeline_functions}"
    )

    if args.thread > 1:
        client = start_multithread_client()
        futures = list()
        matplotlib.pyplot.switch_backend("Agg")

    # Create summary files
    for atlas_path, atlas_info in clean_atlas_seurat.items():
        atlas_name = atlas_info[0]

        # Load adata only if resume is not selected
        # and if csv_summary_path do not exist
        csv_summary_path = os.path.join(
            folders.get_folder(args.path, folders.SUMMARY),
            atlas_name + SUMMARY_EXTENSION,
        )
        logger.debug(f"Search {csv_summary_path}")
        if args.resume and os.path.exists(csv_summary_path):
            seurat = None
        else:
            seurat = atlas_seurat.read_atlas(atlas_path, atlas_info)
        if seurat is not None:
            # Clean adata
            # adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            logger.info(
                f"Run checkatlas pipeline for {atlas_name} Seurat atlas"
            )
            # Run pipeline functions
            for function in pipeline_functions:
                if args.thread == 1:
                    function(seurat, atlas_path, atlas_info, args)
                else:
                    future_name = function.__name__ + "_" + atlas_name
                    future = client.submit(
                        function,
                        seurat,
                        atlas_path,
                        atlas_info,
                        args,
                        key=future_name,
                    )
                    futures.append(future)

            if args.thread > 1:
                # Wait for all thread to end
                wait(futures)


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
