import argparse  # pragma: no cover
import logging
import os.path

import yaml

from . import checkatlas  # pragma: no cover
from . import atlas
from . import atlas_seurat
from . import folders

logger = logging.getLogger("checkatlas-workflow")
FUNCTION_TYPE = ['summary', 'qc', 'metric_cluster', 'metric_annot', 'metric_dimred']


def workflow():
    # Set up logging
    logging.basicConfig(format="|--- %(levelname)-8s %(message)s")

    parser = argparse.ArgumentParser(
        prog="checkatlas-workflow",
        usage="checkatlas-workflow summary my_workflow.yaml",
        description="CheckAtlas-workflow is the script ran by checkatlas for each detected "
                    "atlas in Scanpy, Seurat, and CellRanger format.",
        epilog="!!! You are not supposed to run checkatlas-workflow by yourself."
               " It is the main program (checkatlas) who ran it using Nextflow !!!",
    )

    # All Program arguments
    # main_options = parser.add_argument_group("Main arguments")
    parser.add_argument(
        "type",
        type=str,
        help="Provide the type of function to run in the workflow."
             "Different type of function : summary, qc, metric_cluster, metric_annot, metric_dimred",
        default="",
    )

    parser.add_argument(
        "workflow",
        type=str,
        help="Provide a workflow file with all checkatlas-workflow arguments."
             "The workflow config file are created by checkatlas directly"
             "and saved in temp/ folder.",
        default="",
    )
    parser.add_argument(
        "--atlas_name",
        type=str,
        help="Name of the atlas.",
        default=".",
    )
    parser.add_argument(
        "--atlas_type",
        type=str,
        help="Required argument: The type of atlas "
             "Scanpy, CellRanger and Seurat."
             "Values : scanpy, seurat, cellranger",
        default="scanpy",
    )
    parser.add_argument(
        "--atlas_path",
        type=str,
        help="Path of the atlas.",
        default=".",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print out all debug messages.",
    )

    # Pipeline arguments
    pipeline_options = parser.add_argument_group("Manage checkatlas pipeline")
    pipeline_options.add_argument(
        "-na",
        "--NOADATA",
        action="store_true",
        help="Do not Create adata summary tables.",
    )
    pipeline_options.add_argument(
        "-nq",
        "--NOQC",
        action="store_true",
        help="Do not produce any quality control figures or tables.",
    )
    pipeline_options.add_argument(
        "-nr",
        "--NOREDUCTION",
        action="store_true",
        help="Do not Produce UMAP and t-SNE figures.",
    )
    pipeline_options.add_argument(
        "-nm",
        "--NOMETRIC",
        action="store_true",
        help="Do not calculate any metric.",
    )
    pipeline_options.add_argument(
        "-nc", "--NOMULTIQC", action="store_true", help="Do not run multiqc."
    )

    # Arguments linked to QC
    qc_options = parser.add_argument_group("QC options")
    qc_options.add_argument(
        "--qc_display",
        nargs="+",
        type=str,
        default=[
            "violin_plot",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ],
        help="List of QC to display. "
             "Available qc = violin_plot, total_counts, "
             "n_genes_by_counts, pct_counts_mt. "
             "Default: --qc_display violin_plot total_counts "
             "n_genes_by_counts pct_counts_mt",
    )

    # Arguments linked to metric
    metric_options = parser.add_argument_group("Metric options")
    metric_options.add_argument(
        "--obs_cluster",
        nargs="+",
        type=str,
        default=atlas.OBS_CLUSTERS,
        help="List of obs from the adata file to "
             "use in the clustering metric calculus."
             "Example: --obs_cluster celltype leuven seurat_clusters",
    )
    metric_options.add_argument(
        "--metric_cluster",
        nargs="+",
        type=str,
        default=["silhouette", "davies_bouldin"],
        help="List of clustering metrics to calculate."
             "To get a complete list of metrics, look here:"
             "https://github.com/becavin-lab/checkatlas/"
             "blob/main/checkatlas/metrics/metrics.py"
             " Example: --metric_cluster silhouette davies_bouldin",
    )
    metric_options.add_argument(
        "--metric_annot",
        nargs="+",
        type=str,
        default=["rand_index"],
        help="List of clustering metrics to calculate."
             "To get a complete list of metrics, look here:"
             "https://github.com/becavin-lab/checkatlas/blob/"
             "main/checkatlas/metrics/metrics.py"
             " Example: --metric_annot rand_index",
    )
    metric_options.add_argument(
        "--metric_dimred",
        nargs="+",
        type=str,
        default=["kruskal_stress"],
        help="List of dimensionality reduction metrics to calculate."
             "To get a complete list of metrics, look here:"
             "https://github.com/becavin-lab/checkatlas/blob/"
             "main/checkatlas/metrics/metrics.py"
             " Example: --metric_dimred kruskal_stress",
    )

    # Parse all args
    args = parser.parse_args()

    # If a config file was provided, load the new args !
    # Normally a config file is always provided
    if args.workflow != "":
        logger.info(
            "Read workflow config file {} to get new checkatlas configs".format(
                args.workflow
            )
        )
        args = load_arguments(args, args.workflow)

    # Set logger level
    if args.debug:
        logger.setLevel(getattr(logging, "DEBUG"))
    else:
        logger.setLevel(getattr(logging, "INFO"))

    logger.debug(f"Program arguments: {args}")

    if args.atlas_type=='Seurat' :
        pipeline_functions = checkatlas.get_pipeline_functions(atlas_seurat, args)
        run_seurat(args, pipeline_functions)
    else:
        pipeline_functions = checkatlas.get_pipeline_functions(atlas, args)
        run_scanpy(args, pipeline_functions)


def run_scanpy(args, pipeline_functions):
    """
    Run Checkatlas pipeline for Scanpy atlas
    :param args:
    :param pipeline_functions:
    :return:
    """
    # Load adata only if resume is not selected
    # and if csv_summary_path do not exist
    csv_summary_path = os.path.join(
        folders.get_folder(args.path, folders.SUMMARY),
        args.atlas_name + checkatlas.SUMMARY_EXTENSION,
        )
    logger.debug(f"Search {csv_summary_path}")
    ### !!! atlas_info should be removed in the future !
    atlas_info = [args.atlas_name, args.atlas_type, ".rds", args.atlas_path]
    if args.resume and os.path.exists(csv_summary_path):
        adata = None
    else:
        adata = atlas.read_atlas(args.atlas_path, atlas_info)
    if adata is not None:
        # Clean adata
        adata = atlas.clean_scanpy_atlas(adata, atlas_info)
        logger.info(
            f"Run checkatlas pipeline for {args.atlas_name} Scanpy atlas"
        )
        # FUNCTION_TYPE = ['summary', 'qc', 'dimred', 'metric_cluster', 'metric_annot', 'metric_dimred']
        if args.type == "summary":
            if not args.NOADATA:
                atlas.create_anndata_table(adata, args.atlas_path, atlas_info, args)
            if not args.NOREDUCTION:
                atlas.create_umap_fig(adata, args.atlas_path, atlas_info, args)
                atlas.create_tsne_fig(adata, args.atlas_path, atlas_info, args)
            atlas.create_summary_table(adata, args.atlas_path, atlas_info, args)
        elif args.type == "qc" and not args.NOQC:
            if "violin_plot" in args.qc_display:
                atlas.create_qc_plots(adata, args.atlas_path, atlas_info, args)
            if (
                    "total-counts" in args.qc_display
                    or "n_genes_by_counts" in args.qc_display
                    or "pct_counts_mt" in args.qc_display
            ):
                atlas.create_qc_tables(adata, args.atlas_path, atlas_info, args)
        elif args.type == "metric_cluster" and not args.NOMETRIC:
            if len(args.metric_cluster) > 0:
                atlas.metric_cluster(adata, args.atlas_path, atlas_info, args)
            else:
                logger.debug(
                    "No clustering metric was specified in --metric_cluster"
                )
        elif args.type == "metric_annot" and not args.NOMETRIC:
            if len(args.metric_annot) > 0:
                atlas.metric_annot(adata, args.atlas_path, atlas_info, args)
            else:
                logger.debug(
                    "No annotation metric was specified in --metric_annot"
                )
        elif args.type == "metric_annot" and not args.NOMETRIC:
            if len(args.metric_dimred) > 0:
                atlas.metric_dimred(adata, args.atlas_path, atlas_info, args)
        else:
            logger.debug("No dim red metric was specified in --metric_dimred")


def run_seurat(args, pipeline_functions):
    """
    Run Checkatlas pipeline for Scanpy atlas
    :param args:
    :param pipeline_functions:
    :return:
    """
    # Load adata only if resume is not selected
    # and if csv_summary_path do not exist
    csv_summary_path = os.path.join(
        folders.get_folder(args.atlas_path, folders.SUMMARY),
        args.atlas_name + checkatlas.SUMMARY_EXTENSION,
        )
    logger.debug(f"Search {csv_summary_path}")
    ### !!! atlas_info should be removed in the future !
    atlas_info = [args.atlas_name, args.atlas_type, ".rds", args.atlas_path]
    if args.resume and os.path.exists(csv_summary_path):
        seurat = None
    else:
        seurat = atlas_seurat.read_atlas(args.atlas_path, atlas_info)
    if seurat is not None:
        # Clean adata
        # adata = atlas.clean_scanpy_atlas(adata, atlas_info)
        logger.info(
            f"Run checkatlas pipeline for {args.atlas_name} Seurat atlas"
        )
        if args.type == "summary":
            if not args.NOADATA:
                atlas_seurat.create_anndata_table(seurat, args.atlas_path, atlas_info, args)
            if not args.NOREDUCTION:
                atlas_seurat.create_umap_fig(seurat, args.atlas_path, atlas_info, args)
                atlas_seurat.create_tsne_fig(seurat, args.atlas_path, atlas_info, args)
            atlas_seurat.create_summary_table(seurat, args.atlas_path, atlas_info, args)
        elif args.type == "qc" and not args.NOQC:
            if "violin_plot" in args.qc_display:
                atlas_seurat.create_qc_plots(seurat, args.atlas_path, atlas_info, args)
            if (
                    "total-counts" in args.qc_display
                    or "n_genes_by_counts" in args.qc_display
                    or "pct_counts_mt" in args.qc_display
            ):
                atlas_seurat.create_qc_tables(seurat, args.atlas_path, atlas_info, args)
        elif args.type == "metric_cluster" and not args.NOMETRIC:
            if len(args.metric_cluster) > 0:
                atlas_seurat.metric_cluster(seurat, args.atlas_path, atlas_info, args)
            else:
                logger.debug(
                    "No clustering metric was specified in --metric_cluster"
                )
        elif args.type == "metric_annot" and not args.NOMETRIC:
            if len(args.metric_annot) > 0:
                atlas_seurat.metric_annot(seurat, args.atlas_path, atlas_info, args)
            else:
                logger.debug(
                    "No annotation metric was specified in --metric_annot"
                )
        elif args.type == "metric_annot" and not args.NOMETRIC:
            if len(args.metric_dimred) > 0:
                atlas_seurat.metric_dimred(seurat, args.atlas_path, atlas_info, args)
        else:
            logger.debug("No dim red metric was specified in --metric_dimred")


def load_arguments(args, yaml_name):
    """
    Load all args from a yaml file.

    :param args:
    :param yaml_name:
    :return: args
    """
    with open(yaml_name, "r") as config_file:
        yaml_args = yaml.load(config_file, Loader=yaml.FullLoader)
        arg_dict = args.__dict__
        for key, value in yaml_args.items():
            if isinstance(value, list):
                arg_dict[key] = list()
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
        return args


if __name__ == "__main__":  # pragma: no cover
    workflow()