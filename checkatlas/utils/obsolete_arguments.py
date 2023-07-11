import argparse

import yaml

from .. import atlas
from ..metrics import annot, cluster, dimred


def create_parser():
    parser = argparse.ArgumentParser(
        prog="checkatlas",
        usage="checkatlas [OPTIONS] your_search_folder/",
        description="CheckAtlas is a one liner tool to check the "
        "quality of your single-cell atlases. For "
        "every atlas, it produces the quality control "
        "tables and figures which can be then processed "
        "by multiqc. CheckAtlas is able to load Scanpy, "
        "Seurat, and CellRanger files.",
        epilog="Enjoy the checkatlas functionality!",
    )

    # All Program arguments
    # main_options = parser.add_argument_group("Main arguments")
    parser.add_argument(
        "path",
        type=str,
        help="Required argument: Your folder containing "
        "Scanpy, CellRanger and Seurat atlasesv",
        default=".",
    )

    parser.add_argument(
        "--plot_celllimit",
        type=int,
        default=10000,
        help="Set the maximum number of cells"
        "to plot in QC, UMAP, t-SNE, etc...."
        "If plot_celllimit=0, no limit will"
        "be applied.",
    )

    parser.add_argument(
        "--nextflow",
        type=int,
        default=0,
        help="Activate Nextflow and specify the number of threads to use. \n"
        "Example: --nextflow 8 (for 8 threads)",
    )

    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print out all debug messages.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Checkatlas {get_version()}",
        help="Display checkatlas version.",
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
        # default=["silhouette", "davies_bouldin"],
        default=["davies_bouldin"],
        help="Specify the list of clustering metrics to calculate.\n"
        "   Example: --metric_cluster silhouette davies_bouldin\n"
        f"   List of cluster metrics: {cluster.__all__}",
    )
    metric_options.add_argument(
        "--metric_annot",
        nargs="+",
        type=str,
        default=[],
        # default=["rand_index"],
        help=f"Specify the list of clustering metrics to calculate."
        f"   Example: --metric_annot rand_index"
        f"   List of annotation metrics: {annot.__all__}",
    )
    metric_options.add_argument(
        "--metric_dimred",
        nargs="+",
        type=str,
        # default=["kruskal_stress"],
        default=[],
        help="Specify the list of dimensionality reduction "
        "metrics to calculate.\n"
        "   Example: --metric_dimred kruskal_stress\n"
        f"   List of dim. red. metrics: {dimred.__all__}",
    )
    metric_options.add_argument(
        "--TEST_ALLMETRICS",
        action="store_true",
        help="Run the pipeline with all metrics available.\n"
        f"   List of cluster metrics: {cluster.__all__}\n"
        f"   List of annotation metrics: {annot.__all__}\n"
        f"   List of dim. red. metrics: {dimred.__all__}",
    )
    return parser


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


def save_arguments(args, yaml_name):
    """
    Save all args to a yaml file. Only use this
    function to create example yaml config files.

    :param args:
    :param yaml_name:
    :return:
    """
    with open(yaml_name, "w") as config_file:
        yaml.dump(args.__dict__, config_file)


def get_version():
    """
    Get version of checkatlas from checkatlas/VERSION file
    :return: checkatlas version
    """
    return "Need to Fix version reading!"
    # script_path = os.path.dirname(os.path.realpath(__file__))
    # version_file = os.path.join(script_path, "VERSION")
    # with open(version_file, "r") as version:
    #     return version.readlines()[0].strip()
