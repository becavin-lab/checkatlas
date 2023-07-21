import argparse
from importlib.resources import files

from checkatlas import atlas, checkatlas
from checkatlas.metrics import annot, cluster, dimred


def create_parser():
    parser = argparse.ArgumentParser(
        prog="checkatlas",
        usage="checkatlas [OPTIONS] process atlas_name your_search_folder/",
        description="CheckAtlas is a one liner tool to check the "
        "quality of your single-cell atlases. For "
        "every atlas, it produces the quality control "
        "tables and figures which can be then processed "
        "by multiqc. CheckAtlas is able to load Scanpy, "
        "Seurat, and CellRanger files.",
        epilog="Enjoy the checkatlas functionality!",
    )

    # All Program arguments
    # All Program arguments
    # main_options = parser.add_argument_group("Main arguments")
    parser.add_argument(
        "process",
        type=str,
        help="Required argument: Type of process to run"
        f" among {checkatlas.PROCESS_TYPE}",
        default="",
    )

    parser.add_argument(
        "atlas_name",
        type=str,
        help="Required argument: The name of the atlas to process."
        "Atlas_name should be found in one of the samplesheet provided to"
        "nf-checkatlas, or directly created by checkatlas.list_all_atlases()"
        " function",
        default=".",
    )

    parser.add_argument(
        "path",
        type=str,
        help="Required argument: Your folder containing "
        "Scanpy, CellRanger and Seurat atlasesv",
        default=".",
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
        version=f"Checkatlas, version {get_version()}",
        help="Display checkatlas version.",
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
    qc_options.add_argument(
        "--plot_celllimit",
        type=int,
        default=0,
        help="Set the maximum number of cells"
        "to plot in QC, UMAP, t-SNE, etc...."
        "If plot_celllimit=0, no limit will"
        "be applied.",
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
        default=["davies_bouldin"],
        # default=["silhouette", "davies_bouldin"],
        # default=[""],
        help="Specify the list of clustering metrics to calculate.\n"
        "   Example: --metric_cluster silhouette davies_bouldin\n"
        f"   List of cluster metrics: {cluster.__all__}",
    )
    metric_options.add_argument(
        "--metric_annot",
        nargs="+",
        type=str,
        # default=[],
        default=["rand_index"],
        help=f"Specify the list of clustering metrics to calculate."
        f"   Example: --metric_annot rand_index"
        f"   List of annotation metrics: {annot.__all__}",
    )
    metric_options.add_argument(
        "--metric_dimred",
        nargs="+",
        type=str,
        default=[],
        # default=["kruskal_stress"],
        help="Specify the list of dimensionality reduction "
        "metrics to calculate.\n"
        "   Example: --metric_dimred kruskal_stress\n"
        f"   List of dim. red. metrics: {dimred.__all__}",
    )
    return parser


def get_version():
    """
    Get version of checkatlas from checkatlas/VERSION file
    :return: checkatlas version
    """
    version_file = files(__package__).joinpath("VERSION")
    with open(version_file) as file:
        version = file.readline()
        return version
