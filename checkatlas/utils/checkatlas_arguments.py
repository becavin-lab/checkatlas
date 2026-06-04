import argparse
from importlib.resources import files

try:
    from .. import atlas, check
    from ..metrics import annot, cluster, dimred
except ImportError:
    from checkatlas import atlas, check
    from checkatlas.metrics import annot, cluster, dimred


def create_parser():
    parser = argparse.ArgumentParser(
        prog="checkatlas",
        usage="checkatlas [OPTIONS] process path/ [--atlas_name ATLAS_NAME]",
        formatter_class=argparse.RawTextHelpFormatter,
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
        help=(
            "Process to run. Choose from:\n"
            "  run            — launch the full Nextflow pipeline on a folder\n"
            "  preprocess     — load and preprocess a single atlas\n"
            "  summary        — generate summary tables and UMAP/t-SNE figures\n"
            "  qc             — compute QC tables and plots\n"
            "  metric_cluster — compute clustering metrics\n"
            "  metric_annot   — compute annotation metrics\n"
            "  metric_dimred  — compute dimensionality-reduction metrics\n"
            "  metric         — run all metric processes on a single atlas\n"
            "  analyse        — run all analysis processes on a single atlas\n"
        ),
        default="",
    )

    parser.add_argument(
        "path",
        type=str,
        help="Required argument: Your folder containing "
        "Scanpy, CellRanger and Seurat atlases",
        default=".",
    )

    parser.add_argument(
        "--atlas_name",
        type=str,
        help="Required for all processes except 'run': name of the atlas to "
        "process. Atlas_name should be found in one of the samplesheet "
        "provided to nf-checkatlas, or directly created by "
        "checkatlas.list_all_atlases() function",
        default=None,
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
        default=cluster.__all__,
        help="List of clustering metrics to calculate.\n"
        "   By default all available clustering metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_cluster silhouette\n"
        "   Example: --metric_cluster silhouette davies_bouldin\n"
        "   Example: --metric_cluster none\n"
        f"   Available: {cluster.__all__}",
    )
    metric_options.add_argument(
        "--metric_annot",
        nargs="+",
        type=str,
        default=annot.__all__,
        help="List of annotation metrics to calculate.\n"
        "   By default all available annotation metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_annot rand_index\n"
        "   Example: --metric_annot rand_index fowlkes_mallows\n"
        "   Example: --metric_annot none\n"
        f"   Available: {annot.__all__}",
    )
    metric_options.add_argument(
        "--metric_dimred",
        nargs="+",
        type=str,
        default=dimred.__all__,
        help="List of dimensionality reduction metrics to calculate.\n"
        "   By default all available dimensionality reduction metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_dimred kruskal_stress\n"
        "   Example: --metric_dimred kruskal_stress spearman_rho\n"
        "   Example: --metric_dimred none\n"
        f"   Available: {dimred.__all__}",
    )

    # Nextflow reporting options (only used with 'run' process)
    nf_options = parser.add_argument_group(
        "Nextflow reporting options (only used with 'run' process)"
    )
    nf_options.add_argument(
        "--with_report",
        nargs="?",
        const="report.html",
        default=None,
        metavar="FILE",
        help="Generate a Nextflow execution report. "
        "Optionally specify the output filename (default: report.html).",
    )
    nf_options.add_argument(
        "--with_dag",
        nargs="?",
        const="dag.html",
        default=None,
        metavar="FILE",
        help="Generate a Nextflow pipeline DAG. "
        "Optionally specify the output filename (default: dag.html). "
        "Use .dot extension for Graphviz format.",
    )
    nf_options.add_argument(
        "--with_timeline",
        nargs="?",
        const="timeline.html",
        default=None,
        metavar="FILE",
        help="Generate a Nextflow execution timeline. "
        "Optionally specify the output filename (default: timeline.html).",
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
