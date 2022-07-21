import argparse  # pragma: no cover

from . import checkatlas  # pragma: no cover
import logging
import time
import os


def main() -> None:  # pragma: no cover
    """
    The main function executes on commands:
    `python -m checkatlas` and `$ checkatlas `.

    This is checkatlas entry point.

    Arguments are managed here
    Search fo atlases is managed here
    Then checkatlas is ran with the list of atlases found
    """

    # Set up logging
    start_execution_time = time.time()
    logger = logging.getLogger("checkatlas")
    logging.basicConfig(format='%(name)-12s: %(levelname)-8s %(message)s')

    parser = argparse.ArgumentParser(
        prog="checkatlas",
        usage="checkatlas [OPTIONS] your_search_folder/",
        description="CheckAtlas is a one liner tool to check the quality of your \n"
                    "single-cell atlases. For every atlas, it produces the quality "
                    "control tables and figures which can be then processed by multiqc. "
                    "CheckAtlas is able to load Scanpy, Seurat, and CellRanger files.",
        epilog="Enjoy the checkatlas functionality!",
    )
    
    # All Program arguments
    # main_options = parser.add_argument_group("Main arguments")
    parser.add_argument(
        "path",
        type=str,
        help="Path containing Scanpy, CellRanger and Seurat atlases",
        default=".",
    )
    parser.add_argument(
        "-m",
        "--multiqc",
        type=str,
        help="Set Multiqc out folder. Default: CheckAtlas_MultiQC",
        default="CheckAtlas_MultiQC",
    )
    parser.add_argument(
        "-t",
        "--thread",
        default=1,
        help="Number of threads for parallel computing. \nIf -t > 1,"
             " this will activate automatically the parallel computing"
             " mode using Dask.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print out all debug messages",
    )

    # Pipeline arguments
    pipeline_options = parser.add_argument_group("Manage checkatlas pipeline")
    pipeline_options.add_argument("-ns", "--NOSUMMARY", action="store_true",
                        help="Do not create general summary tables")
    pipeline_options.add_argument("-na", "--NOADATA", action="store_true",
                                  help="Do not Create adata summary tables")
    pipeline_options.add_argument("-nq", "--NOQC", action="store_true",
                                  help="Do not produce any quality control figures or tables")
    pipeline_options.add_argument("-nr", "--NOREDUCTION", action="store_true",
                              help="Do not Produce UMAP and t-SNE figures")
    pipeline_options.add_argument("-nm", "--NOMETRIC", action="store_true",
                                  help="Do not calculate any metric.")
    pipeline_options.add_argument("-nc", "--NOMULTIQC", action="store_true",
                              help="Do not run multiqc.")

    # Arguments linked to QC
    qc_options = parser.add_argument_group("QC options")
    qc_options.add_argument("--qc_display", nargs='+', type=str,
                            default=["violin_plot","total-counts",
                                    "n_genes_by_counts","pct_counts_mt"],
                            help="List of QC to display. "
                                 "Available qc = violin_plot, total_counts, "
                                 "n_genes_by_counts, pct_counts_mt. "
                                 "Default: --qc_display violin_plot total_counts "
                                 "n_genes_by_counts pct_counts_mt")

    # Arguments linked to metric
    metric_options = parser.add_argument_group("Metric options")



    # Parse all args
    args = parser.parse_args()
    print("Get arguments of the program:",args)

    # Set logger level
    if args.debug:
        logger.setLevel(getattr(logging, "DEBUG"))
    else:
        logger.setLevel(getattr(logging, "INFO"))

    script_path = os.path.dirname(os.path.realpath(__file__))
    logger.debug("Path_script", script_path)

    checkatlas.run(args, logger)


if __name__ == "__main__":  # pragma: no cover
    main()