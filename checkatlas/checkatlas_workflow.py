import logging
import os

from checkatlas.utils import files

from . import checkatlas
from .utils import checkatlas_workflow_arguments, folders

PROCESS_TYPE = ["list_scanpy", "list_cellranger", "list_seurat", "html_report"]


def main() -> None:  # pragma: no cover
    """
    The main function executes on commands:
    `python -m checkatlas` and `$ checkatlas `.

    This is checkatlas entry point.

    Arguments are managed here
    Search fo atlases is managed here
    Then checkatlas is ran with the list of atlases found

    Returns:
        None

    """
    # Set up logging
    logger = logging.getLogger("checkatlas")
    logging.basicConfig(format="|--- %(levelname)-8s %(message)s")

    parser = checkatlas_workflow_arguments.create_parser()

    # Parse all args
    args = parser.parse_args()

    # Set logger level
    if args.debug:
        logger.setLevel(getattr(logging, "DEBUG"))
    else:
        logger.setLevel(getattr(logging, "INFO"))

    logger.debug(f"Program arguments: {args}")

    process = args.process
    if process == "list_scanpy":
        logger.debug(f"Search path {args.path} with checkatlas_workflow")
        logger.debug(f"Transform path to absolute:{args.path}")
        args.path = os.path.abspath(args.path)
        logger.debug(f"Check checkatlas folders in:{args.path}")
        folders.checkatlas_folders(args.path)
        logger.info("Searching Seurat, Cellranger and Scanpy files")
        checkatlas.list_scanpy_atlases(args.path)
    elif process == "list_cellranger":
        logger.debug(f"Search path {args.path} with checkatlas_workflow")
        logger.debug(f"Transform path to absolute:{args.path}")
        args.path = os.path.abspath(args.path)
        logger.debug(f"Check checkatlas folders in:{args.path}")
        folders.checkatlas_folders(args.path)
        logger.info("Searching Seurat, Cellranger and Scanpy files")
        checkatlas.list_cellranger_atlases(args.path)
    elif process == "list_seurat":
        logger.debug(f"Search path {args.path} with checkatlas_workflow")
        logger.debug(f"Transform path to absolute:{args.path}")
        args.path = os.path.abspath(args.path)
        logger.debug(f"Check checkatlas folders in:{args.path}")
        folders.checkatlas_folders(args.path)
        logger.info("Searching Seurat, Cellranger and Scanpy files")
        checkatlas.list_seurat_atlases(args.path)
    elif process == "html_report":
        logger.info(
            "Generate QC plots html report in "
            f"{files.get_html_qc_report_path(args.path)}"
        )
        checkatlas.generate_fig_html(args.path, "qc")
        logger.info(
            "Generate UMAP plots html report in "
            f"{files.get_html_umap_report_path(args.path)}"
        )
        checkatlas.generate_fig_html(args.path, "reductions")
    else:
        logger.debug("TO DO : Spatial Transcriptomics not yet managed.")
