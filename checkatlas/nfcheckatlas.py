import argparse
import inspect
import logging
import os
from datetime import datetime
from . import checkatlas
from . import atlas, atlas_seurat, checkatlas_workflow, folders, multiqc

logger = logging.getLogger("checkatlas")

def run(args: argparse.Namespace) -> None:
    """
    Main function of checkatlas
    Run all functions for all atlases:
    - Clean files list by getting only Scanpy atlas (converted from Seurat
    if necessary)
    - Extract summary tables
    - Create UMAP and T-sne figures
    - Calculate every metrics

    Args:
        args: List of args for checkatlas program
    Returns:
        None
    """

    logger.debug(f"Transform path to absolute:{args.path}")
    args.path = os.path.abspath(args.path)
    logger.debug(f"Check checkatlas folders in:{args.path}")
    folders.checkatlas_folders(args.path)

    logger.info("Searching Seurat, Cellranger and Scanpy files")
    atlas_list = checkatlas.list_atlases(args.path)
    # First clean atlas list and keep only the h5ad files
    (
        clean_atlas_scanpy,
        clean_atlas_seurat,
        clean_atlas_cellranger,
    ) = checkatlas.clean_list_atlases(atlas_list, args.path)
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

    # Put all atlases together in the list
    clean_atlas = dict(clean_atlas_scanpy)
    clean_atlas.update(clean_atlas_cellranger)
    clean_atlas.update(clean_atlas_seurat)

    if len(clean_atlas_cellranger) > 0:
        logger.debug("Install Seurat if needed")
        atlas_seurat.check_seurat_install()

    # Run all checkatlas analysis
    if args.nextflow == 0:
        logger.info(
            "--nextflow option not found: Run checkatlas workflow "
            "without Nextflow"
        )
        run_checkatlas(clean_atlas, args)
    else:
        clean_atlas.update(clean_atlas_seurat)
        logger.info(
            "--nextflow option found: Run checkatlas workflow with Nextflow"
        )
        logger.info(f"Use {args.nextflow} threads")
        run_checkatlas_nextflow(clean_atlas, args)

    if not args.NOMULTIQC:
        logger.info("Run MultiQC")
        multiqc.run_multiqc(args)


def run_checkatlas_nextflow(clean_atlas, args) -> None:
    """
    Run the checkatlas pipeline by using Nextflow.
    checkatlas_workflow.nf will be run with specific
    parameters.

    Args:
        clean_atlas: List of atlases
        args: List of args for checkatlas program

    Returns:
        None
    """
    checkatlas_workflow.create_checkatlas_worflows(clean_atlas, args)
    script_path = os.path.dirname(os.path.realpath(__file__))
    nextflow_main = os.path.join(script_path, "checkatlas_workflow.nf")
    yaml_files = os.path.join(
        folders.get_folder(args.path, folders.TEMP), "*.yaml"
    )

    # getting the current date and time
    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%Y%d%m-%H%M%S")
    report_file = os.path.join(
        folders.get_workingdir(args.path),
        f"Nextflow_report-{current_time}.html",
    )
    timeline_file = os.path.join(
        folders.get_workingdir(args.path),
        f"Nextflow_timeline-{current_time}.html",
    )
    working_dir_nextflow = folders.get_folder(args.path, folders.NEXTFLOW)
    nextflow_cmd = (
        f"nextflow run -w "
        f"{working_dir_nextflow}"
        f" {nextflow_main} -queue-size {args.nextflow} --files "
        f'"{yaml_files}" -with-report {report_file}'
        f" -with-timeline {timeline_file}"
    )
    logger.debug(f"Execute: {nextflow_cmd}")
    script_path = os.path.dirname(os.path.realpath(__file__))
    nextflow_main = os.path.join(script_path, "checkatlas_workflow.nf")
    # Run Nextflow
    os.system(nextflow_cmd)
    logger.debug(f"Nextflow report saved in {report_file}")
    logger.debug(f"Nextflow timeline saved in {timeline_file}")


def run_checkatlas(clean_atlas, args) -> None:
    """
    Run Checkatlas pipeline for all Scanpy and Cellranger objects
    Args:
        clean_atlas: List of atlas
        args: List of args for checkatlas program
    Returns:
        None
    """

    # List all functions to run
    pipeline_functions_scanpy = checkatlas.get_pipeline_functions(atlas, args)
    pipeline_functions_seurat = checkatlas.get_pipeline_functions(atlas_seurat, args)
    logger.debug(
        f"List of functions which will be ran "
        f"for each Seurat atlas: {pipeline_functions_scanpy}"
    )
    # Go through all atls
    for atlas_path, atlas_info in clean_atlas.items():
        atlas_name = atlas_info[0]
        # Load adata only if resume is not selected
        # and if csv_summary_path do not exist
        csv_summary_path = os.path.join(
            folders.get_folder(args.path, folders.SUMMARY),
            atlas_name + checkatlas.SUMMARY_EXTENSION,
        )
        if args.resume and os.path.exists(csv_summary_path):
            logger.debug(
                f"Skip {atlas_name} summary file already "
                f"exists: {csv_summary_path}"
            )
        else:
            if atlas_info[1] == "Seurat":
                seurat = atlas_seurat.read_atlas(atlas_path, atlas_info)
                logger.info(
                    f"Run checkatlas pipeline for {atlas_name} Seurat atlas"
                )
                # Run pipeline functions
                for function in pipeline_functions_seurat:
                    function(seurat, atlas_path, atlas_info, args)
            else:
                adata = atlas.read_atlas(atlas_path, atlas_info)
                # Clean adata
                adata = atlas.clean_scanpy_atlas(adata, atlas_info)
                logger.info(
                    f"Run checkatlas pipeline for {atlas_name} Scanpy atlas"
                )
                # Run pipeline functions
                for function in pipeline_functions_scanpy:
                    function(adata, atlas_path, atlas_info, args)


if __name__ == "__main__":
    path = "/Users/christophebecavin/Documents/testatlas/"
    # atlas_path = "/Users/christophebecavin/Documents/testatlas/"
    atlas_info = ["test_version", "Scanpy", ".h5ad", "data/test_version.h5ad"]
    # folders.checkatlas_folders(path)
    # atlas_list = list_atlases(path)
