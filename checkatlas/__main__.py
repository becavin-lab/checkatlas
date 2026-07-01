import logging
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(1, os.path.join(sys.path[0], ".."))

try:
    from . import atlas, cellranger, check, seurat
    from .metrics import annot, batch_correction, cluster, dimred
    from .utils import checkatlas_arguments, folders
except ImportError:
    from checkatlas import atlas, cellranger, check, seurat
    from checkatlas.metrics import annot, batch_correction, cluster, dimred
    from checkatlas.utils import checkatlas_arguments, folders


def _metric_nf_arg(metrics: list, all_metrics: list) -> str:
    """Prepare a metric list for the Nextflow --metric_* command-line argument.

    ["none"]  -> "none"  : tells Nextflow to skip this metric category entirely.
    []        -> all metrics joined by spaces (empty list means run everything).
    [m1, m2]  -> "m1 m2" : run only the specified metrics.
    """
    if metrics == ["none"]:
        return "none"
    if not metrics:
        return " ".join(all_metrics)
    return " ".join(metrics)


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

    parser = checkatlas_arguments.create_parser()

    # Parse all args
    args = parser.parse_args()

    if args.process != check.PROCESS_TYPE[0] and args.atlas_name is None:
        parser.error(f"--atlas_name is required for process '{args.process}'")

    # Set logger level
    if args.debug:
        logger.setLevel(getattr(logging, "DEBUG"))
    else:
        logger.setLevel(getattr(logging, "INFO"))

    logger.debug(f"Program arguments: {args}")

    # ── Create checkatlas folder tree before any processing ──────
    folders.checkatlas_folders(args.path)

    process = args.process
    if process == check.PROCESS_TYPE[0]:  # Run pipeline
        args.path = os.path.abspath(args.path)
        logger.info(f"Running Nextflow pipeline on {args.path}")
        nf_script = Path(__file__).resolve().parent / "nextflow" / "main.nf"
        cmd = [
            "nextflow",
            "run",
            str(nf_script),
            f"--path={args.path}",
            f"--outdir={args.path}",
            f"--qc_display={' '.join(args.qc_display)}",
            f"--plot_celllimit={args.plot_celllimit}",
            f"--n_jobs={getattr(args, 'n_jobs', 48)}",
            f"--obs_cluster={' '.join(args.obs_cluster)}",
            f"--metric_cluster={_metric_nf_arg(args.metric_cluster, cluster.__all__)}",
            f"--metric_annot={_metric_nf_arg(args.metric_annot, annot.__all__)}",
            f"--metric_batch_correction={_metric_nf_arg(getattr(args, 'metric_batch_correction', batch_correction.__all__), batch_correction.__all__)}",
            f"--metric_dimred={_metric_nf_arg(args.metric_dimred, dimred.__all__)}",
            f"--checkatlas_debug={'true' if args.debug else 'false'}",
        ]
        if args.with_report is not None:
            cmd += ["-with-report", args.with_report]
        if args.with_dag is not None:
            cmd += ["-with-dag", args.with_dag]
        if args.with_timeline is not None:
            cmd += ["-with-timeline", args.with_timeline]
        logger.debug(f"Command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    else:
        #   ######    Run Checkatlas for a single atlas  #########
        # Look up <atlas_name>.h5ad / .rds / .qs directly at the given path
        atlas_name = args.atlas_name
        base_path = os.path.join(args.path, atlas_name)
        h5ad_path = base_path + ".h5ad"
        rds_path = base_path + ".rds"
        qs_path = base_path + ".qs"
        atlas_info = None

        if os.path.isfile(h5ad_path):
            atlas_info = {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_TYPE_KEY: atlas.ANNDATA_TYPE,
                check.ATLAS_EXTENSION_KEY: ".h5ad",
                check.ATLAS_PATH_KEY: os.path.abspath(h5ad_path),
            }
        elif os.path.isfile(rds_path):
            if seurat is None:
                logger.error("Seurat module not available (rpy2/R not installed)")
                sys.exit(1)
            atlas_info = {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_TYPE_KEY: seurat.SEURAT_TYPE,
                check.ATLAS_EXTENSION_KEY: ".rds",
                check.ATLAS_PATH_KEY: os.path.abspath(rds_path),
            }
        elif os.path.isfile(qs_path):
            if seurat is None:
                logger.error("Seurat module not available (rpy2/R not installed)")
                sys.exit(1)
            atlas_info = {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_TYPE_KEY: seurat.SEURAT_TYPE,
                check.ATLAS_EXTENSION_KEY: ".qs",
                check.ATLAS_PATH_KEY: os.path.abspath(qs_path),
            }
        else:
            logger.error(f"Cannot find {atlas_name}.h5ad, .rds, or .qs in {args.path}")
            sys.exit(1)

        logger.debug(f"Found atlas: {atlas_info}")
        # Run process
        process = args.process
        atlas_type = atlas_info[check.ATLAS_TYPE_KEY]
        if (
            atlas_type == atlas.ANNDATA_TYPE
            or atlas_type == cellranger.CELLRANGER_TYPE_CURRENT
            or atlas_type == cellranger.CELLRANGER_TYPE_OBSOLETE
        ):
            adata = atlas.preprocess_atlas(atlas_info, args)
            if process == check.PROCESS_TYPE[2]:  # summary
                atlas.create_summary_table(adata, atlas_info, args)
                atlas.create_anndata_table(adata, atlas_info, args)
                atlas.create_umap_fig(adata, atlas_info, args)
                atlas.create_tsne_fig(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[3]:  # qc
                atlas.create_qc_tables(adata, atlas_info, args)
                atlas.create_qc_plots(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[4]:  # cluster metrics
                atlas.create_metric_cluster(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[5]:  # annotation metrics
                atlas.create_metric_annot(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[6]:  # batch-correction metrics
                atlas.create_metric_batch_correction(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[7]:  # dimred metrics
                atlas.create_metric_dimred(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[8]:  #  analyse
                atlas.create_summary_table(adata, atlas_info, args)
                atlas.create_anndata_table(adata, atlas_info, args)
                atlas.create_umap_fig(adata, atlas_info, args)
                atlas.create_tsne_fig(adata, atlas_info, args)
                atlas.create_qc_tables(adata, atlas_info, args)
                atlas.create_qc_plots(adata, atlas_info, args)
                atlas.create_metric_cluster(adata, atlas_info, args)
                atlas.create_metric_annot(adata, atlas_info, args)
                atlas.create_metric_batch_correction(adata, atlas_info, args)
                atlas.create_metric_dimred(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[9]:  #  metric
                atlas.create_metric_cluster(adata, atlas_info, args)
                atlas.create_metric_annot(adata, atlas_info, args)
                atlas.create_metric_batch_correction(adata, atlas_info, args)
                atlas.create_metric_dimred(adata, atlas_info, args)
            elif process == check.SCFM_PROCESS_TYPE:  # scfm
                from .scfm.config import from_args
                from .scfm.orchestrator import ensure_per_task_tsvs
                from .scfm.pipeline import run_scfm_from_cache

                scfm_config = from_args(args)
                ensure_per_task_tsvs(adata, atlas_info, args, scfm_config)
                result = run_scfm_from_cache(adata, scfm_config, args=args)
                logger.info(
                    "scfm QC complete for %s. Outputs in %s",
                    atlas_name,
                    result.get("outdir", ""),
                )

        elif seurat is not None and atlas_type == seurat.SEURAT_TYPE:
            if process == check.PROCESS_TYPE[2]:  # summary
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_summary_table(seurat_data, atlas_info, args)
                seurat.create_anndata_table(seurat_data, atlas_info, args)
                seurat.create_umap_fig(seurat_data, atlas_info, args)
                seurat.create_tsne_fig(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[3]:  # qc
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_qc_tables(seurat_data, atlas_info, args)
                seurat.create_qc_plots(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[4]:  # cluster metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_cluster(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[5]:  # annotation metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_annot(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[6]:  # batch-correction metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_batch_correction(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[7]:  # dimred metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_dimred(seurat_data, atlas_info, args)
        else:
            logger.debug("TO DO : Spatial Transcriptomics not yet managed.")


if __name__ == "__main__":
    main()
