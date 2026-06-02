import logging
import os
import subprocess
import sys
from pathlib import Path

# Patch Python 3.13 multiprocessing resource_tracker to handle joblib/loky
# resource types (file, folder, semlock) that the tracker doesn't recognize.
# Without this, cleanup floods stderr with hundreds of ValueError messages.
def _patch_resource_tracker():
    try:
        import multiprocessing.resource_tracker as _rt
        import multiprocessing.resource_tracker as _mp_rt

        # Build a replacement main() that handles unknown rtypes gracefully
        _KNOWN = frozenset(("semlock", "file", "folder", "shared_memory"))

        def _safe_main():
            import atexit
            import os
            import signal
            import sys
            import threading
            import warnings
            from multiprocessing.util import info, warning

            keep_running = True
            locked_files = set()

            def _cleanup(signum=-1, frame=None):
                nonlocal keep_running
                keep_running = False

            signal.signal(signal.SIGINT, _cleanup)
            signal.signal(signal.SIGTERM, _cleanup)

            # Read fd from argv (copied from original main)
            fd = int(sys.argv[1]) if len(sys.argv) > 1 else None
            if fd is None:
                return

            while keep_running:
                try:
                    # Read a line from the pipe
                    line = os.read(fd, 300)
                    if not line:
                        break
                    line = line.decode("ascii", errors="replace").strip()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    cmd, name, rtype = parts[0], parts[1], parts[2]

                    if rtype not in _KNOWN:
                        continue  # silently skip unknown types

                    if cmd == "REGISTER":
                        locked_files.add(name)
                    elif cmd == "UNREGISTER":
                        try:
                            locked_files.remove(name)
                        except KeyError:
                            pass
                except (OSError, EOFError):
                    break

            # Cleanup remaining
            for name in list(locked_files):
                try:
                    if os.path.exists(name):
                        if os.path.isdir(name):
                            os.rmdir(name)
                        else:
                            os.unlink(name)
                except OSError:
                    pass

        _rt.main = _safe_main
    except Exception:
        pass

_patch_resource_tracker()

sys.path.insert(1, os.path.join(sys.path[0], ".."))

try:
    from . import atlas, cellranger, check, seurat
    from .utils import checkatlas_arguments, folders
except ImportError:
    from checkatlas import atlas, cellranger, check, seurat
    from checkatlas.utils import checkatlas_arguments, folders


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

    # Set logger level
    if args.debug:
        logger.setLevel(getattr(logging, "DEBUG"))
    else:
        logger.setLevel(getattr(logging, "INFO"))

    logger.debug(f"Program arguments: {args}")

    # ── Create checkatlas folder tree before any processing ──────
    folders.checkatlas_folders(args.path)

    process = args.process
    if process == check.PROCESS_TYPE[0]: # Run pipeline
        args.path = os.path.abspath(args.path)
        logger.info(f"Running Nextflow pipeline on {args.path}")
        nf_script = Path(__file__).resolve().parent.parent / "main.nf"
        cmd = [
            "nextflow", "run", str(nf_script),
            f"--path={args.path}",
            f"--outdir={args.path}",
        ]
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
            atlas_info = {
                check.ATLAS_NAME_KEY: atlas_name,
                check.ATLAS_TYPE_KEY: seurat.SEURAT_TYPE,
                check.ATLAS_EXTENSION_KEY: ".rds",
                check.ATLAS_PATH_KEY: os.path.abspath(rds_path),
            }
        elif os.path.isfile(qs_path):
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
            adata = atlas.preprocess_atlas(atlas_info)
            if process == check.PROCESS_TYPE[2]: # summary
                atlas.create_summary_table(adata, atlas_info, args)
                atlas.create_anndata_table(adata, atlas_info, args)
                atlas.create_umap_fig(adata, atlas_info, args)
                atlas.create_tsne_fig(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[3]:  # qc  
                atlas.create_qc_tables(adata, atlas_info, args)
                atlas.create_qc_plots(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[4]: # cluster metrics
                atlas.create_metric_cluster(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[5]: # annotation metrics
                atlas.create_metric_annot(adata, atlas_info, args)
            elif process == check.PROCESS_TYPE[6]: # dimred metrics
                atlas.create_metric_dimred(adata, atlas_info, args)
        elif atlas_type == seurat.SEURAT_TYPE:
            if process == check.PROCESS_TYPE[2]: # summary 
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_summary_table(seurat_data, atlas_info, args)
                seurat.create_anndata_table(seurat_data, atlas_info, args)
                seurat.create_umap_fig(seurat_data, atlas_info, args)
                seurat.create_tsne_fig(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[3]:  # qc
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_qc_tables(seurat_data, atlas_info, args)
                seurat.create_qc_plots(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[4]: # cluster metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_cluster(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[5]: # annotation metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_annot(seurat_data, atlas_info, args)
            elif process == check.PROCESS_TYPE[6]: # dimred metrics
                seurat_data = seurat.read_atlas(atlas_info)
                seurat.create_metric_dimred(seurat_data, atlas_info, args)
        else:
            logger.debug("TO DO : Spatial Transcriptomics not yet managed.")


if __name__ == "__main__":
    main()
