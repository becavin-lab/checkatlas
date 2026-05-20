import logging
import os
import sys

sys.path.insert(1, os.path.join(sys.path[0], ".."))

try:
    from . import atlas, cellranger, check, seurat
    from .utils import checkatlas_arguments
except ImportError:
    from checkatlas import atlas, cellranger, check, seurat
    from checkatlas.utils import checkatlas_arguments


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

    #   ######    Run Checkatlas   #########
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
        if process == check.PROCESS_TYPE[0]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_summary_table(adata, atlas_info, args)
            atlas.create_anndata_table(adata, atlas_info, args)
            atlas.create_umap_fig(adata, atlas_info, args)
            atlas.create_tsne_fig(adata, atlas_info, args)
        elif process == check.PROCESS_TYPE[1]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_qc_tables(adata, atlas_info, args)
            atlas.create_qc_plots(adata, atlas_info, args)
        elif process == check.PROCESS_TYPE[2]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_cluster(adata, atlas_info, args)
        elif process == check.PROCESS_TYPE[3]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_annot(adata, atlas_info, args)
        elif process == check.PROCESS_TYPE[4]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_dimred(adata, atlas_info, args)
    elif atlas_type == seurat.SEURAT_TYPE:
        if process == check.PROCESS_TYPE[0]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_summary_table(seurat_data, atlas_info, args)
            seurat.create_anndata_table(seurat_data, atlas_info, args)
            seurat.create_umap_fig(seurat_data, atlas_info, args)
            seurat.create_tsne_fig(seurat_data, atlas_info, args)
        elif process == check.PROCESS_TYPE[1]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_qc_tables(seurat_data, atlas_info, args)
            seurat.create_qc_plots(seurat_data, atlas_info, args)
        elif process == check.PROCESS_TYPE[2]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_cluster(seurat_data, atlas_info, args)
        elif process == check.PROCESS_TYPE[3]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_annot(seurat_data, atlas_info, args)
        elif process == check.PROCESS_TYPE[4]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_dimred(seurat_data, atlas_info, args)
    else:
        logger.debug("TO DO : Spatial Transcriptomics not yet managed.")


if __name__ == "__main__":
    main()
