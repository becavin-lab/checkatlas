import logging

from . import atlas, cellranger, checkatlas, seurat
from .utils import checkatlas_arguments


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
    (
        clean_scanpy_list,
        clean_cellranger_list,
        clean_seurat_list,
    ) = checkatlas.read_list_atlases(args.path)
    clean_scanpy_list = clean_scanpy_list.to_dict("index")
    clean_cellranger_list = clean_cellranger_list.to_dict("index")
    clean_seurat_list = clean_seurat_list.to_dict("index")

    # get atlas_info
    if args.atlas_name in clean_scanpy_list:
        atlas_info = clean_scanpy_list[args.atlas_name]
    elif args.atlas_name in clean_cellranger_list:
        atlas_info = clean_cellranger_list[args.atlas_name]
    elif args.atlas_name in clean_seurat_list:
        atlas_info = clean_seurat_list[args.atlas_name]
    else:
        logger.error(f"Cannot found {args.atlas_name}")
    logger.debug(f"Found atlas: {atlas_info}")
    # Run process
    process = args.process
    atlas_type = atlas_info[checkatlas.ATLAS_TYPE_KEY]
    print(atlas_info)
    if (
        atlas_type == atlas.ANNDATA_TYPE
        or atlas_type == cellranger.CELLRANGER_TYPE_CURRENT
        or atlas_type == cellranger.CELLRANGER_TYPE_OBSOLETE
    ):
        if process == checkatlas.PROCESS_TYPE[0]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_summary_table(adata, atlas_info, args)
            atlas.create_anndata_table(adata, atlas_info, args)
            atlas.create_umap_fig(adata, atlas_info, args)
            atlas.create_tsne_fig(adata, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[1]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_qc_tables(adata, atlas_info, args)
            atlas.create_qc_plots(adata, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[2]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_cluster(adata, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[3]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_annot(adata, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[4]:
            adata = atlas.read_atlas(atlas_info)
            adata = atlas.clean_scanpy_atlas(adata, atlas_info)
            atlas.create_metric_dimred(adata, atlas_info, args)
    elif atlas_type == seurat.SEURAT_TYPE:
        if process == checkatlas.PROCESS_TYPE[0]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_summary_table(seurat_data, atlas_info, args)
            seurat.create_anndata_table(seurat_data, atlas_info, args)
            seurat.create_umap_fig(seurat_data, atlas_info, args)
            seurat.create_tsne_fig(seurat_data, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[1]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_qc_tables(seurat_data, atlas_info, args)
            seurat.create_qc_plots(seurat_data, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[2]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_cluster(seurat_data, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[3]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_annot(seurat_data, atlas_info, args)
        elif process == checkatlas.PROCESS_TYPE[4]:
            seurat_data = seurat.read_atlas(atlas_info)
            seurat.create_metric_dimred(seurat_data, atlas_info, args)
    else:
        logger.debug("TO DO : Spatial Transcriptomics not yet managed.")
