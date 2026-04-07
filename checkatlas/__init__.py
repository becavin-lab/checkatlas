from . import atlas, cellranger, check
try:
    from . import seurat
except (ImportError, OSError) as e:
    import logging
    logging.getLogger("checkatlas").warning(
        f"Could not import seurat module (rpy2/R not available): {e}"
    )
    seurat = None

__all__ = ["atlas", "seurat", "check", "cellranger"]
