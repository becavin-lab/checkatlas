from types import ModuleType
from typing import Optional

from . import atlas, cellranger, check

seurat: Optional[ModuleType] = None
try:
    from . import seurat
except (ImportError, OSError) as e:
    import logging

    logging.getLogger("checkatlas").warning(
        f"Could not import seurat module (rpy2/R not available): {e}"
    )

__all__ = ["atlas", "seurat", "check", "cellranger"]
