import logging
import os

WORKING_DIR = "checkatlas_files"
SUMMARY = "summary"
ANNDATA = "adata"
QC = "qc"
QC_FIG = "violin"
UMAP = "umap"
TSNE = "tsne"
CLUSTER = "cluster"
ANNOTATION = "annotation"
DIMRED = "dimred"
SPECI = "specificity"
TEMP = "temp"
NEXTFLOW = "temp/nextflow"

MULTIQC_FOLDER = "CheckAtlas_MultiQC"

DICT_FOLDER = {
    SUMMARY: SUMMARY,
    ANNDATA: ANNDATA,
    QC: QC,
    QC_FIG: QC_FIG,
    UMAP: UMAP,
    TSNE: TSNE,
    CLUSTER: CLUSTER,
    ANNOTATION: ANNOTATION,
    DIMRED: DIMRED,
    SPECI: SPECI,
    TEMP: TEMP,
}

logger = logging.getLogger("checkatlas")


def get_workingdir(path: str) -> str:
    """Return the working_dir = path of search
    + working_dir
    with working_dir = checkatlas_files/

    Args:
        path (str): Search path for atlas given by user

    Returns:
        str: os.path.join(path, working_dir)
    """
    return os.path.join(path, WORKING_DIR)


def get_folder(path: str, key_folder: str) -> str:
    """Get the folder path giving the search path and
    the folder key in DICT_FOLDER

    Args:
        path (str): Search path for atlas given by user
        key_folder (str): key folder in the DICT_FOLDER
            example: ANNDATA, SUMMARY, UMAP

    Returns:
        str: the folder path
    """
    return os.path.join(get_workingdir(path), DICT_FOLDER[key_folder])


def checkatlas_folders(path: str) -> None:
    """Check in path if the different checkatlas folders exists.<br>
    Create them if needed.
    All folders are given by DICT_FOLDER

    Args:
        path (str): Search path for atlas given by user

    Returns:
        None: None
    """
    global_path = get_workingdir(path)
    if not os.path.exists(global_path):
        os.mkdir(global_path)

    for key_folder in DICT_FOLDER.keys():
        temp_path = os.path.join(global_path, key_folder)
        if not os.path.exists(temp_path):
            logger.debug(f"Create folder: {temp_path}")
            os.mkdir(temp_path)


def get_multiqc_folder(path: str) -> str:
    return os.path.join(path, MULTIQC_FOLDER)
