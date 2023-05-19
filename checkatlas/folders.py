import logging
import os

working_dir = "checkatlas_files"
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
    NEXTFLOW: NEXTFLOW,
}

logger = logging.getLogger("checkatlas")


def get_workingdir(path):
    """
    Return the working_dir = path of search
    + working_dir
    with working_dir = checkatlas_files/
    Args:
        path: Search path for atlas given by user
    Returns:
        str: os.path.join(path, working_dir)
    """
    return os.path.join(path, working_dir)


def get_folder(path, key_folder):
    """
    Get the folder path giving the search path and
    the folder key in DICT_FOLDER
    Args:
        path: Search path for atlas given by user
        key_folder: key folder in the DICT_FOLDER
            example: ANNDATA, SUMMARY, UMAP

    Returns:
        str: the folder path

    """
    return os.path.join(get_workingdir(path), DICT_FOLDER[key_folder])


def checkatlas_folders(path):
    """
    Check in path if the different checkatlas folders exists.<br>
    Create them if needed.
    All folders are given by DICT_FOLDER
    Args:
        path: Search path for atlas given by user
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
