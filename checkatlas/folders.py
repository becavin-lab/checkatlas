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

dict_folder = {
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
}

logger = logging.getLogger("checkatlas")


def get_workingdir(path):
    return os.path.join(path, working_dir)


def get_folder(path, key_folder):
    return os.path.join(get_workingdir(path), dict_folder[key_folder])


def checkatlas_folders(path):
    """
    Check in path if the different checkatlas folders exists.<br>
    Create them if needed.
    :param path:
    :return: None
    """
    global_path = get_workingdir(path)
    if not os.path.exists(global_path):
        os.mkdir(global_path)

    for key_folder in dict_folder.keys():
        temp_path = os.path.join(global_path, key_folder)
        if not os.path.exists(temp_path):
            logger.debug(f"Create folder: {temp_path}")
            os.mkdir(temp_path)
