import os

import pandas as pd

from checkatlas import checkatlas

from . import folders


def get_file_path(
    atlas_name: str, folder: str, extension: str, path: str
) -> str:
    """_summary_

    Args:
        atlas_name (str): _description_
        args (argparse.Namespace): _description_

    Returns:
        str: _description_
    """
    csv_path = os.path.join(
        folders.get_folder(path, folder),
        atlas_name + extension,
    )
    return csv_path


def get_table_scanpy_path(checkatlas_path: str) -> str:
    file_path = os.path.join(
        folders.get_workingdir(checkatlas_path), "List_scanpy.csv"
    )
    return file_path


def get_table_cellranger_path(checkatlas_path: str) -> str:
    file_path = os.path.join(
        folders.get_workingdir(checkatlas_path), "List_cellranger.csv"
    )
    return file_path


def get_table_seurat_path(checkatlas_path: str) -> str:
    file_path = os.path.join(
        folders.get_workingdir(checkatlas_path), "List_seurat.csv"
    )
    return file_path


def save_list_scanpy(clean_scanpy_list: list, checkatlas_path: str) -> None:
    df_summary = pd.DataFrame(columns=checkatlas.ATLAS_TABLE_HEADER)
    for table_info in clean_scanpy_list:
        df_summary.loc[
            table_info[checkatlas.ATLAS_NAME_KEY]
        ] = table_info.values()
    df_summary.to_csv(
        get_table_scanpy_path(checkatlas_path), index=False, sep=","
    )


def save_list_cellranger(
    clean_cellranger_list: list, checkatlas_path: str
) -> None:
    df_summary = pd.DataFrame(columns=checkatlas.ATLAS_TABLE_HEADER)
    for table_info in clean_cellranger_list:
        df_summary.loc[
            table_info[checkatlas.ATLAS_NAME_KEY]
        ] = table_info.values()
    df_summary.to_csv(
        get_table_cellranger_path(checkatlas_path), index=False, sep=","
    )


def save_list_seurat(clean_seurat_list: list, checkatlas_path: str) -> None:
    df_summary = pd.DataFrame(columns=checkatlas.ATLAS_TABLE_HEADER)
    for table_info in clean_seurat_list:
        df_summary.loc[
            table_info[checkatlas.ATLAS_NAME_KEY]
        ] = table_info.values()
    df_summary.to_csv(
        get_table_seurat_path(checkatlas_path), index=False, sep=","
    )


def get_html_qc_report_path(checkatlas_path: str):
    return os.path.join(
        folders.get_workingdir(checkatlas_path),
        "Checkatlas_QC_report_mqc.html",
    )


def get_html_umap_report_path(checkatlas_path: str):
    return os.path.join(
        folders.get_workingdir(checkatlas_path),
        "Checkatlas_UMAP_report_mqc.html",
    )


def get_html_tsne_report_path(checkatlas_path: str):
    return os.path.join(
        folders.get_workingdir(checkatlas_path),
        "Checkatlas_tSNE_report_mqc.html",
    )
