import os
from . import folders
from checkatlas import checkatlas
import argparse

def get_summary_table_path(atlas_name:str, path: str) -> str:
    """_summary_

    Args:
        atlas_name (str): _description_
        args (argparse.Namespace): _description_

    Returns:
        str: _description_
    """    
    csv_path = os.path.join(
        folders.get_folder(path, folders.SUMMARY),
        atlas_name + checkatlas.SUMMARY_EXTENSION,
    )
    return csv_path


def get_anndata_table_path(atlas_name:str, path: str) -> str:
    """_summary_

    Args:
        atlas_name (str): _description_
        args (argparse.Namespace): _description_

    Returns:
        str: _description_
    """    
    csv_path = os.path.join(
        folders.get_folder(path, folders.ANNDATA),
        atlas_name + checkatlas.ADATA_EXTENSION,
    )
    return csv_path