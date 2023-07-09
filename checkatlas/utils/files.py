import os

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
