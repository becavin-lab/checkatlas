import logging
import os

import pandas as pd

from . import atlas, cellranger, seurat
from .utils import files as chk_files
from .utils import folders

"""
checkatlas base module.
This is the principal module of the checkatlas project.

"""

PROCESS_TYPE = [
    "summary",
    "qc",
    "metric_cluster",
    "metric_annot",
    "metric_dimred",
]

TSV_EXTENSION = ".tsv"
QC_FIG_EXTENSION = "_qc.png"
UMAP_EXTENSION = "_umap.png"
TSNE_EXTENSION = "_tsne.png"

ATLAS_NAME_KEY = "Atlas_name"
ATLAS_TYPE_KEY = "Atlas_type"
ATLAS_EXTENSION_KEY = "Atlas_extension"
ATLAS_PATH_KEY = "Atlas_path"
ATLAS_TABLE_HEADER = [
    ATLAS_NAME_KEY,
    ATLAS_TYPE_KEY,
    ATLAS_EXTENSION_KEY,
    ATLAS_PATH_KEY,
]

openpng_html_script = """
<script>
function openPNG(evt, pngName, tablinks_id, tabcontent_id) {
  var i, tabcontent, tablinks;
  tabcontent = document.getElementsByClassName(tabcontent_id);
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }
  tablinks = document.getElementsByClassName(tablinks_id);
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }
  document.getElementById(pngName).style.display = "block";
  evt.currentTarget.className += " active";
}
</script>
"""

logger = logging.getLogger("checkatlas")


def list_scanpy_atlases(checkatlas_path: str) -> None:
    # Get all files with matching extension
    atlas_list = list()
    for root, dirs, files in os.walk(checkatlas_path):
        for file in files:
            if file.endswith(atlas.ANNDATA_EXTENSION):
                atlas_list.append(os.path.join(root, file))

    # Filter the lists keepng only atlases
    clean_scanpy_list = list()
    for atlas_path in atlas_list:
        atlas_info = atlas.detect_scanpy(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_scanpy_list.append(atlas_info)
    logger.info(
        f"Found {len(clean_scanpy_list)} potential "
        f"seurat files with .rds extension"
    )
    # Save the list of atlas taken into account
    chk_files.save_list_scanpy(clean_scanpy_list, checkatlas_path)


def list_cellranger_atlases(checkatlas_path: str) -> None:
    # Get all files with matching extension
    EXTENSIONS = [
        cellranger.CELLRANGER_FILE,
        cellranger.CELLRANGER_MATRIX_FILE,
    ]
    atlas_list = list()
    for root, dirs, files in os.walk(checkatlas_path):
        for file in files:
            for extension in EXTENSIONS:
                if file.endswith(extension):
                    atlas_list.append(os.path.join(root, file))

    # Filter the lists keepng only atlases
    clean_cellranger_list = list()
    for atlas_path in atlas_list:
        atlas_info = cellranger.detect_cellranger(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a cellranger output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_cellranger_list.append(atlas_info)
    logger.info(
        f"Found {len(clean_cellranger_list)} potential "
        f"seurat files with .rds extension"
    )
    # Save the list of atlas taken into account
    chk_files.save_list_cellranger(clean_cellranger_list, checkatlas_path)


def list_seurat_atlases(checkatlas_path: str) -> None:
    # Get all files with matching extension
    atlas_list = list()
    for root, dirs, files in os.walk(checkatlas_path):
        for file in files:
            if file.endswith(seurat.SEURAT_EXTENSION):
                atlas_list.append(os.path.join(root, file))

    # Filter the lists keepng only atlases
    clean_seurat_list = list()
    for atlas_path in atlas_list:
        atlas_info = seurat.detect_seurat(atlas_path)
        if len(atlas_info) != 0:
            # detect if its a seurat output
            logger.debug(
                f"Include Atlas: {atlas_info[ATLAS_NAME_KEY]}"
                f" of type {atlas_info[ATLAS_TYPE_KEY]}"
                f"from {atlas_info[ATLAS_PATH_KEY]}"
            )
            clean_seurat_list.append(atlas_info)
    logger.info(
        f"Found {len(clean_seurat_list)} potential "
        f"seurat files with .rds extension"
    )
    # Save the list of atlas taken into account
    chk_files.save_list_seurat(clean_seurat_list, checkatlas_path)


def read_list_atlases(checkatlas_path: str) -> tuple:
    clean_scanpy_list = pd.DataFrame()
    clean_cellranger_list = pd.DataFrame()
    clean_seurat_list = pd.DataFrame()
    if os.path.exists(chk_files.get_table_scanpy_path(checkatlas_path)):
        clean_scanpy_list = pd.read_csv(
            chk_files.get_table_scanpy_path(checkatlas_path)
        )
        clean_scanpy_list.index = clean_scanpy_list[ATLAS_NAME_KEY]
    if os.path.exists(chk_files.get_table_cellranger_path(checkatlas_path)):
        clean_cellranger_list = pd.read_csv(
            chk_files.get_table_cellranger_path(checkatlas_path)
        )
        clean_cellranger_list.index = clean_cellranger_list[ATLAS_NAME_KEY]
    if os.path.exists(chk_files.get_table_seurat_path(checkatlas_path)):
        clean_seurat_list = pd.read_csv(
            chk_files.get_table_seurat_path(checkatlas_path)
        )
        clean_seurat_list.index = clean_seurat_list[ATLAS_NAME_KEY]
    return clean_scanpy_list, clean_cellranger_list, clean_seurat_list


def generate_fig_html(checkatlas_path: str, type_viz: str):
    logger.info("Generate html reports with all figs for {type}")
    if type_viz == "qc":
        # data is a dictionnary wth key = data_name and value = fig_path
        qc_dict = dict()
        qc_fig_path = folders.get_folder(checkatlas_path, folders.QC_FIG)
        for qc_fig in os.listdir(qc_fig_path):
            if qc_fig.endswith(".png"):
                atlas_name = os.path.splitext(os.path.basename(qc_fig))[0]
                qc_dict[atlas_name] = os.path.join("violin", qc_fig)
        if len(qc_dict) > 0:
            html = create_img_html_content(type_viz, qc_dict)
            with open(
                chk_files.get_html_qc_report_path(checkatlas_path), "w"
            ) as html_report:
                html_report.write(html)
    elif type_viz == "reductions":
        umap_dict = dict()
        umap_fig_path = folders.get_folder(checkatlas_path, folders.UMAP)
        for umap_fig in os.listdir(umap_fig_path):
            if umap_fig.endswith(".png"):
                atlas_name = os.path.splitext(os.path.basename(umap_fig))[0]
                umap_dict[atlas_name] = os.path.join("umap", umap_fig)
        if len(umap_dict) > 0:
            html = create_img_html_content(type_viz, umap_dict)
            with open(
                chk_files.get_html_umap_report_path(checkatlas_path), "w"
            ) as html_report:
                html_report.write(html)

        tsne_dict = dict()
        tsne_fig_path = folders.get_folder(checkatlas_path, folders.TSNE)
        for tsne_fig in os.listdir(tsne_fig_path):
            if tsne_fig.endswith(".png"):
                atlas_name = os.path.splitext(os.path.basename(tsne_fig))[0]
                tsne_dict[atlas_name] = os.path.join("tsne", tsne_fig)
        if len(tsne_dict) > 0:
            html = create_img_html_content(type_viz, tsne_dict)
            with open(
                chk_files.get_html_tsne_report_path(checkatlas_path), "w"
            ) as html_report:
                html_report.write(html)
    else:
        logger.error(f"Type of vsualization not recognized {type_viz}")


def create_img_html_content(type_viz, data):
    html_content = f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>Custom content
            with all {type_viz} plots found
            in your atlases</head>
            <body><br><br>
            <div class="tab">\n"""
    counter = 0
    tablinks = type_viz + "_tablinks"
    tabcontent = type_viz + "_tabcontent"
    for atlas_name, fig_path in data.items():
        class_tab = tablinks
        if counter == 0:
            class_tab = tablinks + " active"
        html_content += add_selection_img_button(
            class_tab, type_viz, atlas_name, tablinks, tabcontent
        )
        counter += 1
    html_content += """</div>"""

    counter = 0
    for atlas_name, fig_path in data.items():
        style = "display: none;"
        if counter == 0:
            style = "display: block;"
        html_content += add_div_img(
            fig_path, type_viz, atlas_name, style, tabcontent
        )
        counter += 1

    html_content += openpng_html_script
    html_content += """
        </body>
        </html>
        """

    return html_content


def add_selection_img_button(
    class_tab, type_viz, atlas_name, tablinks, tabcontent
):
    return (
        """  <button class=\""""
        + class_tab
        + """\" onclick="openPNG(event, '"""
        + type_viz
        + "_"
        + atlas_name
        + """',
                            '"""
        + tablinks
        + "','"
        + tabcontent
        + """')">"""
        + atlas_name
        + """</button>\n"""
    )


def add_div_img(path_fig, type_viz, atlas_name, style, tabcontent):
    return (
        """<div id=\""""
        + type_viz
        + "_"
        + atlas_name
        + """\" class=\""""
        + tabcontent
        + """\" style =\""""
        + style
        + """\"><img src=\""""
        + path_fig
        + """\" >\n</div>\n"""
    )


if __name__ == "__main__":
    path = "/data/analysis/data_becavin/checkatlas_test/tuto"
    atlas_info = {
        "Atlas_name": "pbmc_6k_v1_v1",
        "Atlas_type": "Cellranger < v3",
        "Atlas_extension": ".mtx",
        "Atlas_path": "/data/analysis/data_becavin/"
        "checkatlas_test/tuto/data5/"
        "pbmc_6k_v1_v1/outs/"
        "filtered_matrices_mex/hg19/matrix.mtx",
    }

    """ atlas_info = {'Atlas_name': 'pbmc_5k_v3_v3',
                   'Atlas_type': 'Cellranger >= v3',
                     'Atlas_extension': '.h5',
                       'Atlas_path':
                       '/data/analysis/data_becavin/'
                       'checkatlas_test/tuto/data5/'
                       'pbmc_5k_v3_v3/outs/'
                       '5k_pbmc_v3_filtered_feature_bc_matrix.h5'}

    atlas_info = {'Atlas_name': 'pbmc_5k_v3_v7',
                   'Atlas_type': 'Cellranger >= v3',
                     'Atlas_extension': '.h5',
                       'Atlas_path':
                       '/data/analysis/data_becavin/'
                       'checkatlas_test/tuto/data5/'
                       'pbmc_5k_v3_v7/outs/'
                       'SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5'}

    atlas_info = {'Atlas_name': 'pbmc_3k_multiome',
                   'Atlas_type': 'Cellranger >= v3',
                     'Atlas_extension': '.h5',
                       'Atlas_path':
                       '/data/analysis/data_becavin/'
                       'checkatlas_test/tuto/data5/'
                       'pbmc_3k_multiome/outs/'
                       'pbmc_unsorted_3k_filtered_feature_bc_matrix.h5'} """

    adata = atlas.read_atlas(atlas_info)
    # atlas_list = list_atlases(path)
