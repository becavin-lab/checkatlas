import os

import scanpy as sc

from checkatlas import atlas, cellranger, seurat

"""
Module in which we create the data for all tests
Anndata, CellRanger, Seurat, MuonData, SpatialData
This module should be ran as a stand-alone :
python create-test_data.py
"""

DATA_PATH = os.path.join(os.getcwd(), "tests", "data")
ADATA_TEST_PATH = os.path.join(DATA_PATH, "pbmc3k_scanpy.h5ad")
CELLRANGER_FOLDER_PATH = os.path.join(DATA_PATH, "pbmc3k_cellranger")
CELLRANGER_TEST_PATH = os.path.join(
    CELLRANGER_FOLDER_PATH, "outs", cellranger.CELLRANGER_FILE
)
SEURAT_TEST_PATH = os.path.join(DATA_PATH, "pbmc3k_seurat.rds")


def get_scanpy_atlas_info():
    atlas_info = {
        "Atlas_name": "pbmc3k_scanpy",
        "Atlas_type": atlas.ANNDATA_TYPE,
        "Atlas_extension": atlas.ANNDATA_EXTENSION,
        "Atlas_path": ADATA_TEST_PATH,
    }
    return atlas_info


def get_cellranger_atlas_info():
    atlas_info = {
        "Atlas_name": "pbmc3k_cellranger",
        "Atlas_type": cellranger.CELLRANGER_TYPE_CURRENT,
        "Atlas_extension": ".h5",
        "Atlas_path": CELLRANGER_TEST_PATH,
    }
    return atlas_info


def get_seurat_atlas_info():
    atlas_info = {
        "Atlas_name": "pbmc3k_seurat",
        "Atlas_type": seurat.SEURAT_TYPE,
        "Atlas_extension": seurat.SEURAT_EXTENSION,
        "Atlas_path": SEURAT_TEST_PATH,
    }
    return atlas_info


def create_adata_data():
    """Create adata for all tests using
    DL from cellxgene : pbmc3k.h5ad
    TSNE, and leiden are calculated.
    Than sampling to 100 cells and 1000 genes
    adata is saved in ADATA_TEST_PATH
    """
    url_data = (
        "https://github.com/chanzuckerberg/"
        "cellxgene/raw/main/example-dataset/pbmc3k.h5ad"
    )
    download_cmd = f"curl --location -o {ADATA_TEST_PATH} {url_data}"
    print(f"Run download : {download_cmd}")
    os.system(download_cmd)
    adata = adata = sc.read_h5ad(ADATA_TEST_PATH)
    # Create qC values
    qc_genes = []
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("mt")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("ribo")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_genes,
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # run leiden clustering
    sc.tl.leiden(adata)
    # sample to 100 cells and 1000 genes
    sc.pp.subsample(adata, n_obs=100)
    adata.var[1:1000]
    adata = adata[:, adata.var[1:1000].index]
    adata.write_h5ad(ADATA_TEST_PATH)
    print("Adata for tests saved in : " f"{ADATA_TEST_PATH}")


def create_cellranger_data():
    """
    Create cellranger file for all tests functions
    Download from 10xgenomics website : 5k_pbmc_v3
    """
    os.mkdir(CELLRANGER_FOLDER_PATH)
    os.mkdir(os.path.join(CELLRANGER_FOLDER_PATH, "outs"))
    url_data = (
        "https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/"
        "5k_pbmc_v3_filtered_feature_bc_matrix.h5"
    )
    download_cmd = "curl --location -o " + f"{CELLRANGER_TEST_PATH} {url_data}"
    print(f"Run download : {download_cmd}")
    os.system(download_cmd)
    print(f"Cellranger data for tests saved in : {CELLRANGER_TEST_PATH}")


def create_seurat_data():
    """
    Create Seurat data for all tests function

    """

    url_data = (
        "https://www.dropbox.com/s/63gnlw45jf7cje8/" "pbmc3k_final.rds?dl=1"
    )
    download_cmd = "curl --location -o " + f"{SEURAT_TEST_PATH} {url_data}"
    print(f"Run download : {download_cmd}")
    os.system(download_cmd)
    # R need to be installed to run that line
    r_command = "Rscript tests/data/create_test_data.R"
    print(f"Run {r_command} inside python or directly from R")
    # os.system(r_command)
    print(f"Seurat data for tests saved in : {SEURAT_TEST_PATH}")


def create_muon_data():
    print("TO DO")


def create_spatial_data():
    print("TO DO")


if __name__ == "__main__":
    create_adata_data()
    # create_cellranger_data()
    # create_seurat_data()
    create_muon_data()
    create_spatial_data()
