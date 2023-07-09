import os
import pytest
from anndata import AnnData
import argparse
from checkatlas import atlas
from checkatlas import folders
from checkatlas import checkatlas

from .data import datasets

given = pytest.mark.parametrize


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, AnnData)])
def test_read_scanpy_atlas(path_input, expected):
    adata = atlas.read_atlas(path_input)
    assert type(adata) == expected


@given("path_input,expected", [(datasets.CELLRANGER_TEST_PATH, AnnData)])
def test_read_cellranger_atlas(path_input, expected):
    adata = atlas.read_cellranger(path_input)
    assert type(adata) == expected


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, AnnData)])
def test_clean_scanpy(path_input, expected):
    """Test if var-names is unique after cleaning

    Args:
        path_input (_type_): _description_
        expected (_type_): _description_
    """    
    adata = atlas.read_atlas(path_input)
    adata_clean = atlas.clean_scanpy_atlas(adata, path_input)
    list_var = adata_clean.var_names
    assert len(set(list_var)) == len(list_var)


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, 
                                ['n_genes_by_counts', 'total_counts'])])
def test_viable_obs_qc(path_input, expected):
    adata = atlas.read_atlas(path_input)
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc_display")
    args = parser.parse_args(['--qc_display',[
            "violin_plot",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ]])
    obs_keys = atlas.get_viable_obs_qc(adata, args)
    assert obs_keys == expected


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, 
                                ['leiden', 'louvain'])])
def test_viable_obs_annot(path_input, expected):
    adata = atlas.read_atlas(path_input)
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(['--obs_cluster',atlas.OBS_CLUSTERS])
    obs_keys = atlas.get_viable_obs_annot(adata, args)
    print(obs_keys)
    assert obs_keys == expected


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, 
                                True)])
def test_summary_table(path_input, expected):
    adata = atlas.read_atlas(path_input)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(['--path',checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_summary_table(adata, checkatlas_path, args)
    atlas_name = checkatlas.get_atlas_name(path_input)
    csv_path = os.path.join(
        folders.get_folder(args.path, folders.SUMMARY),
        atlas_name + checkatlas.SUMMARY_EXTENSION,
    )
    print("exiists",os.path.exists(csv_path))
    print(csv_path)
    assert os.path.exists(csv_path)