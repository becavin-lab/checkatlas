from enum import unique
from traceback import print_list
import pytest
from anndata import AnnData
import argparse
import checkatlas
import checkatlas.atlas as atlas

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
                                ['X_draw_graph_fr', 'X_pca', 'X_tsne', 'X_umap'])])
def test_viable_obsm(path_input, expected):
    adata = atlas.read_atlas(path_input)
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(['--obs_cluster',atlas.OBS_CLUSTERS])
    obs_keys = atlas.get_viable_obsm(adata, args)
    print(obs_keys)
    assert obs_keys == expected