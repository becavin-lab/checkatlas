import pytest
import scanpy as sc

from checkatlas.atlas import clean_scanpy_atlas

given = pytest.mark.parametrize


def test_atlas_object():
    adata = sc.datasets.pbmc68k_reduced()
    atlas_info = ["PBMC68k", "Scanpy", ".h5ad", "Scanpy module"]
    assert clean_scanpy_atlas(adata, atlas_info)


# test_atlas_object()
# @given("fn", [atlas(), list_atlases()])
# def test_parameterized(fn):
#     assert "hello from" in fn()
#
#
# def test_list_atlases():
#     assert list_atlases(".") == "hello from base function"
#
#
# def test_atlas():
#     assert atlas() == "hello from BaseClass"
