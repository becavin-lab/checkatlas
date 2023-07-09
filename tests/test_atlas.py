import pytest
from anndata import AnnData

import checkatlas.atlas as atlas

from .data import datasets

given = pytest.mark.parametrize


@given("path_input,expected", [(datasets.ADATA_TEST_PATH, AnnData)])
def test_read_scanpy_atlas(path_input, expected):
    adata = atlas.read_atlas(path_input)
    print(type(adata))
    print(type(AnnData))

    assert type(adata) == expected


@given("path_input,expected", [(datasets.CELLRANGER_TEST_PATH, AnnData)])
def test_read_cellranger_atlas(path_input, expected):
    adata = atlas.read_cellranger(path_input)
    print(type(adata))
    print(type(AnnData))

    assert type(adata) == expected
