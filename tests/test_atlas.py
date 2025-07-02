import argparse
import os

import pytest
from anndata import AnnData

from checkatlas import atlas, cellranger, checkatlas
from checkatlas.utils import files, folders

from . import datasets

given = pytest.mark.parametrize


@given("atlas_info,expected", [(datasets.get_scanpy_atlas_info(), AnnData)])
def test_read_scanpy_atlas(atlas_info, expected):
    print("Run test")
    # atlas_seurat.check_seurat_install()
    print(atlas_info)
    adata = atlas.read_atlas(atlas_info)
    assert isinstance(adata, expected)


@given(
    "atlas_info,expected", [(datasets.get_cellranger_atlas_info(), AnnData)]
)
def test_read_cellranger_atlas(atlas_info, expected):
    adata = cellranger.read_cellranger_current(atlas_info)
    assert isinstance(adata, expected)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_clean_scanpy(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    adata_clean = atlas.clean_scanpy_atlas(adata, atlas_info)
    list_var = adata_clean.var_names
    assert len(set(list_var)) == len(list_var)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_viable_obs_qc(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    expected = ["n_genes_by_counts", "total_counts"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc_display")
    args = parser.parse_args(
        [
            "--qc_display",
            [
                "violin_plot",
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
            ],
        ]
    )
    obs_keys = atlas.get_viable_obs_qc(adata, args)
    assert obs_keys == expected


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_viable_obs_annot(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    expected = ["leiden", "louvain"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(["--obs_cluster", atlas.OBS_CLUSTERS])
    obs_keys = atlas.get_viable_obs_annot(adata, args)
    assert obs_keys == expected


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_viable_obsm(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    expected = ["X_draw_graph_fr", "X_pca", "X_tsne", "X_umap"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(["--obs_cluster", atlas.OBS_CLUSTERS])
    obs_keys = atlas.get_viable_obsm(adata, args)
    assert obs_keys == expected


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_summary_table(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_summary_table(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.SUMMARY, checkatlas.TSV_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_adata_table(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_anndata_table(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.ANNDATA, checkatlas.TSV_EXTENSION, args.path
    )
    print(csv_path)
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_qc_table(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--qc_display")
    parser.add_argument("--plot_celllimit", type=int)
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        [
            "--path",
            checkatlas_path,
            "--qc_display",
            [
                "violin_plot",
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
            ],
            "--plot_celllimit",
            "10000",
        ]
    )
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_qc_tables(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.QC, checkatlas.TSV_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_qc_plots(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_qc_plots(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.QC_FIG, checkatlas.QC_FIG_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_umap_plots(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        ["--path", checkatlas_path, "--obs_cluster", atlas.OBS_CLUSTERS]
    )

    folders.checkatlas_folders(checkatlas_path)
    atlas.create_umap_fig(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.UMAP, checkatlas.UMAP_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_tsne_plots(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        ["--path", checkatlas_path, "--obs_cluster", atlas.OBS_CLUSTERS]
    )

    folders.checkatlas_folders(checkatlas_path)
    atlas.create_tsne_fig(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name, folders.TSNE, checkatlas.TSNE_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_cluster_metric(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    parser.add_argument("--metric_cluster")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        [
            "--path",
            checkatlas_path,
            "--metric_cluster",
            ["davies_bouldin"],
            "--obs_cluster",
            atlas.OBS_CLUSTERS,
        ]
    )
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_metric_cluster(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name,
        folders.CLUSTER,
        checkatlas.TSV_EXTENSION,
        args.path,
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_annot_metric(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    parser.add_argument("--metric_annot")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        [
            "--path",
            checkatlas_path,
            "--metric_annot",
            ["rand_index"],
            "--obs_cluster",
            atlas.OBS_CLUSTERS,
        ]
    )
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_metric_annot(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name,
        folders.ANNOTATION,
        checkatlas.TSV_EXTENSION,
        args.path,
    )
    assert os.path.exists(csv_path)


@given("atlas_info", [datasets.get_scanpy_atlas_info()])
def test_dimred_metric(atlas_info):
    adata = atlas.read_atlas(atlas_info)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    parser.add_argument("--metric_dimred")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        [
            "--path",
            checkatlas_path,
            "--metric_dimred",
            ["kruskal_stress"],
            "--obs_cluster",
            atlas.OBS_CLUSTERS,
        ]
    )
    folders.checkatlas_folders(checkatlas_path)
    atlas.create_metric_dimred(adata, atlas_info, args)
    atlas_name = atlas_info[checkatlas.ATLAS_NAME_KEY]
    csv_path = files.get_file_path(
        atlas_name,
        folders.DIMRED,
        checkatlas.TSV_EXTENSION,
        args.path,
    )
    assert os.path.exists(csv_path)
