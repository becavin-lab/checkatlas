import argparse
import os

import pytest
from rpy2.robjects.methods import RS4

from checkatlas import atlas, atlas_seurat, checkatlas
from checkatlas.utils import files, folders

from . import datasets

given = pytest.mark.parametrize


@given("atlas_path,expected", [(datasets.SEURAT_TEST_PATH, RS4)])
def test_read_seurat_atlas(atlas_path, expected):
    atlas_seurat.check_seurat_install()
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    assert type(seurat_data) == expected


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_viable_obs_qc(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    expected = ["nCount_RNA", "nFeature_RNA"]
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
    obs_keys = atlas_seurat.get_viable_obs_qc(seurat_data, args)
    assert obs_keys == expected


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_viable_obs_annot(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    expected = ["RNA_snn_res.0.5", "seurat_clusters"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(["--obs_cluster", atlas.OBS_CLUSTERS])
    obs_keys = atlas_seurat.get_viable_obs_annot(seurat_data, args)
    assert obs_keys == expected


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_viable_obsm(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    expected = ["pca", "umap", "tsne"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_cluster")
    args = parser.parse_args(["--obs_cluster", atlas.OBS_CLUSTERS])
    obs_keys = atlas_seurat.get_viable_obsm(seurat_data, args)
    print(obs_keys)
    assert obs_keys == expected


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_summary_table(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.create_summary_table(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.SUMMARY, checkatlas.SUMMARY_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_seurat_data_table(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.create_anndata_table(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.ANNDATA, checkatlas.ADATA_EXTENSION, args.path
    )
    print(csv_path)
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_qc_table(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
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
    atlas_seurat.create_qc_tables(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.QC, checkatlas.QC_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_qc_plots(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(["--path", checkatlas_path])
    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.create_qc_plots(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.QC_FIG, checkatlas.QC_FIG_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_umap_plots(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        ["--path", checkatlas_path, "--obs_cluster", atlas.OBS_CLUSTERS]
    )

    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.create_umap_fig(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.UMAP, checkatlas.UMAP_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_tsne_plots(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(
        ["--path", checkatlas_path, "--obs_cluster", atlas.OBS_CLUSTERS]
    )

    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.create_tsne_fig(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name, folders.TSNE, checkatlas.TSNE_EXTENSION, args.path
    )
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_cluster_metric(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
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
    atlas_seurat.metric_cluster(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(
        atlas_name,
        folders.CLUSTER,
        checkatlas.METRIC_CLUSTER_EXTENSION,
        args.path,
    )
    assert os.path.exists(csv_path)


""" @given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_annot_metric(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    parser.add_argument("--metric_annot")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(['--path',checkatlas_path,
                              '--metric_annot',
                                ["rand_index"],
                                "--obs_cluster",atlas.OBS_CLUSTERS])
    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.metric_annot(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(atlas_name, folders.ANNOTATION,
                                   checkatlas.METRIC_ANNOTATION_EXTENSION,
                                   args.path)
    assert os.path.exists(csv_path)


@given("atlas_path", [(datasets.SEURAT_TEST_PATH)])
def test_dimred_metric(atlas_path):
    seurat_data = atlas_seurat.read_atlas(atlas_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--obs_cluster")
    parser.add_argument("--metric_dimred")
    checkatlas_path = os.getcwd()
    args = parser.parse_args(['--path',checkatlas_path,
                              '--metric_dimred',
                                ["kruskal_stress"],
                                "--obs_cluster",atlas.OBS_CLUSTERS])
    folders.checkatlas_folders(checkatlas_path)
    atlas_seurat.metric_dimred(seurat_data, atlas_path, args)
    atlas_name = checkatlas.get_atlas_name(atlas_path)
    csv_path = files.get_file_path(atlas_name, folders.DIMRED,
                                   checkatlas.METRIC_DIMRED_EXTENSION,
                                   args.path)
    assert os.path.exists(csv_path) """
