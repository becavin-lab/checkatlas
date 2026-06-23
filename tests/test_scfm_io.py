"""Tests for the scfm io / cache layer.

These tests cover the wide-to-long conversion that bridges the
on-disk per-task TSVs (``checkatlas_files/<task>/<atlas>.tsv``) and
the long-format metric table consumed by the diagnostic engine.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from checkatlas.scfm import io
from checkatlas.utils import folders


def _write_tsv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _cluster_wide() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Clust_Sample": "blood_X_pca_scvi_leiden_donorassay_full",
                "obs": "scvi_leiden_donorassay_full",
                "silhouette": 0.21,
                "silhouette_running_time": 0.123,
                "calinski_harabasz": 10983.16,
                "calinski_harabasz_running_time": 0.081,
            },
            {
                "Clust_Sample": "blood_X_scvi_scvi_leiden_donorassay_full",
                "obs": "scvi_leiden_donorassay_full",
                "silhouette": 0.25,
                "silhouette_running_time": 0.456,
                "calinski_harabasz": 12000.0,
                "calinski_harabasz_running_time": 0.05,
            },
        ]
    )


def _annot_wide() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Annot_Sample": "blood_X_pca__scvi_batch",
                "Reference": "X_pca",
                "obs": "_scvi_batch",
                "iLISI": 0.017,
                "iLISI_running_time": 0.103,
                "kbet": 0.996,
                "kbet_running_time": 0.059,
            },
            {
                "Annot_Sample": "blood_X_pca_assay",
                "Reference": "X_pca",
                "obs": "assay",
                "iLISI": 0.005,
                "iLISI_running_time": 0.099,
                "kbet": 0.116,
                "kbet_running_time": 0.059,
            },
            {
                "Annot_Sample": "blood_X_scvi_assay",
                "Reference": "X_scvi",
                "obs": "assay",
                "iLISI": 0.50,
                "iLISI_running_time": 0.123,
                "kbet": 0.30,
                "kbet_running_time": 0.05,
            },
        ]
    )


def _dimred_wide() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Dimred_Sample": "blood_X_pca",
                "obsm": "X_pca",
                "dCor": 0.0,
                "dCor_running_time": 0.482,
                "kruskal_stress": 0.0,
                "kruskal_stress_running_time": 23.189,
            },
            {
                "Dimred_Sample": "blood_X_scvi",
                "obsm": "X_scvi",
                "dCor": 0.0,
                "dCor_running_time": 0.454,
                "kruskal_stress": 0.0,
                "kruskal_stress_running_time": 5.549,
            },
        ]
    )


def test_wide_cluster_to_long():
    df = _cluster_wide()
    long = io.wide_to_long(df, "cluster", "blood")
    assert list(long["Atlas Name"].unique()) == ["blood"]
    assert set(long.columns) >= {
        "Atlas Name",
        "Task",
        "Metric Name",
        "Embedding",
        "Reference/Input 1",
        "Prediction/Input 2",
        "Value",
        "Time (s)",
    }
    assert set(long["Task"].unique()) == {"cluster"}
    assert set(long["Embedding"].unique()) == {"X_pca", "X_scvi"}
    assert set(long["Metric Name"].unique()) == {
        "silhouette",
        "calinski_harabasz",
    }
    assert long["Value"].notna().all()


def test_wide_annot_to_long():
    df = _annot_wide()
    long = io.wide_to_long(df, "annotation", "blood")
    assert set(long["Task"].unique()) == {"annot"}
    assert set(long["Embedding"].unique()) == {"X_pca", "X_scvi"}
    # obs = "_scvi_batch" rows should not be NaN
    scvi_batch_rows = long[
        long["Prediction/Input 2"] == "_scvi_batch"
    ]
    assert len(scvi_batch_rows) >= 1
    assert scvi_batch_rows["Value"].notna().all()


def test_wide_dimred_to_long():
    df = _dimred_wide()
    long = io.wide_to_long(df, "dimred", "blood")
    assert set(long["Task"].unique()) == {"dimred"}
    assert set(long["Embedding"].unique()) == {"X_pca", "X_scvi"}
    assert set(long["Reference/Input 1"].unique()) == {"X"}
    assert long["Prediction/Input 2"].eq("").all()


def test_wide_to_long_skips_nan_values():
    df = _cluster_wide()
    # Add a row where silhouette is NaN
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                [
                    {
                        "Clust_Sample": "blood_X_pca_label_empty",
                        "obs": "label_empty",
                        "silhouette": np.nan,
                        "silhouette_running_time": 0.0,
                        "calinski_harabasz": 100.0,
                        "calinski_harabasz_running_time": 0.05,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    long = io.wide_to_long(df, "cluster", "blood")
    # The NaN silhouette row should not produce a long row for silhouette,
    # but should still produce one for calinski_harabasz.
    label_empty_rows = long[long["Prediction/Input 2"] == "label_empty"]
    assert (label_empty_rows["Metric Name"] == "calinski_harabasz").all()
    assert "silhouette" not in set(label_empty_rows["Metric Name"])


def test_wide_to_long_empty_input():
    out = io.wide_to_long(pd.DataFrame(), "cluster", "blood")
    assert out.empty


def test_load_per_task_tsvs_missing(tmp_path: Path):
    tsvs = io.load_per_task_tsvs(str(tmp_path), "nonexistent")
    assert tsvs == {"cluster": None, "annotation": None, "dimred": None}


def test_load_per_task_tsvs_present(tmp_path: Path):
    base = str(tmp_path)
    _write_tsv(_cluster_wide(), folders.get_folder(base, "cluster") + "/blood.tsv")
    _write_tsv(_annot_wide(), folders.get_folder(base, "annotation") + "/blood.tsv")
    _write_tsv(_dimred_wide(), folders.get_folder(base, "dimred") + "/blood.tsv")
    tsvs = io.load_per_task_tsvs(base, "blood")
    assert tsvs["cluster"] is not None
    assert tsvs["annotation"] is not None
    assert tsvs["dimred"] is not None
    assert len(tsvs["cluster"]) == 2
    assert len(tsvs["annotation"]) == 3
    assert len(tsvs["dimred"]) == 2


def test_load_per_task_tsvs_partial(tmp_path: Path):
    base = str(tmp_path)
    _write_tsv(_cluster_wide(), folders.get_folder(base, "cluster") + "/blood.tsv")
    tsvs = io.load_per_task_tsvs(base, "blood")
    assert tsvs["cluster"] is not None
    assert tsvs["annotation"] is None
    assert tsvs["dimred"] is None


def test_run_scfm_from_cache_with_only_scfm_metrics(tmp_path: Path):
    """End-to-end: with no per-task TSVs on disk, ``run_scfm_from_cache``
    must still produce the 7 standard output files and the 9 verdicts.
    """
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.pipeline import run_scfm_from_cache

    import anndata as ad

    n = 60
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=rng.standard_normal((n, 30)),
        obs=pd.DataFrame(
            {
                "celltype": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma"], n)
                ),
                "celltype_scfm_pred": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma"], n)
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2", "b3"], n)
                ),
                "tech": pd.Categorical(
                    rng.choice(["t1", "t2", "t3", "t4"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, 20)).astype(np.float32)
    adata.obsm["X_scvi_dropin"] = adata.obsm["X_pca"].copy()
    adata.obsm["X_geneformer"] = (
        adata.obsm["X_pca"] + 0.5 * rng.standard_normal((n, 20))
    ).astype(np.float32)

    cfg = SCFMConfig(
        atlas_name="blood",
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca", "X_scvi_dropin"),
        ref_label="celltype",
        predicted_label="celltype_scfm_pred",
        batch_key="batch",
        domain_key="tech",
        n_seeds=2,
        min_domain_size=5,
    )

    import argparse

    args = argparse.Namespace(path=str(tmp_path), n_jobs=2)

    result = run_scfm_from_cache(adata, cfg, args=args)

    # All 7 output files exist (no scfm_ prefix, atlas name is in
    # the parent path).
    expected = {
        "verdicts.tsv",
        "composite.tsv",
        "per_metric.tsv",
        "inputs.tsv",
        "grade_legend.md",
        "resolved_weights.json",
        "resolved_thresholds.yaml",
    }
    assert expected.issubset(set(result["out_paths"].keys())), (
        f"missing output files: {expected - set(result['out_paths'].keys())}"
    )
    # New schema: 9 headline verdicts + 9 per-combo verdicts = 18.
    # Every verdict covers a problem_id in 1..9.
    assert len(result["verdicts"]) >= 9
    for v in result["verdicts"]:
        assert v.problem_id in range(1, 10)
    # The first 9 verdicts are the headline (worst grade per problem)
    assert result["verdicts"][0].combo_id == "headline"
    # Output dir is <tmp>/checkatlas_files/scfm/blood
    assert result["outdir"].endswith(
        os.path.join("checkatlas_files", "scfm", "blood")
    )


def test_run_scfm_from_cache_all_combinations_mode(tmp_path: Path):
    """End-to-end: with no --scfm_* flags set, ``run_scfm_from_cache``
    falls into the 'all combinations' mode: it auto-runs the 3
    per-task engines and produces 9 verdicts × N combos plus a
    headline row. The output TSVs get a ``combo_id`` column and a
    ``remark`` column.
    """
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.pipeline import run_scfm_from_cache

    import anndata as ad

    n = 80
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=rng.standard_normal((n, 30)),
        obs=pd.DataFrame(
            {
                "celltype": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma", "delta"], n)
                ),
                "celltype_scfm_pred": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma", "delta"], n)
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2", "b3"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, 20)).astype(np.float32)
    adata.obsm["X_geneformer"] = (
        adata.obsm["X_pca"] + 0.5 * rng.standard_normal((n, 20))
    ).astype(np.float32)

    # No --scfm_* flags set; orchestrator runs all 3 engines
    # and the pipeline enumerates all detected (ref, pred, batch,
    # scfm, baseline) combinations.
    cfg = SCFMConfig(
        atlas_name="x",
        scaling_fractions=(),  # skip cal_scfm to keep the test fast
        n_seeds=0,
    )

    import argparse

    args = argparse.Namespace(
        path=str(tmp_path),
        n_jobs=2,
        metric_cluster=["silhouette"],
        metric_annot=["adj_rand_index"],
        metric_dimred=["dCor"],
        scfm_fast=False,
        scfm_max_combos=10,
    )

    result = run_scfm_from_cache(adata, cfg, args=args)

    # All 7 output files exist (back-compat names)
    expected = {
        "verdicts.tsv",
        "composite.tsv",
        "per_metric.tsv",
        "inputs.tsv",
        "grade_legend.md",
        "resolved_weights.json",
        "resolved_thresholds.yaml",
    }
    assert expected.issubset(set(result["out_paths"].keys())), (
        f"missing output files: {expected - set(result['out_paths'].keys())}"
    )

    # The headline row is at the top of verdicts (combo_id == 'headline')
    verdicts = result["verdicts"]
    assert verdicts[0].combo_id == "headline"
    assert len(verdicts) >= 1 + 9  # headline (9) + at least one combo (9)

    # All combo verdicts have a non-empty combo_id and a remark
    for v in verdicts:
        if v.combo_id == "headline":
            continue
        assert v.combo_id
        assert v.combo is not None

    # The composite has a headline row + at least one per-combo row
    composite_df = pd.read_csv(result["out_paths"]["composite.tsv"], sep="\t")
    assert "headline" in composite_df["combo_id"].values
    n_per_combo = (composite_df["combo_id"] != "headline").sum()
    assert n_per_combo >= 1

    # The per-metric TSV has a combo_id column
    pm_df = pd.read_csv(result["out_paths"]["per_metric.tsv"], sep="\t")
    assert "combo_id" in pm_df.columns
    assert "remark" in pm_df.columns

    # The verdicts TSV has a remark column
    vd_df = pd.read_csv(result["out_paths"]["verdicts.tsv"], sep="\t")
    assert "remark" in vd_df.columns
    assert "combo_id" in vd_df.columns

    # The headline_scores are in inputs.tsv
    inputs_df = pd.read_csv(result["out_paths"]["inputs.tsv"], sep="\t")
    assert "headline_fmf" in inputs_df.columns
    assert "headline_bf" in inputs_df.columns
    assert "headline_pr" in inputs_df.columns


def test_run_scfm_from_cache_scfm_fast_mode_single_combo(
    tmp_path: Path,
):
    """End-to-end: with --scfm_fast and no --scfm_* flags, the
    pipeline should produce exactly one combo (back-compat with
    today's single-combo output).
    """
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.pipeline import run_scfm_from_cache

    import anndata as ad

    n = 60
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=rng.standard_normal((n, 30)),
        obs=pd.DataFrame(
            {
                "celltype": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma"], n)
                ),
                "celltype_scfm_pred": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma"], n)
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2", "b3"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, 20)).astype(np.float32)
    adata.obsm["X_geneformer"] = (
        adata.obsm["X_pca"] + 0.5 * rng.standard_normal((n, 20))
    ).astype(np.float32)

    cfg = SCFMConfig(
        atlas_name="x",
        scaling_fractions=(),  # skip cal_scfm to keep the test fast
        n_seeds=0,
    )

    import argparse

    args = argparse.Namespace(
        path=str(tmp_path),
        n_jobs=2,
        metric_cluster=["silhouette"],
        metric_annot=["adj_rand_index"],
        metric_dimred=["dCor"],
        scfm_fast=True,  # opt out of all-combinations sweep
        scfm_max_combos=10,
    )

    result = run_scfm_from_cache(adata, cfg, args=args)

    # 7 output files still produced
    expected = {
        "verdicts.tsv",
        "composite.tsv",
        "per_metric.tsv",
        "inputs.tsv",
        "grade_legend.md",
        "resolved_weights.json",
        "resolved_thresholds.yaml",
    }
    assert expected.issubset(set(result["out_paths"].keys()))

    # Headline + 1 combo = 1 × 9 + 9 = 18 verdict rows
    per_combo_count = sum(
        1 for v in result["verdicts"] if v.combo_id != "headline"
    ) // 9
    assert per_combo_count == 1, (
        f"expected exactly 1 combo in --scfm_fast mode, got {per_combo_count}"
    )
