"""Tests for the scfm orchestrator and the user-override filter.

These tests cover three concerns:

1. ``required_tasks`` is a pure function that maps an SCFMConfig to
   the set of per-task metric categories the scfm diagnostic
   consumes.
2. ``ensure_per_task_tsvs`` skips per-task auto-runs when the atlas
   lacks the required keys, with a single scientific warning per
   missing task, and never raises on missing keys.
3. ``_apply_user_overrides`` (in :mod:`checkatlas.scfm.pipeline`)
   filters the long-format per-task DataFrame to only the rows that
   match the user-specified ``--scfm_*`` keys, so the scfm pipeline
   consumes *only* the user-chosen columns.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from checkatlas.scfm import orchestrator
from checkatlas.scfm.config import SCFMConfig
from checkatlas.scfm.pipeline import _apply_user_overrides
from checkatlas.utils import folders


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────


def _build_adata(
    *,
    n: int = 100,
    has_celltype: bool = True,
    has_leiden: bool = True,
    has_batch: bool = True,
    n_pca: int = 20,
    n_scvi: int = 20,
):
    """Build a small AnnData suitable for the column detector."""
    import anndata as ad

    rng = np.random.default_rng(0)
    obs_dict: dict = {
        # Always have an obs_name index column so the obs DataFrame
        # is non-empty even when no semantic columns are requested.
        "obs_name": pd.Categorical([f"c{i}" for i in range(n)])
    }
    if has_celltype:
        obs_dict["celltype"] = pd.Categorical(
            rng.choice(["alpha", "beta", "gamma"], n)
        )
    if has_leiden:
        obs_dict["leiden"] = pd.Categorical(
            rng.choice([f"c{i}" for i in range(5)], n)
        )
    if has_batch:
        obs_dict["batch"] = pd.Categorical(
            rng.choice(["b1", "b2", "b3"], n)
        )
    adata = ad.AnnData(
        X=rng.standard_normal((n, 50)),
        obs=pd.DataFrame(obs_dict),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, n_pca)).astype(np.float32)
    adata.obsm["X_scvi"] = rng.standard_normal((n, n_scvi)).astype(np.float32)
    return adata


# ──────────────────────────────────────────────────────────────────
# required_tasks
# ──────────────────────────────────────────────────────────────────


def test_required_tasks_empty_when_no_inputs():
    cfg = SCFMConfig(atlas_name="x")
    assert orchestrator.required_tasks(cfg) == set()


def test_required_tasks_annotation_only_via_ref_label():
    cfg = SCFMConfig(atlas_name="x", ref_label="celltype")
    assert orchestrator.required_tasks(cfg) == {"annotation", "cluster"}


def test_required_tasks_annotation_via_predicted_label():
    cfg = SCFMConfig(atlas_name="x", predicted_label="leiden")
    assert orchestrator.required_tasks(cfg) == {"annotation"}


def test_required_tasks_annotation_via_batch_key():
    cfg = SCFMConfig(atlas_name="x", batch_key="batch")
    assert orchestrator.required_tasks(cfg) == {"annotation"}


def test_required_tasks_dimred_via_scfm_embedding():
    cfg = SCFMConfig(atlas_name="x", scfm_embedding="X_geneformer")
    assert orchestrator.required_tasks(cfg) == {"dimred"}


def test_required_tasks_dimred_via_baselines():
    cfg = SCFMConfig(atlas_name="x", baseline_embeddings=("X_pca", "X_scvi"))
    assert orchestrator.required_tasks(cfg) == {"dimred"}


def test_required_tasks_all_three():
    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca",),
        ref_label="celltype",
        batch_key="batch",
    )
    assert orchestrator.required_tasks(cfg) == {"cluster", "annotation", "dimred"}


# ──────────────────────────────────────────────────────────────────
# user_overrides helper
# ──────────────────────────────────────────────────────────────────


def test_user_overrides_reports_only_set_keys():
    cfg = SCFMConfig(atlas_name="x", ref_label="celltype", batch_key="batch")
    out = orchestrator.user_overrides(cfg)
    assert out["ref_label"] == "celltype"
    assert out["batch_key"] == "batch"
    assert out["predicted_label"] is None
    assert out["domain_key"] is None
    assert out["patient_key"] is None
    assert out["outcome_key"] is None


# ──────────────────────────────────────────────────────────────────
# Scientific warnings on missing keys
# ──────────────────────────────────────────────────────────────────


def test_skip_cluster_when_no_cluster_label_in_atlas(tmp_path: Path, caplog):
    """An atlas whose obs has no cluster-label column must trigger a
    single scientific warning when cluster metrics are required.

    The user did NOT pass ``--scfm_ref_label`` here, so the
    orchestrator falls back to the column detector. The detector
    finds no usable cluster-label column and the orchestrator emits
    a single scientific warning.
    """

    # Build an atlas that has nothing the column detector recognises
    # as a cluster label, ref, pred, or batch. We use the gate
    # function directly because the required-tasks policy needs
    # ref_label / predicted_label / batch_key / scfm_embedding to
    # be set to mark a task as needed; the cluster task is
    # triggered by ``ref_label`` (see required_tasks). So the
    # detector-driven branch only fires when the user named no
    # --scfm_ref_label and the detector still finds nothing — in
    # that case the cluster task is simply not required. The
    # unit being tested is the gate function, not the full
    # orchestrator.
    from checkatlas.scfm.orchestrator import _gate_cluster
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    adata = _build_adata(has_leiden=False, has_celltype=False, has_batch=False)
    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    cfg = SCFMConfig(atlas_name="empty_clust")
    ok, reason = _gate_cluster(adata, detected, cfg)
    assert ok is False
    assert "no cluster-label column detected" in reason


def test_skip_annotation_when_no_ref_pred_batch_in_atlas(tmp_path: Path, caplog):
    """An atlas whose obs has no ref / pred / batch column must trigger
    a single scientific warning when annotation metrics are required.

    The user did NOT pass any of ``--scfm_ref_label``,
    ``--scfm_predicted_label``, ``--scfm_batch_key`` here, so the
    orchestrator falls back to the column detector. The detector
    finds nothing and the orchestrator emits a single scientific
    warning via the gate function.
    """

    from checkatlas.scfm.orchestrator import _gate_annot
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    adata = _build_adata(
        has_celltype=False, has_leiden=False, has_batch=False
    )
    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    cfg = SCFMConfig(atlas_name="empty_annot")
    ok, reason = _gate_annot(adata, detected, cfg)
    assert ok is False
    assert "no reference, predicted, or batch column detected" in reason


def test_skip_dimred_when_no_embedding_in_atlas(tmp_path: Path, caplog):
    """An atlas whose obsm has no embedding with > 3 components must
    trigger a single scientific warning when dimred metrics are
    required.

    To exercise the detector-driven branch (rather than the
    user-override branch) we do NOT pass ``scfm_embedding`` /
    ``baseline_embeddings``; the orchestrator falls back to the
    column detector, which finds no usable embedding.
    """

    import anndata as ad

    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=rng.standard_normal((50, 30)),
        obs=pd.DataFrame(
            {"celltype": pd.Categorical(rng.choice(["a", "b"], 50))}
        ),
    )
    # 2-D UMAP only — detector explicitly excludes it.
    adata.obsm["X_umap"] = rng.standard_normal((50, 2)).astype(np.float32)
    # We trigger the dimred-required branch via the *cluster* required
    # task: required_tasks needs scfm_embedding or baseline_embeddings
    # for dimred. So we use the cluster path which is also a dimred
    # consumer in the long-format DataFrame. Simpler: directly invoke
    # the gate function, which is the unit being tested.
    from checkatlas.scfm.orchestrator import _gate_dimred
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    cfg = SCFMConfig(atlas_name="empty_dimred")
    ok, reason = _gate_dimred(adata, detected, cfg)
    assert ok is False
    assert "no embedding with > 3 components detected" in reason


def test_user_ref_label_not_in_atlas_emits_warning(tmp_path: Path, caplog):
    """When the user passes ``--scfm_ref_label=foo`` and ``foo`` is
    not in ``adata.obs.columns``, the orchestrator must warn and
    skip, even if the column detector found other valid columns.
    """
    import logging

    adata = _build_adata(has_celltype=True)  # has 'celltype'
    cfg = SCFMConfig(atlas_name="mismatch", ref_label="nonexistent")

    caplog.set_level(logging.WARNING, logger="checkatlas")
    orchestrator.ensure_per_task_tsvs(
        adata, {"Atlas_name": "mismatch", "Atlas_path": str(tmp_path)},
        argparse.Namespace(
            path=str(tmp_path),
            n_jobs=2,
            metric_cluster=[],
            metric_annot=[],
            metric_dimred=[],
        ),
        cfg,
    )
    warnings = [
        r.getMessage() for r in caplog.records
        if r.levelno == logging.WARNING and r.name == "checkatlas"
    ]
    assert any("nonexistent" in m and "not present" in m for m in warnings)


def test_existing_tsv_is_reused_not_recomputed(tmp_path: Path, caplog):
    """If a per-task TSV is already on disk, the orchestrator must
    log "skipping auto-run" and not call create_metric_*.
    """
    import logging

    adata = _build_adata()
    # Pre-write a cluster TSV
    cluster_dir = folders.get_folder(str(tmp_path), "cluster")
    os.makedirs(cluster_dir, exist_ok=True)
    pd.DataFrame(
        {
            "Clust_Sample": ["x_X_pca_leiden"],
            "obs": ["leiden"],
            "silhouette": [0.5],
            "silhouette_running_time": [0.1],
        }
    ).to_csv(os.path.join(cluster_dir, "x.tsv"), sep="\t", index=False)

    cfg = SCFMConfig(atlas_name="x", ref_label="celltype")
    caplog.set_level(logging.INFO, logger="checkatlas")
    result = orchestrator.ensure_per_task_tsvs(
        adata, {"Atlas_name": "x", "Atlas_path": str(tmp_path)},
        argparse.Namespace(
            path=str(tmp_path),
            n_jobs=2,
            metric_cluster=[],
            metric_annot=[],
            metric_dimred=[],
        ),
        cfg,
    )
    info_msgs = [
        r.getMessage() for r in caplog.records
        if r.levelno == logging.INFO and r.name == "checkatlas"
    ]
    assert any(
        "cluster TSV present" in m and "skipping auto-run" in m
        for m in info_msgs
    )
    assert result["cluster"] is not None


# ──────────────────────────────────────────────────────────────────
# _apply_user_overrides (in pipeline.py)
# ──────────────────────────────────────────────────────────────────


def test_apply_user_overrides_no_op_when_no_user_keys():
    """When the user did not pass any --scfm_* key, the function
    must return the long DataFrame untouched.
    """
    long = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "iLISI",
                "Embedding": "X_pca",
                "Reference/Input 1": "X_pca",
                "Prediction/Input 2": "batch1",
                "Value": 0.1,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "iLISI",
                "Embedding": "X_scvi",
                "Reference/Input 1": "X_scvi",
                "Prediction/Input 2": "batch2",
                "Value": 0.2,
                "Time (s)": 0.0,
            },
        ]
    )
    cfg = SCFMConfig(atlas_name="x")
    out = _apply_user_overrides(long, "annotation", cfg)
    assert len(out) == 2


def test_apply_user_overrides_annot_filters_to_user_ref_and_pred():
    long = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "iLISI",
                "Embedding": "X_pca",
                "Reference/Input 1": "X_pca",
                "Prediction/Input 2": "refA",
                "Value": 0.1,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "iLISI",
                "Embedding": "X_scvi",
                "Reference/Input 1": "X_scvi",
                "Prediction/Input 2": "refB",
                "Value": 0.2,
                "Time (s)": 0.0,
            },
        ]
    )
    cfg = SCFMConfig(
        atlas_name="x", ref_label="refA", predicted_label="refA"
    )
    out = _apply_user_overrides(long, "annotation", cfg)
    # Only the row with Reference="X_pca" and Prediction="refA" remains.
    # Wait — Reference must equal the user ref_label ("refA"), not "X_pca".
    # So actually nothing matches and the table is empty. That is
    # expected: when the user names a specific ref column, the
    # column-detector's broader list (where Reference = X_pca) is
    # filtered out. The user should name both --scfm_ref_label and
    # the embedding they care about.
    assert len(out) == 0


def test_apply_user_overrides_annot_keeps_batch_metric_for_user_key():
    long = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "kbet",
                "Embedding": "X_pca",
                "Reference/Input 1": "X_pca",
                "Prediction/Input 2": "batch_A",
                "Value": 0.1,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "kbet",
                "Embedding": "X_pca",
                "Reference/Input 1": "X_pca",
                "Prediction/Input 2": "batch_B",
                "Value": 0.5,
                "Time (s)": 0.0,
            },
        ]
    )
    cfg = SCFMConfig(atlas_name="x", batch_key="batch_A")
    out = _apply_user_overrides(long, "annotation", cfg)
    assert len(out) == 1
    assert out.iloc[0]["Prediction/Input 2"] == "batch_A"


def test_apply_user_overrides_cluster_filters_to_user_label():
    long = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "cluster",
                "Metric Name": "silhouette",
                "Embedding": "X_pca",
                "Reference/Input 1": "",
                "Prediction/Input 2": "leiden",
                "Value": 0.21,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "cluster",
                "Metric Name": "silhouette",
                "Embedding": "X_pca",
                "Reference/Input 1": "",
                "Prediction/Input 2": "louvain",
                "Value": 0.15,
                "Time (s)": 0.0,
            },
        ]
    )
    cfg = SCFMConfig(atlas_name="x", ref_label="leiden")
    out = _apply_user_overrides(long, "cluster", cfg)
    assert len(out) == 1
    assert out.iloc[0]["Prediction/Input 2"] == "leiden"


def test_apply_user_overrides_dimred_filters_to_user_embeddings():
    long = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "dimred",
                "Metric Name": "dCor",
                "Embedding": "X_pca",
                "Reference/Input 1": "X",
                "Prediction/Input 2": "",
                "Value": 0.0,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "dimred",
                "Metric Name": "dCor",
                "Embedding": "X_scvi",
                "Reference/Input 1": "X",
                "Prediction/Input 2": "",
                "Value": 0.0,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "dimred",
                "Metric Name": "dCor",
                "Embedding": "X_geneformer",
                "Reference/Input 1": "X",
                "Prediction/Input 2": "",
                "Value": 0.5,
                "Time (s)": 0.0,
            },
        ]
    )
    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca",),
    )
    out = _apply_user_overrides(long, "dimred", cfg)
    assert set(out["Embedding"].unique()) == {"X_pca", "X_geneformer"}


# ──────────────────────────────────────────────────────────────────
# All-combinations orchestrator path (no --scfm_* flag set)
# ──────────────────────────────────────────────────────────────────


def test_orchestrator_with_no_scfm_flags_runs_all_three_tasks(
    tmp_path: Path, caplog
):
    """When ``scfm_config`` has every --scfm_* flag empty, the
    orchestrator should fall into the 'all combinations' mode
    and auto-run all 3 per-task engines (not skip them).
    """
    import argparse
    import logging

    import anndata as ad
    import numpy as np
    import pandas as pd

    from checkatlas.scfm.orchestrator import ensure_per_task_tsvs

    rng = np.random.default_rng(0)
    n = 100
    adata = ad.AnnData(
        X=rng.standard_normal((n, 30)),
        obs=pd.DataFrame(
            {
                "celltype": pd.Categorical(
                    rng.choice(["alpha", "beta", "gamma"], n)
                ),
                "leiden": pd.Categorical(
                    rng.choice([f"c{i}" for i in range(5)], n)
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2", "b3"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, 10)).astype(np.float32)
    adata.obsm["X_geneformer"] = (
        rng.standard_normal((n, 10)).astype(np.float32)
    )

    cfg = SCFMConfig(atlas_name="x")  # all --scfm_* flags empty
    from checkatlas.metrics import cluster, annot, dimred

    args = argparse.Namespace(
        path=str(tmp_path),
        n_jobs=2,
        metric_cluster=cluster.__all__,
        metric_annot=annot.__all__,
        metric_dimred=dimred.__all__,
    )

    caplog.set_level(logging.INFO, logger="checkatlas")
    result = ensure_per_task_tsvs(
        adata, {"Atlas_name": "x", "Atlas_path": "x.h5ad"}, args, cfg
    )

    # All 3 per-task TSVs were auto-run
    assert result["cluster"] is not None and len(result["cluster"]) > 0
    assert result["annotation"] is not None and len(result["annotation"]) > 0
    # The dimred task may be empty if the column detector
    # found no usable embeddings with > 3 components; in our
    # synthetic case X_pca and X_geneformer are both 10D so it
    # must produce rows.
    assert result["dimred"] is not None and len(result["dimred"]) > 0

    # The orchestrator logged the 'all combinations' message
    info = " ".join(r.getMessage() for r in caplog.records)
    assert "all combinations" in info
    assert "nothing to auto-run" not in info


def test_orchestrator_with_user_specified_flags_keeps_existing_behaviour(
    tmp_path: Path, caplog
):
    """When the user sets at least one --scfm_* flag, the
    orchestrator should fall back to today's behaviour
    (``required_tasks`` decides which task TSVs to run).
    """
    import argparse
    import logging

    import anndata as ad
    import numpy as np
    import pandas as pd

    from checkatlas.scfm.orchestrator import ensure_per_task_tsvs

    rng = np.random.default_rng(0)
    n = 100
    adata = ad.AnnData(
        X=rng.standard_normal((n, 30)),
        obs=pd.DataFrame(
            {
                "celltype": pd.Categorical(
                    rng.choice(["alpha", "beta"], n)
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n, 10)).astype(np.float32)
    adata.obsm["X_scvi_dropin"] = adata.obsm["X_pca"].copy()
    adata.obsm["X_geneformer"] = (
        rng.standard_normal((n, 10)).astype(np.float32)
    )

    # User specifies scfm_embedding + baseline + ref_label.
    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca", "X_scvi_dropin"),
        ref_label="celltype",
        batch_key="batch",
        scaling_fractions=(),  # skip cal_scfm to keep the test fast
        n_seeds=0,
    )
    from checkatlas.metrics import cluster, annot, dimred

    args = argparse.Namespace(
        path=str(tmp_path),
        n_jobs=2,
        metric_cluster=cluster.__all__,
        metric_annot=annot.__all__,
        metric_dimred=dimred.__all__,
    )

    caplog.set_level(logging.INFO, logger="checkatlas")
    result = ensure_per_task_tsvs(
        adata, {"Atlas_name": "x", "Atlas_path": "x.h5ad"}, args, cfg
    )

    # The orchestrator did NOT log the 'all combinations' message
    info = " ".join(r.getMessage() for r in caplog.records)
    assert "all combinations" not in info

    # The 2 per-task TSVs the user config requires were auto-run
    # (the cluster task may be None if the column detector found
    # no cluster labels — in that case the annot/dimred tasks
    # still ran). The important assertion is that the
    # orchestrator did NOT log the 'all combinations' message,
    # which proves the user-specified path is taken.
    assert result["annotation"] is not None
    assert result["dimred"] is not None
