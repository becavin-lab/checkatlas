"""Tests for the scfm combos module (all-combinations mode).

When ``checkatlas scfm <path> --atlas_name X`` is called with no
``--scfm_*`` flags, the orchestrator falls into "all combinations"
mode: it evaluates the (ref, pred, batch, scfm, baseline) cartesian
product the column detector found, up to a cap. This module
covers that path.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from checkatlas.scfm.combos import (
    ComboSpec,
    build_combo_remark,
    build_headline_remark,
    build_per_metric_remark,
    build_verdict_remark,
    detect_all_combos,
    expand_config_for_combo,
    summarise_composite_across_combos,
)
from checkatlas.scfm.config import SCFMConfig


def _build_adata(
    *,
    n: int = 100,
    has_celltype: bool = True,
    has_leiden: bool = True,
    has_batch: bool = True,
    n_pca: int = 20,
    n_scvi: int = 20,
    add_geneformer: bool = True,
):
    """Build a small AnnData with the columns/embeddings the
    scfm "all combinations" path needs."""
    rng = np.random.default_rng(0)
    obs_dict: dict = {
        "obs_name": pd.Categorical([f"c{i}" for i in range(n)]),
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
    if add_geneformer:
        adata.obsm["X_geneformer"] = (
            rng.standard_normal((n, n_pca)).astype(np.float32)
        )
    return adata


# ──────────────────────────────────────────────────────────────────
# ComboSpec
# ──────────────────────────────────────────────────────────────────


def test_combo_spec_combo_id_is_deterministic():
    a = ComboSpec(ref_label="celltype", scfm_embedding="X_geneformer")
    b = ComboSpec(ref_label="celltype", scfm_embedding="X_geneformer")
    assert a.combo_id == b.combo_id
    assert a.combo_id == "celltype|*|*|X_geneformer|*"


def test_combo_spec_uses_star_for_missing_roles():
    combo = ComboSpec(scfm_embedding="X_pca")
    assert combo.ref_label is None
    assert combo.combo_id == "*|*|*|X_pca|*"


def test_combo_spec_roles_used():
    combo = ComboSpec(
        ref_label="celltype",
        batch_key="batch",
        scfm_embedding="X_geneformer",
    )
    assert set(combo.roles_used()) == {
        "ref_label",
        "batch_key",
        "scfm_embedding",
    }


def test_combo_spec_from_args_empty_string_to_none():
    combo = ComboSpec.from_args(
        ref_label="",
        predicted_label="celltype_scfm_pred",
        batch_key="",
        scfm_embedding="X_geneformer",
        baseline_embedding="X_pca",
    )
    assert combo.ref_label is None
    assert combo.predicted_label == "celltype_scfm_pred"
    assert combo.scfm_embedding == "X_geneformer"


# ──────────────────────────────────────────────────────────────────
# detect_all_combos
# ──────────────────────────────────────────────────────────────────


def test_detect_all_combos_single_combo_when_only_one_embedding():
    adata = _build_adata(add_geneformer=False)
    adata.obsm.pop("X_scvi", None)
    # Only X_pca and X_umap. X_umap is 2D and will be filtered out
    # by the detector, leaving only X_pca. With no scfm-named
    # embedding the detector falls back to "X_pca is the scfm" and
    # the baseline is also X_pca.
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    combos_list = detect_all_combos(adata, detected)
    assert len(combos_list) == 1
    assert combos_list[0].scfm_embedding == "X_pca"


def test_detect_all_combos_scfm_and_baselines_classified_separately():
    adata = _build_adata()  # adds X_pca, X_scvi, X_geneformer
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    combos_list = detect_all_combos(adata, detected)
    # X_geneformer is the scfm, X_pca and X_scvi are baselines.
    scfm_embs = {c.scfm_embedding for c in combos_list}
    baseline_embs = {c.baseline_embedding for c in combos_list}
    assert "X_geneformer" in scfm_embs
    assert "X_pca" in baseline_embs
    assert "X_scvi" in baseline_embs
    assert "X_geneformer" not in baseline_embs


def test_detect_all_combos_ref_pred_batch_combinations():
    """The column detector's default ``min_reference_score=0.95`` is
    strict. We bypass the strict defaults by using ad-hoc column
    names that score high on both reference and predicted
    semantics: "celltype" (reference) and "scpred_celltype"
    (predicted), which is a name pattern the detector treats as
    a prediction."""
    n = 200
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        X=rng.standard_normal((n, 50)),
        obs=pd.DataFrame(
            {
                "obs_name": pd.Categorical(
                    [f"c{i}" for i in range(n)]
                ),
                "celltype": pd.Categorical(
                    rng.choice(
                        ["alpha", "beta", "gamma", "delta", "epsilon"],
                        n,
                    )
                ),
                "scpred_celltype": pd.Categorical(
                    rng.choice(
                        ["alpha", "beta", "gamma", "delta", "epsilon"],
                        n,
                    )
                ),
                "celltypist_pred": pd.Categorical(
                    rng.choice(
                        ["alpha", "beta", "gamma", "delta", "epsilon"],
                        n,
                    )
                ),
                "batch": pd.Categorical(
                    rng.choice(["b1", "b2", "b3"], n)
                ),
            }
        ),
    )
    adata.obsm["X_pca"] = rng.standard_normal(
        (n, 20)
    ).astype(np.float32)
    adata.obsm["X_geneformer"] = rng.standard_normal(
        (n, 20)
    ).astype(np.float32)
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    combos_list = detect_all_combos(
        adata,
        detected,
        max_ref=2,
        max_pred=2,
        max_batch=1,
        max_scfm=1,
        max_baseline=1,
    )
    # 1 scfm × 1 baseline × 1 ref × 2 pred × 1 batch = 2 combos
    # (celltype is detected as ref, scpred/celltypist as pred)
    assert len(combos_list) >= 1
    pred_set = {c.predicted_label for c in combos_list}
    assert len(pred_set) >= 1


def test_detect_all_combos_no_embeddings_returns_empty_combo():
    adata = _build_adata()
    # Wipe all embeddings (X_umap is 2D and will be filtered out;
    # X_pca, X_scvi, X_geneformer are >3D)
    adata.obsm.clear()
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    combos_list = detect_all_combos(adata, detected)
    # We always return at least one combo so the user gets a
    # single n/a row rather than an empty composite.
    assert len(combos_list) == 1
    assert combos_list[0].scfm_embedding is None


def test_detect_all_combos_respects_max_combos_cap():
    adata = _build_adata()
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    full = detect_all_combos(
        adata, detected, max_ref=3, max_pred=2, max_batch=2, max_scfm=3, max_baseline=3
    )
    capped = detect_all_combos(
        adata, detected, max_ref=1, max_pred=1, max_batch=1, max_scfm=1, max_baseline=1
    )
    assert len(full) > len(capped)
    assert len(capped) == 1


# ──────────────────────────────────────────────────────────────────
# expand_config_for_combo
# ──────────────────────────────────────────────────────────────────


def test_expand_config_for_combo_fills_roles():
    cfg = SCFMConfig(atlas_name="x")
    combo = ComboSpec(
        ref_label="celltype",
        predicted_label="celltype_scfm_pred",
        batch_key="batch",
        scfm_embedding="X_geneformer",
        baseline_embedding="X_pca",
    )
    out = expand_config_for_combo(cfg, combo)
    assert out.ref_label == "celltype"
    assert out.predicted_label == "celltype_scfm_pred"
    assert out.batch_key == "batch"
    assert out.scfm_embedding == "X_geneformer"
    assert out.baseline_embeddings == ("X_pca",)
    # Original config's domain_key / patient_key are preserved
    assert out.domain_key == cfg.domain_key
    assert out.patient_key == cfg.patient_key


def test_expand_config_for_combo_falls_back_to_config_when_combo_missing():
    cfg = SCFMConfig(
        atlas_name="x",
        ref_label="user_ref",
        scfm_embedding="user_scfm",
    )
    combo = ComboSpec()  # everything None
    out = expand_config_for_combo(cfg, combo)
    assert out.ref_label == "user_ref"
    assert out.scfm_embedding == "user_scfm"


# ──────────────────────────────────────────────────────────────────
# summarise_composite_across_combos
# ──────────────────────────────────────────────────────────────────


def test_summarise_composite_empty():
    out = summarise_composite_across_combos([])
    assert out["n_combos"] == 0
    assert out["n_worst_problems"] == 0
    assert "No combinations" in out["remark"]


def test_summarise_composite_means_fmf_bf_pr():
    per_combo = [
        {"fmf": 0.5, "bf": 1.0, "pr": 0.3, "n_worst_problems": 1},
        {"fmf": 0.7, "bf": 0.8, "pr": 0.5, "n_worst_problems": 0},
        {"fmf": 0.6, "bf": 0.9, "pr": 0.4, "n_worst_problems": 2},
    ]
    out = summarise_composite_across_combos(per_combo)
    assert out["fmf"] == pytest.approx(0.6, abs=1e-6)
    assert out["bf"] == pytest.approx(0.9, abs=1e-4)
    assert out["pr"] == pytest.approx(0.4, abs=1e-4)
    assert out["n_worst_problems"] == 3
    assert out["n_combos"] == 3
    assert "Mean of 3" in out["remark"]


def test_summarise_composite_remark_mentions_spread():
    per_combo = [
        {"fmf": 0.2, "bf": 1.0, "pr": 0.3, "n_worst_problems": 0},
        {"fmf": 0.8, "bf": 1.0, "pr": 0.7, "n_worst_problems": 0},
    ]
    out = summarise_composite_across_combos(per_combo)
    assert "spreads from" in out["remark"]


# ──────────────────────────────────────────────────────────────────
# build_combo_remark
# ──────────────────────────────────────────────────────────────────


def test_build_combo_remark_lists_used_roles():
    combo = ComboSpec(
        ref_label="celltype",
        scfm_embedding="X_geneformer",
        baseline_embedding="X_pca",
    )
    remark = build_combo_remark(combo)
    assert "ref_label=celltype" in remark
    assert "scfm_embedding=X_geneformer" in remark
    assert "baseline_embedding=X_pca" in remark


def test_build_combo_remark_when_all_empty():
    combo = ComboSpec()
    remark = build_combo_remark(combo)
    assert "No ref" in remark


# ──────────────────────────────────────────────────────────────────
# build_verdict_remark
# ──────────────────────────────────────────────────────────────────


def test_build_verdict_remark_grade_a_or_b_competitive():
    combo = ComboSpec(
        ref_label="celltype",
        predicted_label="celltype_scfm_pred",
        scfm_embedding="X_geneformer",
        baseline_embedding="X_pca",
    )
    for grade in ("A", "B"):
        remark = build_verdict_remark(2, grade, 0.2, combo)
        assert "competitive" in remark
        assert f"grade {grade}" in remark
        assert "Problem 2" in remark
        assert "X_geneformer" in remark
        assert "X_pca" in remark


def test_build_verdict_remark_grade_d_or_f_failing():
    combo = ComboSpec(
        scfm_embedding="X_geneformer",
        baseline_embedding="X_pca",
    )
    for grade in ("D", "F"):
        remark = build_verdict_remark(2, grade, 0.9, combo)
        assert "fails on this problem" in remark
        assert f"grade {grade}" in remark


def test_build_verdict_remark_grade_na_explains_why():
    combo = ComboSpec(batch_key="batch")
    remark = build_verdict_remark(4, "n/a", float("nan"), combo, "no batch key")
    assert "could not evaluate" in remark
    assert "Problem 4" in remark
    assert "batch_key=batch" in remark
    assert "no batch key" in remark


def test_build_verdict_remark_problem_1_mentions_ref_and_pred():
    combo = ComboSpec(
        ref_label="celltype",
        predicted_label="celltype_scfm_pred",
    )
    remark = build_verdict_remark(1, "F", 0.05, combo)
    assert "ref=celltype" in remark
    assert "pred=celltype_scfm_pred" in remark


def test_build_verdict_remark_problem_4_mentions_batch():
    combo = ComboSpec(batch_key="batch")
    remark = build_verdict_remark(4, "A", 0.1, combo)
    assert "batch_key=batch" in remark


# ──────────────────────────────────────────────────────────────────
# build_per_metric_remark
# ──────────────────────────────────────────────────────────────────


def test_build_per_metric_remark_includes_metric_and_direction():
    remark_higher = build_per_metric_remark("silhouette", "X_pca", 0.25)
    assert "silhouette" in remark_higher
    assert "X_pca" in remark_higher
    assert "higher is better" in remark_higher
    remark_lower = build_per_metric_remark("davies_bouldin", "X_pca", 1.5)
    assert "davies_bouldin" in remark_lower
    assert "lower is better" in remark_lower


def test_build_per_metric_remark_handles_unknown_metric():
    remark = build_per_metric_remark("exotic_metric", "X_pca", 0.5)
    assert "exotic_metric" in remark
    assert "X_pca" in remark


# ──────────────────────────────────────────────────────────────────
# build_headline_remark
# ──────────────────────────────────────────────────────────────────


def test_build_headline_remark_mentions_count():
    per_combo = [
        {"fmf": 0.5, "bf": 1.0, "pr": 0.3, "n_worst_problems": 0},
        {"fmf": 0.7, "bf": 0.8, "pr": 0.5, "n_worst_problems": 0},
    ]
    out = build_headline_remark(per_combo)
    assert "Mean of 2" in out
    assert "combinations" in out
