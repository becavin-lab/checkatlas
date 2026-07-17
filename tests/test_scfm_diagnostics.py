"""Tests for the scFM QC layer.

These tests cover the *interpretation* layer (rules, grades, composite,
diagnostics) and the report writers. They use synthetic long-format
metric DataFrames; no scanpy / sklearn / anndata required, so the
tests run in any environment that has numpy + pandas.

The metric-collection layer (scaling, stability, rare_types,
cross_domain) is exercised separately in
``tests/test_scfm_metrics.py`` and skips gracefully when its
optional dependencies (scanpy, sklearn, anndata) are missing.
"""

import json
import os

import numpy as np
import pandas as pd
import pytest


# ──────────────────────────────────────────────────────────────────
# grades
# ──────────────────────────────────────────────────────────────────


def test_grade_letter_bands():
    from checkatlas.scfm.grades import letter

    assert letter(0.0) == "A"
    assert letter(0.10) == "A"
    assert letter(0.149) == "A"
    assert letter(0.15) == "B"
    assert letter(0.34) == "B"
    assert letter(0.35) == "C"
    assert letter(0.59) == "C"
    assert letter(0.60) == "D"
    assert letter(0.79) == "D"
    assert letter(0.80) == "F"
    assert letter(1.0) == "F"


def test_grade_letter_nan():
    from checkatlas.scfm.grades import letter

    assert letter(float("nan")) == "n/a"
    assert letter(None) == "n/a"
    assert letter("not-a-number") == "n/a"


def test_grade_legend_is_markdown():
    from checkatlas.scfm.grades import grade_legend

    md = grade_legend()
    assert "# scFM problem grade legend" in md
    assert "| A |" in md
    assert "| F |" in md


# ──────────────────────────────────────────────────────────────────
# rules
# ──────────────────────────────────────────────────────────────────


def test_default_rules_cover_all_nine_problems():
    from checkatlas.scfm.rules import DEFAULT_RULES

    assert set(DEFAULT_RULES.keys()) == {1, 2, 3, 4, 5, 6, 7, 8, 9}
    for pid, rule in DEFAULT_RULES.items():
        assert "name" in rule
        assert "metric_thresholds" in rule
        assert len(rule["metric_thresholds"]) > 0
        for entry in rule["metric_thresholds"]:
            assert len(entry) == 4
            metric, op, thr, ref = entry
            assert op in ("<", "<=", ">", ">=")
            assert isinstance(thr, (int, float))
            assert ref


def test_load_thresholds_default_when_no_path():
    from checkatlas.scfm.rules import DEFAULT_RULES, load_thresholds

    rules = load_thresholds(None)
    assert rules.keys() == DEFAULT_RULES.keys()


def test_load_thresholds_default_when_file_missing():
    from checkatlas.scfm.rules import DEFAULT_RULES, load_thresholds

    rules = load_thresholds("/tmp/this-file-does-not-exist.yaml")
    assert rules.keys() == DEFAULT_RULES.keys()


def test_write_resolved_thresholds(tmp_path):
    from checkatlas.scfm.rules import (
        DEFAULT_RULES,
        write_resolved_thresholds,
    )

    out = tmp_path / "thr.yaml"
    path = write_resolved_thresholds(DEFAULT_RULES, str(out))
    assert os.path.exists(path)
    assert path == str(out)
    # The fallback plain-text format is used when PyYAML is missing,
    # but the file must exist in any case.
    assert os.path.getsize(path) > 0


# ──────────────────────────────────────────────────────────────────
# composite
# ──────────────────────────────────────────────────────────────────


def test_default_weights_has_three_sections():
    from checkatlas.scfm.composite import DEFAULT_WEIGHTS

    assert set(DEFAULT_WEIGHTS.keys()) == {"fmf", "bf", "pr"}
    for section in ("fmf", "bf", "pr"):
        for k, v in DEFAULT_WEIGHTS[section].items():
            assert isinstance(v, (int, float))


def test_load_weights_user_override(tmp_path):
    from checkatlas.scfm.composite import DEFAULT_WEIGHTS, load_weights

    user = {"fmf": {"batch": 5.0}}
    p = tmp_path / "w.json"
    p.write_text(json.dumps(user))
    weights = load_weights(str(p))
    assert weights["fmf"]["batch"] == 5.0
    assert weights["fmf"]["rare_types"] == DEFAULT_WEIGHTS["fmf"]["rare_types"]


def test_load_weights_ignores_unknown_keys(tmp_path):
    from checkatlas.scfm.composite import load_weights

    p = tmp_path / "w.json"
    p.write_text(json.dumps({"fmf": {"batch": 2.0, "nonsense_key": 99.0}}))
    weights = load_weights(str(p))
    assert "nonsense_key" not in weights["fmf"]
    assert weights["fmf"]["batch"] == 2.0


def test_fmf_monotonic_in_score():
    from checkatlas.scfm.composite import DEFAULT_WEIGHTS, fmf

    def verdicts_with_uniform_score(s):
        from checkatlas.scfm.diagnostics import ProblemVerdict

        return [
            ProblemVerdict(
                problem_id=pid,
                problem_name=f"Problem {pid}",
                score=s,
                confidence="low",
                grade="",
            )
            for pid in range(1, 10)
        ]

    scores = [fmf(verdicts_with_uniform_score(s), DEFAULT_WEIGHTS) for s in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert scores == sorted(scores, reverse=True)
    assert scores[0] == pytest.approx(100.0, abs=1e-6)
    assert scores[-1] == pytest.approx(0.0, abs=1e-6)


def test_fmf_geometric_mean_pulls_down_on_one_failure():
    from checkatlas.scfm.composite import DEFAULT_WEIGHTS, fmf
    from checkatlas.scfm.diagnostics import ProblemVerdict

    base = [
        ProblemVerdict(pid, f"p{pid}", 0.0, "low") for pid in range(1, 10)
    ]
    one_bad = list(base)
    one_bad[3] = ProblemVerdict(4, "Batch effects persist", 1.0, "low")
    full_good = fmf(base, DEFAULT_WEIGHTS)
    one_bad_fmf = fmf(one_bad, DEFAULT_WEIGHTS)
    assert full_good > one_bad_fmf
    assert one_bad_fmf < 50.0


def test_bf_in_unit_interval():
    from checkatlas.scfm.composite import bf

    for clisi in (0, 0.5, 1):
        for gc in (0, 0.5, 1):
            for ilisi in (0, 0.5, 1):
                v = bf(
                    clisi_norm=clisi,
                    graph_connectivity_norm=gc,
                    ilisi_penalty_norm=ilisi,
                )
                assert 0.0 <= v <= 1.0


def test_pr_in_unit_interval():
    from checkatlas.scfm.composite import pr

    for stab in (0, 0.5, 1):
        for tok in (0, 0.5, 1):
            for fr in (0, 0.5, 1):
                v = pr(
                    stability_term=stab,
                    tokenization_term=tok,
                    failure_rate=fr,
                )
                assert 0.0 <= v <= 1.0


def test_compute_all_returns_three_keys():
    from checkatlas.scfm.composite import compute_all
    from checkatlas.scfm.diagnostics import ProblemVerdict

    verdicts = [
        ProblemVerdict(pid, f"p{pid}", 0.3, "low") for pid in range(1, 10)
    ]
    out = compute_all(verdicts)
    assert set(out.keys()) == {"fmf", "bf", "pr"}


# ──────────────────────────────────────────────────────────────────
# config
# ──────────────────────────────────────────────────────────────────


def test_scfmconfig_defaults():
    from checkatlas.scfm.config import SCFMConfig

    c = SCFMConfig()
    assert c.atlas_name == ""
    assert c.scfm_embedding == ""
    assert c.baseline_embeddings == ()
    assert c.scaling_fractions == (0.01, 0.10, 0.50, 1.00)
    assert c.n_seeds == 5
    assert c.weights_path is None


def test_scfmconfig_all_embeddings():
    from checkatlas.scfm.config import SCFMConfig

    c = SCFMConfig(
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca", "X_scvi", "X_geneformer"),
    )
    embeds = c.all_embeddings()
    assert embeds[0] == "X_geneformer"
    assert "X_pca" in embeds
    assert "X_scvi" in embeds
    assert len(embeds) == 3  # dedup


def test_from_args_with_namespace():
    from argparse import Namespace
    from checkatlas.scfm.config import from_args

    args = Namespace(
        atlas_name="x",
        scfm_embedding="X_gpt",
        baseline_embeddings=["X_pca"],
        scfm_predicted_label="pred",
        scfm_ref_label="ref",
        scfm_batch_key="b",
        scfm_domain_key="d",
        scfm_patient_key="p",
        scfm_outcome_key="o",
        scfm_scaling_fractions=[0.1, 1.0],
        scfm_n_seeds=3,
        scfm_noise_sigma=0.2,
        scfm_min_domain_size=20,
        scfm_weights=None,
        scfm_thresholds=None,
    )
    c = from_args(args)
    assert c.scfm_embedding == "X_gpt"
    assert c.scaling_fractions == (0.1, 1.0)
    assert c.n_seeds == 3
    assert c.noise_sigma == 0.2


# ──────────────────────────────────────────────────────────────────
# diagnostics
# ──────────────────────────────────────────────────────────────────


def _make_synth_metrics():
    return pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "annot",
                "Metric Name": "isolated_f1_score",
                "Reference/Input 1": "X_scfm",
                "Prediction/Input 2": "ref",
                "Embedding": "X_scfm",
                "Value": 0.30,
                "Time (s)": 0.1,
            },
            {
                "Atlas Name": "x",
                "Task": "batch_correction",
                "Metric Name": "kbet",
                "Embedding": "X_scfm",
                "Value": 0.55,
                "Time (s)": 0.1,
            },
            {
                "Atlas Name": "x",
                "Task": "batch_correction",
                "Metric Name": "iLISI",
                "Embedding": "X_scfm",
                "Value": 0.20,
                "Time (s)": 0.1,
            },
            {
                "Atlas Name": "x",
                "Task": "scfm_stability",
                "Metric Name": "silhouette_cv",
                "Embedding": "X_scfm",
                "Value": 0.04,
                "Time (s)": 0.0,
            },
            {
                "Atlas Name": "x",
                "Task": "scfm_stability",
                "Metric Name": "ari_cv",
                "Embedding": "X_scfm",
                "Value": 0.05,
                "Time (s)": 0.0,
            },
        ]
    )


def test_diagnose_returns_nine_verdicts():
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.diagnostics import diagnose

    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_scfm",
        baseline_embeddings=("X_pca",),
        predicted_label="pred",
        ref_label="ref",
        batch_key="b",
    )
    verdicts = diagnose(_make_synth_metrics(), cfg)
    assert len(verdicts) == 9
    assert [v.problem_id for v in verdicts] == list(range(1, 10))
    for v in verdicts:
        assert v.grade in {"A", "B", "C", "D", "F", "n/a"}


def test_diagnose_missing_inputs_yield_na():
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.diagnostics import diagnose

    cfg = SCFMConfig()  # no labels, no batch key, no domain key
    verdicts = diagnose(pd.DataFrame(), cfg)
    assert all(v.confidence == "n/a" for v in verdicts)


def test_diagnose_problem_4_triggers_on_kbet_and_ilisi():
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.diagnostics import diagnose

    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_scfm",
        batch_key="b",
        ref_label="ref",
        predicted_label="pred",
    )
    df = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "batch_correction",
                "Metric Name": "kbet",
                "Embedding": "X_scfm",
                "Value": 0.70,
            },
            {
                "Atlas Name": "x",
                "Task": "batch_correction",
                "Metric Name": "iLISI",
                "Embedding": "X_scfm",
                "Value": 0.20,
            },
        ]
    )
    verdicts = diagnose(df, cfg)
    p4 = next(v for v in verdicts if v.problem_id == 4)
    assert p4.score > 0
    assert p4.grade in {"C", "D", "F"}


def test_diagnose_problem_8_high_cv():
    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.diagnostics import diagnose

    cfg = SCFMConfig(
        atlas_name="x",
        scfm_embedding="X_scfm",
        batch_key="b",
    )
    df = pd.DataFrame(
        [
            {
                "Atlas Name": "x",
                "Task": "scfm_stability",
                "Metric Name": "silhouette_cv",
                "Embedding": "X_scfm",
                "Value": 0.25,
            },
            {
                "Atlas Name": "x",
                "Task": "scfm_stability",
                "Metric Name": "ari_cv",
                "Embedding": "X_scfm",
                "Value": 0.30,
            },
        ]
    )
    verdicts = diagnose(df, cfg)
    p8 = next(v for v in verdicts if v.problem_id == 8)
    assert p8.score > 0
    assert p8.grade in {"C", "D", "F"}


# ──────────────────────────────────────────────────────────────────
# report
# ──────────────────────────────────────────────────────────────────


def test_report_writes_all_files(tmp_path):
    from checkatlas.scfm.diagnostics import ProblemVerdict
    from checkatlas.scfm.report import write_all

    verdicts = [
        ProblemVerdict(pid, f"p{pid}", 0.3, "low", [f"m{pid}"], "ex", "ref")
        for pid in range(1, 10)
    ]
    metrics_df = _make_synth_metrics()
    out = tmp_path / "out"
    out.mkdir()
    paths = write_all(
        atlas_name="x",
        verdicts=verdicts,
        metrics_df=metrics_df,
        outdir=str(out),
        config_dict={"atlas": "x", "scfm_embedding": "X_scfm"},
    )
    for filename in (
        "verdicts.tsv",
        "composite.tsv",
        "per_metric.tsv",
        "inputs.tsv",
        "grade_legend.md",
        "resolved_weights.json",
        "resolved_thresholds.yaml",
    ):
        assert filename in paths
        assert os.path.exists(paths[filename]), f"{filename} missing"
        assert os.path.getsize(paths[filename]) > 0


def test_report_verdicts_have_nine_rows(tmp_path):
    from checkatlas.scfm.diagnostics import ProblemVerdict
    from checkatlas.scfm.report import write_verdicts

    verdicts = [
        ProblemVerdict(pid, f"p{pid}", 0.5, "low") for pid in range(1, 10)
    ]
    out = tmp_path / "out"
    out.mkdir()
    path = write_verdicts(verdicts, "x", str(out))
    df = pd.read_csv(path, sep="\t")
    assert len(df) == 9
    assert list(df["problem_id"]) == list(range(1, 10))
    assert "grade" in df.columns
    assert "score" in df.columns
    assert "explanation" in df.columns
    assert "reference" in df.columns


def test_report_composite_one_row_per_atlas(tmp_path):
    from checkatlas.scfm.diagnostics import ProblemVerdict
    from checkatlas.scfm.report import write_composite

    verdicts = [
        ProblemVerdict(pid, f"p{pid}", 0.3, "low") for pid in range(1, 10)
    ]
    out = tmp_path / "out"
    out.mkdir()
    path = write_composite(
        "x", verdicts, {"fmf": {}, "bf": {}, "pr": {}}, str(out)
    )
    df = pd.read_csv(path, sep="\t")
    assert len(df) == 1
    assert "fmf" in df.columns
    assert "bf" in df.columns
    assert "pr" in df.columns
    for pid in range(1, 10):
        assert f"problem_{pid}_score" in df.columns
        assert f"problem_{pid}_grade" in df.columns


# ──────────────────────────────────────────────────────────────────
# column detector: scFM embedding auto-detection (Becavin comment 1)
# ──────────────────────────────────────────────────────────────────


def test_column_detector_finds_scfm_embeddings():
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    class _FakeAnnData:
        def __init__(self, keys):
            self.obsm = {k: np.zeros((10, 5)) for k in keys}
            self.var = type("v", (), {"columns": pd.Index([])})()
            self.obs = pd.DataFrame()

    for key in (
        "X_geneformer",
        "X_scgpt",
        "X_uce",
        "X_scfoundation",
        "X_sccello",
        "X_scbert",
    ):
        adata = _FakeAnnData([key])
        detector = CheckAtlasColumnDetector(adata)
        matched = detector.analyze_obsm_semantics()
        names = [k for k, _ in matched]
        assert key in names, f"{key} not auto-detected"
