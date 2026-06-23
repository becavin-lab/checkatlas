"""Diagnostics: from raw metric values to nine ProblemVerdict objects.

This module is the *interpretation* layer. It takes a long-format
metric DataFrame and a ``SCFMConfig``, applies nine rule-based
verdicts (each anchored in a specific paper from the report), and
returns a list of ``ProblemVerdict`` objects with score, grade,
evidence, and citation.

Becavin comment 7: only problems that are *not* already covered by
existing CheckAtlas metrics get their own rule here. Problems 1, 2,
4, 6, 7 are largely re-uses of the existing metric outputs; the
verdicts just summarise them.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd

from .config import SCFMConfig
from .grades import letter
from .rules import load_thresholds

logger = logging.getLogger("checkatlas")


@dataclass
class ProblemVerdict:
    """Result of one diagnostic rule."""

    problem_id: int
    problem_name: str
    score: float
    confidence: str
    evidence: list[str] = field(default_factory=list)
    explanation: str = ""
    reference: str = ""
    grade: str = ""

    def __post_init__(self) -> None:
        if not self.grade:
            self.grade = letter(self.score)


# ──────────────────────────────────────────────────────────────────
# Helpers for extracting metric values from the long-format table
# ──────────────────────────────────────────────────────────────────


def _safe_get(
    df: pd.DataFrame,
    *,
    metric: str,
    embedding: Optional[str] = None,
    task: Optional[str] = None,
    train_domain: Optional[str] = None,
    test_domain: Optional[str] = None,
) -> Optional[float]:
    """Look up a single metric value from the long-format table."""
    if df is None or df.empty:
        return None
    sub = df
    if "Metric Name" in sub.columns:
        sub = sub[sub["Metric Name"] == metric]
    if embedding is not None and "Embedding" in sub.columns:
        sub = sub[sub["Embedding"] == embedding]
    if task is not None and "Task" in sub.columns:
        sub = sub[sub["Task"] == task]
    if train_domain is not None and "Train_Domain" in sub.columns:
        sub = sub[sub["Train_Domain"] == train_domain]
    if test_domain is not None and "Test_Domain" in sub.columns:
        sub = sub[sub["Test_Domain"] == test_domain]
    if sub.empty:
        return None
    try:
        v = float(sub.iloc[0]["Value"])
    except (KeyError, ValueError, TypeError):
        return None
    if v != v:
        return None
    return v


def _conf(score: float, threshold: float, op: str) -> str:
    """Map the magnitude of a breach to a confidence band."""
    try:
        s = float(score)
        t = float(threshold)
    except (TypeError, ValueError):
        return "low"
    if op in ("<", "<="):
        if t == 0:
            return "high" if s < 0 else "low"
        ratio = abs(s - t) / max(abs(t), 1e-9)
    elif op in (">", ">="):
        ratio = abs(s - t) / max(abs(t), 1e-9)
    else:
        ratio = 0.0
    if ratio >= 2.0:
        return "high"
    if ratio >= 1.0:
        return "medium"
    if ratio >= 0.2:
        return "low"
    return "low"


def _verdict_from_threshold(
    value: Optional[float],
    *,
    op: str,
    threshold: float,
    problem_id: int,
    rules: dict[int, dict[str, Any]],
    direction: str = "high_is_bad",
) -> tuple[float, str, str]:
    """Apply a single threshold rule to a single value.

    Returns ``(score, confidence, evidence_metric)`` where ``score``
    is 0..1 (1 = severe breach), capped at 1.0.
    """
    if value is None:
        return float("nan"), "n/a", ""
    rules_entry = rules.get(problem_id, {})
    name = rules_entry.get("name", f"Problem {problem_id}")
    try:
        v = float(value)
        t = float(threshold)
    except (TypeError, ValueError):
        return float("nan"), "n/a", ""
    breach = False
    if op == "<":
        breach = v < t
    elif op == "<=":
        breach = v <= t
    elif op == ">":
        breach = v > t
    elif op == ">=":
        breach = v >= t
    if not breach:
        return 0.0, "high", name
    if direction == "high_is_bad":
        over = (v - t) / max(abs(t), 1e-9)
    else:
        over = (t - v) / max(abs(t), 1e-9)
    score = float(np.clip(0.5 + 0.25 * over, 0.0, 1.0))
    conf = _conf(v, t, op)
    return score, conf, name


# ──────────────────────────────────────────────────────────────────
# The nine rules
# ──────────────────────────────────────────────────────────────────


def _verdict_problem_1(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 1
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Zero-shot performance failure")
    if not config.predicted_label or not config.ref_label:
        return ProblemVerdict(
            p, name, float("nan"), "n/a", [],
            "Predicted and reference labels are required for problem 1.",
            p_rules.get("primary_reference", ""),
        )
    f1 = _safe_get(
        df, metric="isolated_f1_score", embedding=config.scfm_embedding
    )
    ari = _safe_get(
        df, metric="adj_rand_index", embedding=config.scfm_embedding
    )
    scores: list[float] = []
    confidences: list[str] = []
    evidence: list[str] = []
    for metric_name, value, op, thr in [
        ("isolated_f1", f1, "<", 0.50),
        ("adj_rand_index", ari, "<", 0.30),
    ]:
        s, c, n = _verdict_from_threshold(
            value, op=op, threshold=thr, problem_id=p, rules=rules
        )
        if not np.isnan(s):
            scores.append(s)
            confidences.append(c)
            evidence.append(f"{metric_name}={value:.3f} ({n})")
    if not scores:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No ref-vs-pred metrics available for this embedding.",
            p_rules.get("primary_reference", ""),
        )
    return ProblemVerdict(
        p,
        name,
        float(max(scores)),
        "high" if "high" in confidences else ("medium" if "medium" in confidences else "low"),
        evidence,
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_2(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 2
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Baselines outperform scFM")
    if not config.baseline_embeddings or config.scfm_embedding not in df.get("Embedding", pd.Series()).values:
        return ProblemVerdict(p, name, float("nan"), "n/a", [], "Missing baseline or scfm embedding.", p_rules.get("primary_reference", ""))
    metric_cols = [
        "silhouette",
        "average_silhouette_width",
        "isolated_f1_score",
        "adj_rand_index",
        "isolated_f1",
    ]
    deltas: list[float] = []
    evidence: list[str] = []
    for metric_name in metric_cols:
        v_scfm = _safe_get(df, metric=metric_name, embedding=config.scfm_embedding)
        if v_scfm is None:
            continue
        for base in config.baseline_embeddings:
            v_base = _safe_get(df, metric=metric_name, embedding=base)
            if v_base is None:
                continue
            deltas.append(float(v_scfm) - float(v_base))
            evidence.append(f"{metric_name}: {v_scfm:.3f} vs {base}={v_base:.3f}")
    if not deltas:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No paired scFM-vs-baseline metrics available.",
            p_rules.get("primary_reference", ""),
        )
    worst = float(min(deltas))
    score = float(np.clip(0.5 - worst * 5.0, 0.0, 1.0))
    return ProblemVerdict(
        p,
        name,
        score,
        "high" if worst < -0.05 else ("medium" if worst < -0.02 else "low"),
        evidence,
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_3(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 3
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Pretraining scaling laws do not apply")
    plateau = _safe_get(
        df,
        metric="ari",
        embedding=config.scfm_embedding,
        task="scfm_scaling",
    )
    if plateau is None:
        plateau = _safe_get(
            df,
            metric="silhouette",
            embedding=config.scfm_embedding,
            task="scfm_scaling",
        )
    scaling_rows = df[
        (df.get("Task", pd.Series()) == "scfm_scaling")
        & (df.get("Embedding", pd.Series()) == config.scfm_embedding)
    ] if "Task" in df.columns else pd.DataFrame()
    if not scaling_rows.empty and "Fraction" in scaling_rows.columns:
        full = scaling_rows[scaling_rows["Fraction"] >= 0.999]["Value"]
        small = scaling_rows[scaling_rows["Fraction"] <= 0.011]["Value"]
        if not full.empty and not small.empty:
            try:
                plateau = float(full.iloc[0]) / max(abs(float(small.iloc[0])), 1e-9)
            except (ZeroDivisionError, ValueError, TypeError):
                plateau = None
    s, c, n = _verdict_from_threshold(
        plateau, op="<", threshold=1.10, problem_id=p, rules=rules
    )
    return ProblemVerdict(
        p,
        name,
        s,
        c,
        [f"plateau_ratio={plateau:.3f}"] if plateau is not None else [],
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_4(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 4
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Batch effects persist")
    if not config.batch_key:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No --scfm_batch_key provided, batch metrics not evaluable.",
            p_rules.get("primary_reference", ""),
        )
    kbet_v = _safe_get(df, metric="kbet", embedding=config.scfm_embedding)
    ilisi_v = _safe_get(df, metric="iLISI", embedding=config.scfm_embedding)
    pcr_v = _safe_get(df, metric="pcr", embedding=config.scfm_embedding)
    gc_v = _safe_get(
        df, metric="graph_connectivity", embedding=config.scfm_embedding
    )
    scores: list[float] = []
    evidence: list[str] = []
    confidences: list[str] = []
    for metric_name, value, op, thr in [
        ("kbet", kbet_v, ">", 0.50),
        ("ilisi", ilisi_v, "<", 0.30),
        ("pcr", pcr_v, ">", 0.85),
        ("graph_connectivity", gc_v, "<", 0.50),
    ]:
        s, c, n = _verdict_from_threshold(
            value, op=op, threshold=thr, problem_id=p, rules=rules
        )
        if not np.isnan(s) and value is not None:
            scores.append(s)
            confidences.append(c)
            evidence.append(f"{metric_name}={value:.3f} ({n})")
    if not scores:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Batch metrics not computed.",
            p_rules.get("primary_reference", ""),
        )
    score = float(max(scores))
    conf = "high" if "high" in confidences else ("medium" if "medium" in confidences else "low")
    return ProblemVerdict(
        p,
        name,
        score,
        conf,
        evidence,
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_5(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 5
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Tokenization and representation failures")
    gap = _safe_get(
        df, metric="rare_common_gap", embedding=config.scfm_embedding
    )
    dcor = _safe_get(df, metric="dCor", embedding=config.scfm_embedding)
    scores: list[float] = []
    evidence: list[str] = []
    confidences: list[str] = []
    for metric_name, value, op, thr in [
        ("rare_common_gap", gap, ">", 0.15),
        ("dCor_vs_pca", dcor, "<", 0.50),
    ]:
        s, c, n = _verdict_from_threshold(
            value, op=op, threshold=thr, problem_id=p, rules=rules
        )
        if not np.isnan(s) and value is not None:
            scores.append(s)
            confidences.append(c)
            evidence.append(f"{metric_name}={value:.3f} ({n})")
    if not scores:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Tokenization metrics not available.",
            p_rules.get("primary_reference", ""),
        )
    return ProblemVerdict(
        p,
        name,
        float(max(scores)),
        "high" if "high" in confidences else ("medium" if "medium" in confidences else "low"),
        evidence,
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_6(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 6
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Generalization failures across domains")
    if not config.domain_key or "Task" not in df.columns:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No --scfm_domain_key provided, cross-domain metrics not evaluable.",
            p_rules.get("primary_reference", ""),
        )
    cross = df[
        (df["Task"] == "scfm_cross_domain")
        & (df["Embedding"] == config.scfm_embedding)
    ]
    if cross.empty:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Cross-domain metrics not computed (not enough domains).",
            p_rules.get("primary_reference", ""),
        )
    vals = cross["Value"].dropna()
    if vals.empty:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Cross-domain values are all NaN.",
            p_rules.get("primary_reference", ""),
        )
    mean_v = float(vals.mean())
    s, c, n = _verdict_from_threshold(
        mean_v, op="<", threshold=0.30, problem_id=p, rules=rules
    )
    return ProblemVerdict(
        p,
        name,
        s,
        c,
        [f"cross_domain_ari_mean={mean_v:.3f} ({n})"],
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_7(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 7
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Real-world clinical / translational failures")
    if not config.patient_key:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No --scfm_patient_key provided, patient metrics not evaluable.",
            p_rules.get("primary_reference", ""),
        )
    return ProblemVerdict(
        p,
        name,
        0.0,
        "low",
        [],
        (
            "Patient-level metrics (ASW, outcome AUC) are reported "
            "in per_metric.tsv; cross-check against Elmarakeby "
            "2025 thresholds."
        ),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_8(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 8
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Model stability and interpretability issues")
    if "Task" not in df.columns:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Stability table not in metric DataFrame.",
            p_rules.get("primary_reference", ""),
        )
    stab = df[
        (df["Task"] == "scfm_stability")
        & (df["Embedding"] == config.scfm_embedding)
        & (df["Metric Name"].astype(str).str.endswith("_cv"))
    ]
    if stab.empty:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Stability CV values not available.",
            p_rules.get("primary_reference", ""),
        )
    vals = stab["Value"].dropna()
    if vals.empty:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "Stability CV values are all NaN.",
            p_rules.get("primary_reference", ""),
        )
    mean_cv = float(vals.mean())
    s, c, n = _verdict_from_threshold(
        mean_cv, op=">", threshold=0.10, problem_id=p, rules=rules
    )
    return ProblemVerdict(
        p,
        name,
        s,
        c,
        [f"stability_cv_mean={mean_cv:.3f} ({n})"],
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


def _verdict_problem_9(
    df: pd.DataFrame, config: SCFMConfig, rules: dict[int, dict[str, Any]]
) -> ProblemVerdict:
    p = 9
    p_rules = rules.get(p, {})
    name = p_rules.get("name", "Denoising hypothesis (linear manifold)")
    if not config.baseline_embeddings:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No baseline embeddings provided, denoiser verdict not evaluable.",
            p_rules.get("primary_reference", ""),
        )
    metric_cols = [
        "silhouette",
        "average_silhouette_width",
        "isolated_f1_score",
        "adj_rand_index",
        "isolated_f1",
    ]
    better = 0
    total = 0
    evidence: list[str] = []
    for metric_name in metric_cols:
        v_scfm = _safe_get(df, metric=metric_name, embedding=config.scfm_embedding)
        if v_scfm is None:
            continue
        for base in config.baseline_embeddings:
            v_base = _safe_get(df, metric=metric_name, embedding=base)
            if v_base is None:
                continue
            total += 1
            if v_scfm > v_base + 0.01:
                better += 1
            evidence.append(
                f"{metric_name}: scfm={v_scfm:.3f} vs {base}={v_base:.3f}"
            )
    if total == 0:
        return ProblemVerdict(
            p,
            name,
            float("nan"),
            "n/a",
            [],
            "No paired scFM-vs-baseline metrics.",
            p_rules.get("primary_reference", ""),
        )
    competitiveness = better / total
    s, c, n = _verdict_from_threshold(
        competitiveness, op="<", threshold=0.20, problem_id=p, rules=rules
    )
    return ProblemVerdict(
        p,
        name,
        s,
        c,
        [
            f"baseline_competitiveness={competitiveness:.2f} "
            f"({better}/{total} metrics)"
        ],
        p_rules.get("explanation", ""),
        p_rules.get("primary_reference", ""),
    )


# ──────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────


_RULE_FUNCS = [
    _verdict_problem_1,
    _verdict_problem_2,
    _verdict_problem_3,
    _verdict_problem_4,
    _verdict_problem_5,
    _verdict_problem_6,
    _verdict_problem_7,
    _verdict_problem_8,
    _verdict_problem_9,
]


def diagnose(
    metrics_df: pd.DataFrame,
    config: SCFMConfig,
    *,
    thresholds_path: Optional[str] = None,
) -> list[ProblemVerdict]:
    """Apply the nine diagnostic rules to a long-format metric table.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        The combined long-format table from ``cal_scfm`` and the
        existing ``cal_annot/cluster/dimred`` outputs. Required
        columns: ``Metric Name``, ``Embedding``, ``Value``. Optional:
        ``Task``, ``Train_Domain``, ``Test_Domain``, ``Fraction``.
    config : SCFMConfig
    thresholds_path : str, optional
        User-supplied YAML of threshold overrides.

    Returns
    -------
    list[ProblemVerdict]
        One entry per problem (always 9, in order 1..9).
    """
    rules = load_thresholds(thresholds_path or config.thresholds_path)
    verdicts: list[ProblemVerdict] = []
    for fn in _RULE_FUNCS:
        try:
            verdicts.append(fn(metrics_df, config, rules))
        except Exception as exc:  # pragma: no cover
            logger.warning("diagnose: rule %s failed: %s", fn.__name__, exc)
            verdicts.append(
                ProblemVerdict(
                    problem_id=int(fn.__name__.split("_")[-1]),
                    problem_name=fn.__name__,
                    score=float("nan"),
                    confidence="n/a",
                    evidence=[],
                    explanation=f"Rule failed: {exc}",
                    reference="",
                )
            )
    return verdicts
