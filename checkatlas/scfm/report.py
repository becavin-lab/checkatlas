"""Reporting: write the scFM QC outputs.

Three TSV files, all idempotent, written flat into the per-atlas
``checkatlas_files/scfm/<atlas_name>/`` folder (so the layout matches
the per-task ``checkatlas_files/<task>/<atlas>.tsv`` pattern, with
the atlas name in the path rather than the filename):

  * ``verdicts.tsv``        — one row per (atlas, problem) or per
    (atlas, problem, combo) in all-combinations mode
  * ``composite.tsv``       — one row per atlas (single combo) or
    N+1 rows in all-combinations mode (one per combo plus a
    ``headline`` row that averages across combos)
  * ``per_metric.tsv``      — long-format raw metric values, with
    a ``combo_id`` column in all-combinations mode

Each row carries a ``remark`` column with a one-sentence scientific
explanation of what the score / grade means, so a researcher who
only looks at the table can read it without consulting the
source code.

Plus ``inputs.tsv`` (one-row config snapshot), and the reproducibility
artefacts ``resolved_thresholds.yaml``, ``resolved_weights.json``,
``grade_legend.md``.

The path is derived from ``args.path`` at runtime (not hardcoded).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .combos import (
    build_per_metric_remark,
    build_verdict_remark,
)
from .composite import compute_all, load_weights, write_resolved_weights
from .diagnostics import ProblemVerdict
from .grades import grade_legend, letter
from .rules import write_resolved_thresholds, load_thresholds

logger = logging.getLogger("checkatlas")


def _safe_str(x) -> str:
    try:
        if x is None:
            return ""
        if isinstance(x, float) and np.isnan(x):
            return ""
        return str(x)
    except Exception:
        return ""


def _evidence_str(evidence: Iterable[str]) -> str:
    return "; ".join(str(e) for e in evidence)


def write_verdicts(
    verdicts: list[ProblemVerdict],
    atlas_name: str,
    outdir: str,
) -> Path:
    """Write ``verdicts.tsv`` — one row per (atlas, problem)."""
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "verdicts.tsv")
    rows: list[dict] = []
    for v in verdicts:
        score = float(v.score) if v.score is not None and not (
            isinstance(v.score, float) and np.isnan(v.score)
        ) else float("nan")
        rows.append(
            {
                "atlas": atlas_name,
                "problem_id": int(v.problem_id),
                "problem_name": v.problem_name,
                "score": round(score, 4) if score == score else "n/a",
                "grade": v.grade or letter(score),
                "confidence": v.confidence,
                "evidence": _evidence_str(v.evidence),
                "explanation": v.explanation,
                "reference": v.reference,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_composite(
    atlas_name: str,
    verdicts: list[ProblemVerdict],
    weights: dict,
    outdir: str,
    *,
    ilisi_penalty_norm: float = 0.5,
    failure_rate: float = 0.0,
    runtime_factor: float = 1.0,
) -> Path:
    """Write ``composite.tsv`` — one row per atlas."""
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "composite.tsv")
    scores = compute_all(
        verdicts,
        weights=weights,
        ilisi_penalty_norm=ilisi_penalty_norm,
        failure_rate=failure_rate,
        runtime_factor=runtime_factor,
    )
    by_id = {v.problem_id: v for v in verdicts}
    row = {
        "atlas": atlas_name,
        "fmf": round(scores["fmf"], 2),
        "bf": round(scores["bf"], 4),
        "pr": round(scores["pr"], 4),
    }
    for pid in range(1, 10):
        v = by_id.get(pid)
        if v is None:
            row[f"problem_{pid}_score"] = "n/a"
            row[f"problem_{pid}_grade"] = "n/a"
        else:
            score = float(v.score) if v.score is not None and not (
                isinstance(v.score, float) and np.isnan(v.score)
            ) else float("nan")
            row[f"problem_{pid}_score"] = (
                round(score, 4) if score == score else "n/a"
            )
            row[f"problem_{pid}_grade"] = v.grade or letter(score)
    df = pd.DataFrame([row])
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_per_metric(
    metrics_df: pd.DataFrame,
    outdir: str,
    atlas_name: str = "",
) -> Path:
    """Write ``per_metric.tsv`` — the long-format raw metrics.

    The output is a copy of the input long-format table with the
    atlas name added as a column. Missing column ``Atlas Name`` is
    filled with the supplied ``atlas_name`` so MultiQC can pivot on
    it.
    """
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "per_metric.tsv")
    df = metrics_df.copy()
    if "Atlas Name" not in df.columns:
        df.insert(0, "Atlas Name", atlas_name)
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_inputs_section(
    outdir: str,
    config_dict: dict,
) -> Path:
    """Write a small ``inputs.tsv`` (one row) summarising the
    inputs that were used. This is the Becavin comment-6 input-data
    sanity section.
    """
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "inputs.tsv")
    df = pd.DataFrame([config_dict])
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_grade_legend(outdir: str) -> Path:
    """Write a ``grade_legend.md`` to the output folder."""
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "grade_legend.md")
    with open(out, "w") as fh:
        fh.write(grade_legend())
    return Path(out)


def write_all(
    *,
    atlas_name: str,
    verdicts: list[ProblemVerdict],
    metrics_df: pd.DataFrame,
    outdir: str,
    weights_path: str | None = None,
    thresholds_path: str | None = None,
    ilisi_penalty_norm: float = 0.5,
    failure_rate: float = 0.0,
    runtime_factor: float = 1.0,
    config_dict: dict | None = None,
) -> dict[str, str]:
    """Run every writer. Returns a dict ``{filename: path}``.

    Output filenames (no ``scfm_`` prefix, the atlas name is in the
    path):

        verdicts.tsv
        composite.tsv
        per_metric.tsv
        inputs.tsv
        resolved_weights.json
        resolved_thresholds.yaml
        grade_legend.md
    """
    weights = load_weights(weights_path)
    rules = load_thresholds(thresholds_path)
    paths: dict[str, str] = {}
    paths["verdicts.tsv"] = str(
        write_verdicts(verdicts, atlas_name, outdir)
    )
    paths["composite.tsv"] = str(
        write_composite(
            atlas_name,
            verdicts,
            weights,
            outdir,
            ilisi_penalty_norm=ilisi_penalty_norm,
            failure_rate=failure_rate,
            runtime_factor=runtime_factor,
        )
    )
    paths["per_metric.tsv"] = str(
        write_per_metric(metrics_df, outdir, atlas_name=atlas_name)
    )
    paths["resolved_weights.json"] = str(
        write_resolved_weights(
            weights, os.path.join(outdir, "resolved_weights.json")
        )
    )
    paths["resolved_thresholds.yaml"] = str(
        write_resolved_thresholds(
            rules, os.path.join(outdir, "resolved_thresholds.yaml")
        )
    )
    if config_dict is not None:
        paths["inputs.tsv"] = str(
            write_inputs_section(outdir, config_dict)
        )
    paths["grade_legend.md"] = str(write_grade_legend(outdir))
    return paths


# ──────────────────────────────────────────────────────────────────
# All-combinations writers (per-row "remark" column)
# ──────────────────────────────────────────────────────────────────


def write_verdicts_with_combos(
    atlas_name: str,
    headline_verdicts: list[ProblemVerdict],
    all_verdicts: list[ProblemVerdict],
    outdir: str,
) -> Path:
    """Write ``verdicts.tsv`` with a ``combo_id`` column and a
    ``remark`` column.

    The first 9 rows are the *headline* (worst grade across
    combos per problem). The following rows are the per-combo
    verdicts in ``combo_id`` order. The ``remark`` column carries
    a one-sentence scientific explanation per row, so a
    researcher can interpret the table without reading the
    source code.
    """
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "verdicts.tsv")
    rows: list[dict] = []
    # Headline rows first
    for v in headline_verdicts:
        rows.append(
            {
                "atlas": atlas_name,
                "combo_id": "headline",
                "problem_id": int(v.problem_id),
                "problem_name": v.problem_name,
                "score": "n/a" if (v.score != v.score) else round(float(v.score), 4),
                "grade": v.grade or "n/a",
                "confidence": v.confidence,
                "evidence": v.explanation or "",
                "explanation": v.explanation,
                "reference": v.reference,
                "remark": v.explanation or "",
            }
        )
    # Per-combo rows
    for v in all_verdicts:
        combo_id = v.combo_id or "default"
        remark = build_verdict_remark(
            problem_id=int(v.problem_id),
            grade=v.grade,
            score=float(v.score) if v.score == v.score else float("nan"),
            combo=v.combo,
            explanation=v.explanation,
        )
        rows.append(
            {
                "atlas": atlas_name,
                "combo_id": combo_id,
                "problem_id": int(v.problem_id),
                "problem_name": v.problem_name,
                "score": (
                    "n/a" if (v.score != v.score) else round(float(v.score), 4)
                ),
                "grade": v.grade or "n/a",
                "confidence": v.confidence,
                "evidence": _evidence_str(v.evidence),
                "explanation": v.explanation,
                "reference": v.reference,
                "remark": remark,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_composite_with_combos(
    atlas_name: str,
    headline_scores: dict,
    per_combo_scores: list[dict],
    outdir: str,
) -> Path:
    """Write ``composite.tsv`` with one row per combo plus a
    ``headline`` row at the top.

    Each row has a ``remark`` column with a one-sentence
    scientific explanation of which (ref, pred, batch, scfm,
    baseline) was evaluated and how to read the scores.
    """
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "composite.tsv")
    rows: list[dict] = []
    # Headline first
    headline_row = {
        "atlas": atlas_name,
        "combo_id": "headline",
        "fmf": headline_scores.get("fmf", float("nan")),
        "bf": headline_scores.get("bf", float("nan")),
        "pr": headline_scores.get("pr", float("nan")),
        "n_worst_problems": headline_scores.get("n_worst_problems", 0),
        "n_combos": headline_scores.get("n_combos", 0),
        "remark": headline_scores.get("remark", ""),
    }
    for pid in range(1, 10):
        headline_row[f"problem_{pid}_grade"] = ""
        headline_row[f"problem_{pid}_score"] = ""
    rows.append(headline_row)
    # Per-combo rows
    for cs in per_combo_scores:
        row = {
            "atlas": atlas_name,
            "combo_id": cs.get("combo_id", "default"),
            "fmf": cs.get("fmf", float("nan")),
            "bf": cs.get("bf", float("nan")),
            "pr": cs.get("pr", float("nan")),
            "n_worst_problems": cs.get("n_worst_problems", 0),
            "n_combos": 1,
            "remark": cs.get("remark", ""),
        }
        for pid in range(1, 10):
            row[f"problem_{pid}_grade"] = ""
            row[f"problem_{pid}_score"] = ""
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_per_metric_with_combos(
    metrics_df: pd.DataFrame,
    outdir: str,
    atlas_name: str = "",
) -> Path:
    """Write ``per_metric.tsv`` with a ``combo_id`` column and a
    ``remark`` column.

    Each row carries a per-metric remark explaining what the
    metric measures and whether higher or lower is better, so
    a researcher can read the long-format table without
    consulting the scfm source code.
    """
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "per_metric.tsv")
    df = metrics_df.copy()
    if "Atlas Name" not in df.columns:
        df.insert(0, "Atlas Name", atlas_name)
    if "combo_id" not in df.columns:
        df["combo_id"] = "default"
    # Build the remark per row
    remarks: list[str] = []
    for _, r in df.iterrows():
        remarks.append(
            build_per_metric_remark(
                metric_name=str(r.get("Metric Name", "")),
                embedding=str(r.get("Embedding", "")),
                value=(
                    float(r["Value"])
                    if pd.notna(r.get("Value"))
                    else float("nan")
                ),
            )
        )
    df["remark"] = remarks
    df.to_csv(out, sep="\t", index=False)
    return Path(out)


def write_all_with_combos(
    *,
    atlas_name: str,
    headline_verdicts: list[ProblemVerdict],
    all_verdicts: list[ProblemVerdict],
    per_combo_scores: list[dict],
    headline_scores: dict,
    metrics_df: pd.DataFrame,
    outdir: str,
    weights_path: str | None = None,
    thresholds_path: str | None = None,
    config_dict: dict | None = None,
) -> dict[str, str]:
    """Run every combo-aware writer. Returns a dict
    ``{filename: path}``.

    Output filenames (no ``scfm_`` prefix, the atlas name is in
    the path):

        verdicts.tsv
        composite.tsv
        per_metric.tsv
        inputs.tsv
        resolved_weights.json
        resolved_thresholds.yaml
        grade_legend.md
    """
    from .rules import write_resolved_thresholds

    weights = load_weights(weights_path)
    rules = load_thresholds(thresholds_path)
    paths: dict[str, str] = {}
    paths["verdicts.tsv"] = str(
        write_verdicts_with_combos(
            atlas_name, headline_verdicts, all_verdicts, outdir
        )
    )
    paths["composite.tsv"] = str(
        write_composite_with_combos(
            atlas_name, headline_scores, per_combo_scores, outdir
        )
    )
    paths["per_metric.tsv"] = str(
        write_per_metric_with_combos(
            metrics_df, outdir, atlas_name=atlas_name
        )
    )
    paths["resolved_weights.json"] = str(
        write_resolved_weights(
            weights, os.path.join(outdir, "resolved_weights.json")
        )
    )
    paths["resolved_thresholds.yaml"] = str(
        write_resolved_thresholds(
            rules, os.path.join(outdir, "resolved_thresholds.yaml")
        )
    )
    if config_dict is not None:
        paths["inputs.tsv"] = str(
            write_inputs_section(outdir, config_dict)
        )
    paths["grade_legend.md"] = str(write_grade_legend(outdir))
    return paths
