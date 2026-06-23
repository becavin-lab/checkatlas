"""Reporting: write the scFM QC outputs.

Three TSV files, all idempotent, written flat into the per-atlas
``checkatlas_files/scfm/<atlas_name>/`` folder (so the layout matches
the per-task ``checkatlas_files/<task>/<atlas>.tsv`` pattern, with
the atlas name in the path rather than the filename):

  * ``verdicts.tsv``        — one row per (atlas, problem)
  * ``composite.tsv``       — one row per atlas, with FMF, BF, PR
  * ``per_metric.tsv``      — long-format raw metric values

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
