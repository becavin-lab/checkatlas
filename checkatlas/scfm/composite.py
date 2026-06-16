"""Three in-house composite scores.

FMF — FoundationModelFitness (headline, 0..100). Geometric mean of
      per-problem "fitness" components so any one failure pulls the
      whole score down.

BF  — BiologicalFaithfulness (0..1). Net biology-preservation score:
      positive contribution from cLISI / graph_connectivity /
      cell_cycle / HVG overlap, penalty from iLISI (batch leakage).

PR  — ProductionReadiness (0..1). Sigmoid of stability,
      tokenization loss, runtime, and failure rate.

Weights are hard-coded by default. A user-supplied JSON file at
``--scfm_weights`` overrides individual values; missing keys fall
back to defaults. The resolved weights are written to
``checkatlas_files/scfm/resolved_weights.json`` for reproducibility.
"""

from __future__ import annotations

import json
import logging
import math
import os
from typing import Any, Optional

logger = logging.getLogger("checkatlas")


DEFAULT_WEIGHTS: dict[str, Any] = {
    "fmf": {
        "batch": 1.0,
        "rare_types": 1.0,
        "cross_domain": 1.0,
        "scaling": 1.0,
        "stability": 1.0,
        "baseline_competitiveness": 1.0,
        "tokenization": 1.0,
    },
    "bf": {
        "cLISI": 1.0,
        "graph_connectivity": 1.0,
        "cell_cycle": 0.5,
        "hvg_overlap": 0.5,
        "iLISI_penalty": 1.0,
    },
    "pr": {
        "stability": 1.0,
        "tokenization": 1.0,
        "runtime": 0.5,
        "failure_rate": 1.5,
    },
}


def load_weights(weights_path: Optional[str]) -> dict[str, Any]:
    """Merge user-supplied weights over defaults."""
    weights: dict[str, Any] = {
        section: {**defaults} for section, defaults in DEFAULT_WEIGHTS.items()
    }
    if not weights_path or not os.path.exists(weights_path):
        return weights
    try:
        with open(weights_path) as fh:
            user = json.load(fh)
    except Exception as exc:
        logger.warning("composite.load_weights: %s, using defaults", exc)
        return weights
    if not isinstance(user, dict):
        logger.warning("composite.load_weights: not a dict, using defaults")
        return weights
    for section, section_weights in user.items():
        if section not in weights or not isinstance(section_weights, dict):
            logger.warning(
                "composite.load_weights: unknown section %r, skipping", section
            )
            continue
        for k, v in section_weights.items():
            if k not in weights[section]:
                logger.warning(
                    "composite.load_weights: unknown key %r in %r, skipping",
                    k,
                    section,
                )
                continue
            try:
                weights[section][k] = float(v)
            except (TypeError, ValueError):
                logger.warning(
                    "composite.load_weights: non-numeric value for %s/%s",
                    section,
                    k,
                )
    return weights


def write_resolved_weights(
    weights: dict[str, Any], out_path: str
) -> str:
    """Write the resolved (after merge) weights to JSON. Returns the path."""
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(weights, fh, indent=2, sort_keys=True)
    return out_path


# ──────────────────────────────────────────────────────────────────
# FMF: FoundationModelFitness
# ──────────────────────────────────────────────────────────────────


def _safe_score(verdict_score: float) -> float:
    """Convert a verdict score (0..1, higher = more problem) to a
    fitness component (0..1, higher = better). NaN → 0.5 (neutral)."""
    if verdict_score is None or (
        isinstance(verdict_score, float) and math.isnan(verdict_score)
    ):
        return 0.5
    return float(max(0.0, min(1.0, 1.0 - float(verdict_score))))


def fmf(verdicts, weights: dict[str, Any]) -> float:
    """Compute FMF (0..100) from nine ProblemVerdict objects.

    FMF uses a *weighted* geometric mean of seven fitness components
    derived from the nine verdicts. The weight of each component is
    taken from ``weights["fmf"]`` and the missing entries default to
    1.0. A component with weight 0.0 is dropped from the mean. A
    component with value 0.0 (a hard failure) floors the FMF at 0
    so the score honestly reflects the failure.
    """
    by_id = {v.problem_id: v for v in verdicts}
    components: dict[str, float] = {
        "batch": _safe_score(by_id.get(4).score if 4 in by_id else float("nan")),
        "rare_types": _safe_score(
            by_id.get(5).score if 5 in by_id else float("nan")
        ),
        "cross_domain": _safe_score(
            by_id.get(6).score if 6 in by_id else float("nan")
        ),
        "scaling": _safe_score(
            by_id.get(3).score if 3 in by_id else float("nan")
        ),
        "stability": _safe_score(
            by_id.get(8).score if 8 in by_id else float("nan")
        ),
        "baseline_competitiveness": _safe_score(
            by_id.get(9).score if 9 in by_id else float("nan")
        ),
        "tokenization": _safe_score(
            by_id.get(5).score if 5 in by_id else float("nan")
        ),
    }
    w = weights.get("fmf", {})
    if any(v <= 0.0 for v in components.values()):
        return 0.0
    log_acc = 0.0
    total_w = 0.0
    for name, val in components.items():
        wi = float(w.get(name, 1.0))
        if wi <= 0.0:
            continue
        log_acc += wi * math.log(val)
        total_w += wi
    if total_w <= 0.0:
        return 0.0
    gmean = math.exp(log_acc / total_w)
    return float(gmean * 100.0)


# ──────────────────────────────────────────────────────────────────
# BF: BiologicalFaithfulness
# ──────────────────────────────────────────────────────────────────


def bf(
    *,
    clisi_norm: float = 0.5,
    graph_connectivity_norm: float = 0.5,
    cell_cycle_norm: float = 0.5,
    hvg_overlap_norm: float = 0.5,
    ilisi_penalty_norm: float = 0.5,
    weights: dict[str, Any] | None = None,
) -> float:
    """Net biology-preservation score in [0, 1].

    Higher is better. ``ilisi_penalty_norm`` is *subtracted* (it
    represents batch leakage) and the result is clipped to [0, 1].
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS
    w = weights.get("bf", {})
    pos = (
        float(w.get("cLISI", 1.0)) * float(clisi_norm)
        + float(w.get("graph_connectivity", 1.0)) * float(graph_connectivity_norm)
        + float(w.get("cell_cycle", 0.5)) * float(cell_cycle_norm)
        + float(w.get("hvg_overlap", 0.5)) * float(hvg_overlap_norm)
    )
    neg = float(w.get("iLISI_penalty", 1.0)) * float(ilisi_penalty_norm)
    raw = pos - neg
    return float(max(0.0, min(1.0, raw)))


# ──────────────────────────────────────────────────────────────────
# PR: ProductionReadiness
# ──────────────────────────────────────────────────────────────────


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def pr(
    *,
    stability_term: float = 0.5,
    tokenization_term: float = 0.5,
    runtime_factor: float = 1.0,
    failure_rate: float = 0.0,
    weights: dict[str, Any] | None = None,
) -> float:
    """Sigmoid-based production-readiness score in [0, 1]."""
    if weights is None:
        weights = DEFAULT_WEIGHTS
    w = weights.get("pr", {})
    x = (
        float(w.get("stability", 1.0)) * float(stability_term)
        + float(w.get("tokenization", 1.0)) * float(tokenization_term)
        + float(w.get("runtime", 0.5)) * math.log(max(runtime_factor, 1e-9))
        - float(w.get("failure_rate", 1.5)) * float(failure_rate)
    )
    return float(_sigmoid(x))


# ──────────────────────────────────────────────────────────────────
# Top-level: compute all three
# ──────────────────────────────────────────────────────────────────


def compute_all(
    verdicts,
    *,
    weights: dict[str, Any] | None = None,
    ilisi_penalty_norm: float = 0.5,
    failure_rate: float = 0.0,
    runtime_factor: float = 1.0,
) -> dict[str, float]:
    """Compute FMF, BF, PR and return a flat dict.

    The dict is suitable for direct DataFrame construction:

        pd.DataFrame([compute_all(verdicts)])
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS
    return {
        "fmf": fmf(verdicts, weights),
        "bf": bf(ilisi_penalty_norm=ilisi_penalty_norm, weights=weights),
        "pr": pr(
            failure_rate=failure_rate,
            runtime_factor=runtime_factor,
            weights=weights,
        ),
    }
