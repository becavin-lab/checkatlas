"""Literature-grounded threshold table for the nine scFM problems.

Each entry maps a problem ID to:

  * a human-readable name
  * a list of ``(metric_name, operator, threshold, citation)`` rules
  * the primary literature reference
  * a short explanation string

Becavin comment 3: thresholds are user-overridable via a YAML file
(``--scfm_thresholds PATH``). The hard-coded ``DEFAULT_RULES`` below
is the fallback. The resolved thresholds (after merge) are written
to ``checkatlas_files/scfm/resolved_thresholds.yaml`` for
reproducibility.

Threshold values are taken from the papers in
``docs/scfm_benchmarking_report.md`` (Section "Summary Table") and
the discussions in §2 of this implementation plan. They are
*starting points*, not gospel — users can override them with their
own YAML.
"""

from __future__ import annotations

import logging
import os
from typing import Any

logger = logging.getLogger("checkatlas")


DEFAULT_RULES: dict[int, dict[str, Any]] = {
    1: {
        "name": "Zero-shot performance failure",
        "metric_thresholds": [
            ("isolated_f1", "<", 0.50, "Kedzierska 2024"),
            ("adj_rand_index", "<", 0.30, "Kedzierska 2025"),
        ],
        "primary_reference": "Kedzierska 2024, 2025; Boiarsky 2023",
        "explanation": (
            "Zero-shot predictions on the scFM embedding are below "
            "the F1 / ARI bands reported for the failing scFMs in "
            "Kedzierska 2024."
        ),
    },
    2: {
        "name": "Baselines outperform scFM",
        "metric_thresholds": [
            ("silhouette_delta_vs_pca", "<", -0.01, "Liu 2024"),
            ("isolated_f1_delta_vs_pca", "<", -0.01, "Souza & Mehta 2026"),
        ],
        "primary_reference": "Liu 2024; Souza & Mehta 2026; DenAdel 2025",
        "explanation": (
            "For every metric, the scFM embedding scores below the "
            "PCA baseline by more than the margin reported in Liu "
            "2024 and Souza & Mehta 2026."
        ),
    },
    3: {
        "name": "Pretraining scaling laws do not apply",
        "metric_thresholds": [
            ("plateau_ratio", "<", 1.10, "DenAdel 2025"),
        ],
        "primary_reference": "DenAdel 2025; Wang 2025",
        "explanation": (
            "Performance plateaus before 10% of the data, "
            "reproducing the DenAdel 2025 scaling finding."
        ),
    },
    4: {
        "name": "Batch effects persist",
        "metric_thresholds": [
            ("kbet", ">", 0.50, "Wang, Zhang & Zhang 2025"),
            ("ilisi", "<", 0.30, "Wang, Zhang & Zhang 2025"),
            ("pcr", ">", 0.85, "Wang, Zhang & Zhang 2025"),
            ("graph_connectivity", "<", 0.50, "Wang, Zhang & Zhang 2025"),
        ],
        "primary_reference": "Wang, Zhang & Zhang 2025",
        "explanation": (
            "Batch effects dominate the embedding: kBET rejection "
            "is high and iLISI is low, matching the 8 scFMs "
            "evaluated by Wang, Zhang & Zhang 2025."
        ),
    },
    5: {
        "name": "Tokenization and representation failures",
        "metric_thresholds": [
            ("rare_common_gap", ">", 0.15, "Atti & Subramaniam 2025"),
            ("dCor_vs_pca", "<", 0.50, "Atti & Subramaniam 2025"),
        ],
        "primary_reference": "Atti & Subramaniam 2025; DenAdel 2025; Haber 2025",
        "explanation": (
            "Tokenization / representation causes information loss: "
            "rare-type recall is much lower than common-type recall "
            "(scGPT CD14 finding) and dCor to raw X is below the "
            "Atti & Subramaniam 2025 band."
        ),
    },
    6: {
        "name": "Generalization failures across domains",
        "metric_thresholds": [
            ("cross_domain_ari_mean", "<", 0.30, "Souza & Mehta 2026"),
        ],
        "primary_reference": "Souza & Mehta 2026; Boiarsky 2023; Wenteler 2024",
        "explanation": (
            "For at least 50% of (train, test) domain pairs the "
            "ARI on the test domain is below the Souza & Mehta "
            "2026 band."
        ),
    },
    7: {
        "name": "Real-world clinical / translational failures",
        "metric_thresholds": [
            ("asw_patient", "<", 0.05, "Elmarakeby 2025"),
            ("outcome_AUC", "<", 0.60, "Elmarakeby 2025"),
        ],
        "primary_reference": "Elmarakeby 2025; Csendes 2025",
        "explanation": (
            "Patient-level separation and outcome prediction are "
            "below the Elmarakeby 2025 cancer benchmarks."
        ),
    },
    8: {
        "name": "Model stability and interpretability issues",
        "metric_thresholds": [
            ("stability_cv_mean", ">", 0.10, "Liu 2024"),
        ],
        "primary_reference": "Liu 2024; Xu 2026",
        "explanation": (
            "Coefficient of variation across random subsamples is "
            "above the Liu 2024 / Xu 2026 instability band on at "
            "least two metrics."
        ),
    },
    9: {
        "name": "Denoising hypothesis (linear manifold)",
        "metric_thresholds": [
            ("baseline_competitiveness", "<", 0.20, "Souza & Mehta 2026"),
        ],
        "primary_reference": "Souza & Mehta 2026",
        "explanation": (
            "scFM ≤ PCA+δ on at least 80% of metrics, matching the "
            "Souza & Mehta 2026 'denoiser only' finding."
        ),
    },
}


def _yaml_available() -> bool:
    try:
        import yaml  # noqa: F401

        return True
    except ImportError:
        return False


def load_thresholds(thresholds_path: str | None) -> dict[int, dict[str, Any]]:
    """Load threshold table from YAML or fall back to defaults.

    Parameters
    ----------
    thresholds_path : str, optional
        Path to a YAML file. The file may contain a partial table
        (only the problems the user wants to override). Missing
        problems fall back to ``DEFAULT_RULES``. Unknown keys warn
        and are ignored.
    """
    rules = {k: {**v} for k, v in DEFAULT_RULES.items()}

    if not thresholds_path:
        return rules
    if not os.path.exists(thresholds_path):
        logger.warning(
            "scfm.thresholds: file not found at %s, using defaults", thresholds_path
        )
        return rules
    if not _yaml_available():
        logger.warning(
            "scfm.thresholds: PyYAML not installed, using defaults"
        )
        return rules

    import yaml

    with open(thresholds_path) as fh:
        user = yaml.safe_load(fh) or {}

    for key, value in user.items():
        try:
            problem_id = int(key)
        except (TypeError, ValueError):
            logger.warning("scfm.thresholds: skipping non-int key %r", key)
            continue
        if problem_id not in rules:
            logger.warning(
                "scfm.thresholds: unknown problem %d, skipping", problem_id
            )
            continue
        if not isinstance(value, dict):
            logger.warning(
                "scfm.thresholds: problem %d value is not a dict, skipping",
                problem_id,
            )
            continue
        rules[problem_id].update(value)

    return rules


def write_resolved_thresholds(
    rules: dict[int, dict[str, Any]], out_path: str
) -> str:
    """Write the resolved (after merge) thresholds to YAML. Returns the path."""

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    if not _yaml_available():
        # Fall back to JSON-ish plain text
        with open(out_path, "w") as fh:
            for pid, r in sorted(rules.items()):
                fh.write(f"problem {pid}: {r['name']}\n")
                for metric, op, thr, ref in r.get("metric_thresholds", []):
                    fh.write(f"  - ({metric} {op} {thr}) [{ref}]\n")
        return out_path
    import yaml

    serialisable = {}
    for pid, r in sorted(rules.items()):
        serialisable[str(pid)] = {
            "name": r["name"],
            "primary_reference": r.get("primary_reference", ""),
            "explanation": r.get("explanation", ""),
            "metric_thresholds": [
                {"metric": m, "op": op, "threshold": t, "citation": ref}
                for m, op, t, ref in r.get("metric_thresholds", [])
            ],
        }
    with open(out_path, "w") as fh:
        yaml.safe_dump(serialisable, fh, sort_keys=False)
    return out_path
