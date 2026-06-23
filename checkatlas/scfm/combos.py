"""All-combinations orchestration for the scfm layer.

When ``checkatlas scfm <path> --atlas_name X`` is called with **no**
``--scfm_*`` flags, the orchestrator runs all 3 per-task engines and
the diagnostic engine evaluates the (ref, pred, batch, scfm,
baseline) combinations the column detector found. This module is
the single source of truth for that combinatorial expansion.

Public surface
--------------

* :class:`ComboSpec` — one frozen (ref, pred, batch, scfm, baseline)
  tuple, with a deterministic ``combo_id`` for output tagging.
* :func:`detect_all_combos` — build the list of
  :class:`ComboSpec` instances from the column detector result,
  capped to the given maxima.
* :func:`expand_config_for_combo` — return a new ``SCFMConfig``
  with ref/pred/batch/scfm/baseline fields set from a combo.
* :func:`summarise_composite_across_combos` — compute the
  "headline" row of the composite table (mean across combos + the
  worst-problem-grade count).
* :func:`build_combo_remark` — generate the user-facing ``remark``
  text for a single combo's headline row, so the user can see at
  a glance *what* was evaluated and *what it means*.
* :func:`build_verdict_remark` — generate the per-row remark for a
  single (combo, problem) verdict.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Optional

from anndata import AnnData

from .config import SCFMConfig

logger = logging.getLogger("checkatlas")

# Regex-substring matchers for "scfm-shaped" embedding names. The
# diagnostic engine (problem 2 / 9) compares an scfm embedding to a
# baseline, so the scfm role needs a name that looks like a
# foundation-model embedding. If none exists, the orchestrator
# falls back to "the first non-UMAP embedding is the scfm, the rest
# are baselines" (a heuristic that works for Kedzierska-style
# atlases where X_geneformer is the only scfm).
_SCFM_NAME_PATTERNS = (
    "geneformer",
    "scgpt",
    "scgpt_human",
    "scgpt_blood",
    "scgpt_kidney",
    "uce",
    "scfoundation",
    "sccello",
    "scbert",
    "scfm",
    "foundation",
    "cellplm",
    "scprint",
    "scimilarity",
    "genept",
    "langcell",
    "cellfm",
    "transcriptformer",
)

# Maps a problem_id to a short, human-readable "what is this score
# measuring" sentence. Used as the per-verdict remark text.
_PROBLEM_DESCRIPTIONS = {
    1: "Zero-shot cell-type classification accuracy (ref-vs-pred)",
    2: "scFM-vs-baseline embedding separation quality",
    3: "How early the scFM's performance saturates as data grows",
    4: "Whether the scFM removes technical batch effects",
    5: "Information loss in the scFM's tokenization / representation",
    6: "Cross-domain generalisation (out-of-distribution cell types)",
    7: "Clinical / translational utility on patient cohorts",
    8: "Stability of the scFM across random seeds / subsamples",
    9: "Whether the scFM is just denoising (matches PCA+noise)",
}

_GRADE_LEGEND = {
    "A": "excellent — meets or exceeds the gold-standard baseline",
    "B": "good — comparable to the best non-scFM baseline",
    "C": "fair — passes the threshold but with caveats",
    "D": "poor — fails the threshold on this problem",
    "F": "failing — the scFM performs worse than the simplest baseline",
}


@dataclass(frozen=True)
class ComboSpec:
    """One (ref, pred, batch, scfm, baseline) combination.

    The five fields correspond to the roles the scfm diagnostic
    engine needs. Any of them can be ``None`` if the role is not
    used by the current problem; the engine handles n/a
    gracefully.
    """

    ref_label: Optional[str] = None
    predicted_label: Optional[str] = None
    batch_key: Optional[str] = None
    scfm_embedding: Optional[str] = None
    baseline_embedding: Optional[str] = None
    combo_id: str = field(default="")

    def __post_init__(self) -> None:
        if not self.combo_id:
            object.__setattr__(
                self,
                "combo_id",
                self._build_combo_id(),
            )

    def _build_combo_id(self) -> str:
        parts = [
            self.ref_label or "*",
            self.predicted_label or "*",
            self.batch_key or "*",
            self.scfm_embedding or "*",
            self.baseline_embedding or "*",
        ]
        return "|".join(parts)

    @classmethod
    def from_args(
        cls,
        ref_label: str,
        predicted_label: str,
        batch_key: str,
        scfm_embedding: str,
        baseline_embedding: str,
    ) -> "ComboSpec":
        return cls(
            ref_label=ref_label or None,
            predicted_label=predicted_label or None,
            batch_key=batch_key or None,
            scfm_embedding=scfm_embedding or None,
            baseline_embedding=baseline_embedding or None,
        )

    def roles_used(self) -> list[str]:
        """Return the list of role names that are filled in for
        this combo. Used by the diagnostic engine to know which
        problems can be evaluated."""
        out = []
        if self.ref_label:
            out.append("ref_label")
        if self.predicted_label:
            out.append("predicted_label")
        if self.batch_key:
            out.append("batch_key")
        if self.scfm_embedding:
            out.append("scfm_embedding")
        if self.baseline_embedding:
            out.append("baseline_embedding")
        return out


# ──────────────────────────────────────────────────────────────────
# Detection
# ──────────────────────────────────────────────────────────────────


def _looks_scfm(name: str) -> bool:
    n = name.lower()
    return any(p in n for p in _SCFM_NAME_PATTERNS)


def _classify_embeddings(
    embedding_keys: list[str],
) -> tuple[list[str], list[str]]:
    """Split a list of obsm keys into (scfm-shaped, baseline).

    If any key matches an scfm name pattern, those become the
    scfm candidates and the rest become baselines. If none
    matches, the first key is the scfm (best-effort fallback) and
    the rest are baselines. The caller should also handle the
    empty case (no embeddings at all).
    """
    scfm = [k for k in embedding_keys if _looks_scfm(k)]
    rest = [k for k in embedding_keys if k not in scfm]
    if not scfm and embedding_keys:
        # Fallback: assume the first embedding is the scfm. This
        # is what the Kedzierska-style pancreas h5ad wants when
        # the user adds synthetic X_geneformer with no name hint.
        scfm = [embedding_keys[0]]
        rest = embedding_keys[1:]
    return scfm, rest


def detect_all_combos(
    adata: AnnData,
    detected: dict[str, Any],
    *,
    max_ref: int = 3,
    max_pred: int = 2,
    max_batch: int = 2,
    max_scfm: int = 3,
    max_baseline: int = 3,
) -> list[ComboSpec]:
    """Build the list of (ref, pred, batch, scfm, baseline) combos
    from the column detector result.

    The function is *pure* — it depends only on the AnnData and
    the detected dict, not on the SCFMConfig. Caps are applied
    per role to keep the cartesian product bounded.
    """
    ref_candidates = [k for k, _ in detected["annotation"]["reference"][:max_ref]]
    pred_candidates = [k for k, _ in detected["annotation"]["predicted"][:max_pred]]
    batch_candidates = [k for k, _ in detected["batch"][:max_batch]]
    embedding_keys = [k for k, _ in detected["clustering"]["embeddings"]]
    scfm_keys, baseline_keys = _classify_embeddings(embedding_keys)
    scfm_keys = scfm_keys[:max_scfm]
    baseline_keys = baseline_keys[:max_baseline]

    if not scfm_keys:
        # No embeddings at all — the scfm step will produce n/a
        # verdicts. We still return one combo so the user gets a
        # single empty row in the composite and can see why.
        return [ComboSpec()]

    # If no baseline is detectable, every combo uses the same
    # embedding for scfm and baseline roles (problem 2 / 9
    # become n/a because there is nothing to compare against).
    if not baseline_keys:
        baseline_keys = [scfm_keys[0]]

    combos: list[ComboSpec] = []
    seen: set[str] = set()
    for scfm_emb in scfm_keys:
        for base_emb in baseline_keys:
            for ref in ref_candidates or [None]:
                for pred in pred_candidates or [None]:
                    for batch in batch_candidates or [None]:
                        combo = ComboSpec.from_args(
                            ref_label=ref or "",
                            predicted_label=pred or "",
                            batch_key=batch or "",
                            scfm_embedding=scfm_emb,
                            baseline_embedding=base_emb,
                        )
                        if combo.combo_id in seen:
                            continue
                        seen.add(combo.combo_id)
                        combos.append(combo)
    return combos


# ──────────────────────────────────────────────────────────────────
# Config expansion
# ──────────────────────────────────────────────────────────────────


def expand_config_for_combo(
    scfm_config: SCFMConfig, combo: ComboSpec
) -> SCFMConfig:
    """Return a new :class:`SCFMConfig` with ref/pred/batch/scfm/
    baseline fields set from ``combo``. All other fields
    preserved.
    """
    return SCFMConfig(
        atlas_name=scfm_config.atlas_name,
        scfm_embedding=combo.scfm_embedding or scfm_config.scfm_embedding,
        baseline_embeddings=(
            (combo.baseline_embedding,)
            if combo.baseline_embedding
            else scfm_config.baseline_embeddings
        ),
        predicted_label=combo.predicted_label or scfm_config.predicted_label,
        ref_label=combo.ref_label or scfm_config.ref_label,
        batch_key=combo.batch_key or scfm_config.batch_key,
        domain_key=scfm_config.domain_key,
        patient_key=scfm_config.patient_key,
        outcome_key=scfm_config.outcome_key,
        scaling_fractions=scfm_config.scaling_fractions,
        n_seeds=scfm_config.n_seeds,
        noise_sigma=scfm_config.noise_sigma,
        min_domain_size=scfm_config.min_domain_size,
        weights_path=scfm_config.weights_path,
        thresholds_path=scfm_config.thresholds_path,
    )


# ──────────────────────────────────────────────────────────────────
# Composite summary
# ──────────────────────────────────────────────────────────────────


def summarise_composite_across_combos(
    per_combo_scores: list[dict],
) -> dict[str, Any]:
    """Compute the "headline" row of the composite table.

    ``per_combo_scores`` is a list of dicts each shaped like the
    output of :func:`checkatlas.scfm.composite.compute_all`, e.g.
    ``{"fmf": 0.5, "bf": 1.0, "pr": 0.3, "combo_id": "..."}``.

    Returns a single dict with keys ``fmf``, ``bf``, ``pr``,
    ``n_combos``, ``n_worst_problems`` (count of problems with
    grade in {D, F} summed across combos), and ``remark``
    (a one-sentence human-readable summary).
    """
    if not per_combo_scores:
        return {
            "fmf": float("nan"),
            "bf": float("nan"),
            "pr": float("nan"),
            "n_combos": 0,
            "n_worst_problems": 0,
            "remark": (
                "No combinations were evaluated. The atlas may "
                "lack the columns or embeddings the scfm layer "
                "needs; check the orchestrator warnings for "
                "details."
            ),
        }
    fmfs = [s["fmf"] for s in per_combo_scores]
    bfs = [s["bf"] for s in per_combo_scores]
    prs = [s["pr"] for s in per_combo_scores]
    return {
        "fmf": round(sum(fmfs) / len(fmfs), 2),
        "bf": round(sum(bfs) / len(bfs), 4),
        "pr": round(sum(prs) / len(prs), 4),
        "n_combos": len(per_combo_scores),
        "n_worst_problems": sum(
            s.get("n_worst_problems", 0) for s in per_combo_scores
        ),
        "remark": build_headline_remark(per_combo_scores),
    }


def build_headline_remark(per_combo_scores: list[dict]) -> str:
    """Build a one-sentence remark for the headline composite row."""
    if not per_combo_scores:
        return "No scfm combinations could be evaluated."
    n = len(per_combo_scores)
    fmfs = [s["fmf"] for s in per_combo_scores]
    fmf_min, fmf_max = min(fmfs), max(fmfs)
    if fmf_max - fmf_min < 0.05:
        spread = "consistent across combinations"
    else:
        spread = (
            f"FMF spreads from {fmf_min:.2f} to {fmf_max:.2f} across "
            f"combinations — the result depends on which "
            f"(ref, pred, batch) the user picks"
        )
    return (
        f"Mean of {n} (ref, pred, batch, scfm, baseline) "
        f"combinations; {spread}. The per-combo breakdown is "
        f"in the rows below."
    )


# ──────────────────────────────────────────────────────────────────
# User-facing remark text
# ──────────────────────────────────────────────────────────────────


def build_combo_remark(combo: ComboSpec) -> str:
    """Build a user-facing remark for one combo's headline row.

    The remark is meant to be read by a researcher who is
    scanning the composite table: it tells them *which* columns
    and embeddings were evaluated, and which problems this
    particular combo can answer.
    """
    roles = combo.roles_used()
    if not roles:
        return (
            "No ref / pred / batch / embedding keys were "
            "detected; this combo produced n/a verdicts. See "
            "the orchestrator warnings for why."
        )
    parts = [f"{r}={getattr(combo, r)}" for r in roles]
    return (
        "Evaluated: " + ", ".join(parts) + ". "
        "This combo is one of N evaluated for the same atlas; "
        "compare across rows to see how sensitive the scfm "
        "verdict is to the column / embedding choice."
    )


def build_verdict_remark(
    problem_id: int,
    grade: str,
    score: float,
    combo: ComboSpec,
    explanation: str = "",
) -> str:
    """Build a one-sentence remark for a single (combo, problem)
    verdict row.

    The remark is meant to be read by a researcher who is
    scanning the verdicts table: it tells them *what* the score
    measures and how to interpret the grade.
    """
    what = _PROBLEM_DESCRIPTIONS.get(
        problem_id, "scfm diagnostic metric"
    )
    if grade in {"A", "B"}:
        quality = "the scFM is competitive on this problem"
    elif grade == "C":
        quality = "the scFM is borderline on this problem"
    elif grade in {"D", "F"}:
        quality = "the scFM fails on this problem"
    else:
        quality = "this problem could not be evaluated for this combo"
    grade_text = _GRADE_LEGEND.get(grade, "unknown")
    role_hint = ""
    if problem_id == 2 and combo.scfm_embedding and combo.baseline_embedding:
        role_hint = (
            f" (comparing {combo.scfm_embedding} against "
            f"{combo.baseline_embedding})"
        )
    elif problem_id == 4 and combo.batch_key:
        role_hint = f" (batch_key={combo.batch_key})"
    elif problem_id == 1 and combo.ref_label and combo.predicted_label:
        role_hint = (
            f" (ref={combo.ref_label}, pred={combo.predicted_label})"
        )
    if grade == "n/a" or grade == "":
        return (
            f"Problem {problem_id}: {what}. This combo could "
            f"not evaluate it{role_hint}; the reason is: "
            f"{explanation or 'no relevant metric in the data'}. "
            f"See orchestrator warnings for details."
        )
    return (
        f"Problem {problem_id}: {what}{role_hint}. Score "
        f"{score:.2f} (grade {grade} — {grade_text}). "
        f"{quality}."
    )


def build_per_metric_remark(
    metric_name: str, embedding: str, value: float
) -> str:
    """Build a one-line remark for a per-metric long-format row.

    The remark explains *what* the metric measures, so a
    researcher who only looks at the long table can read it
    without reading the scfm source code.
    """
    what = {
        "silhouette": "cluster separation (1=well-separated, 0=overlapping)",
        "davies_bouldin": "cluster separation (lower is better)",
        "calinski_harabasz": "cluster separation (higher is better)",
        "isolated_f1_score": "ref-vs-pred F1 (1=perfect, 0=mismatch)",
        "adj_rand_index": "ref-vs-pred ARI (1=perfect, 0=random)",
        "average_silhouette_width": "silhouette (variant used by scIB)",
        "kbet": "batch-effect rejection rate (lower=more mixed)",
        "ilisi": "iLISI integration score (1=perfect batch mixing)",
        "graph_connectivity": "graph connectivity (1=fully connected)",
        "pcr": "PCR batch score (1=perfect mixing)",
        "dCor": "distance correlation to raw X",
        "kruskal_stress": "Kruskal stress (0=perfect embedding)",
        "spearman_rho": "Spearman correlation to raw X ranks",
        "coknn": "co-kNN preservation",
        "entourage": "entourage score (local neighborhood preservation)",
        "den_pre": "denoising preservation",
        "lcmc": "local continuity meta-metric",
        "avg_jaccard_dis": "average Jaccard distance (1=perfect)",
        "ged": "graph-based embedding distance",
        "stability_cv": "coefficient of variation across subsamples (lower=more stable)",
        "rare_f1": "F1 on rare cell types",
        "common_f1": "F1 on common cell types",
        "rare_common_gap": "F1 gap between rare and common cell types",
        "cross_domain_ari": "cross-domain ARI",
    }.get(metric_name, "scfm diagnostic metric")
    direction = "higher is better" if metric_name in {
        "silhouette", "calinski_harabasz", "isolated_f1_score",
        "adj_rand_index", "average_silhouette_width", "ilisi",
        "graph_connectivity", "dCor", "spearman_rho", "coknn",
        "entourage", "den_pre", "lcmc", "avg_jaccard_dis", "ged",
        "common_f1",
    } else "lower is better"
    return (
        f"{metric_name} on {embedding}: {value:.4f}. "
        f"({what}; {direction}.)"
    )
