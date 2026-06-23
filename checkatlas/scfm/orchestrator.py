"""Pre-scfm orchestration: ensure per-task metric TSVs exist on disk.

The scfm diagnostic engine consumes the wide-format TSVs written by
``checkatlas metric_cluster / metric_annot / metric_dimred`` (see
:mod:`checkatlas.scfm.io`). On a first-time atlas those TSVs do not
exist yet. This module is the *single* place that decides whether to
auto-run the corresponding ``atlas.create_metric_*`` step before
``run_scfm_from_cache`` is called.

Contract
--------

* No recompute logic lives inside :mod:`checkatlas.scfm.pipeline` or
  :mod:`checkatlas.scfm.io` — those modules are pure load-then-diagnose.
* The orchestrator is invoked by the CLI / Nextflow layer
  (see ``checkatlas.__main__``) between the scfm config build and the
  ``run_scfm_from_cache`` call.
* The column detector is the source of truth for "what is in the
  atlas". When the user explicitly passed ``--scfm_ref_label`` /
  ``--scfm_predicted_label`` / ``--scfm_batch_key`` etc., those
  user-specified keys are locked in and used directly, regardless of
  what the column detector would have picked.
* When the orchestrator auto-runs a per-task engine it does NOT
  silence the engine's own per-metric progress / verbose output. The
  user wants to see what metrics are being computed at the backend.
* When the orchestrator skips a task because the atlas lacks the
  required keys, it emits a single scientific warning that names the
  missing keys. No hard error.
"""

from __future__ import annotations

import logging
import time
from typing import Any, Optional

from anndata import AnnData

from .. import atlas
from ..utils.col_detector import CheckAtlasColumnDetector
from . import io as scfm_io
from .config import SCFMConfig

logger = logging.getLogger("checkatlas")

_TASK_CLUSTER = "cluster"
_TASK_ANNOT = "annotation"
_TASK_DIMRED = "dimred"

_TASKS = (_TASK_CLUSTER, _TASK_ANNOT, _TASK_DIMRED)


# ──────────────────────────────────────────────────────────────────
# Pure helpers (no I/O, no AnnData)
# ──────────────────────────────────────────────────────────────────


def required_tasks(scfm_config: SCFMConfig) -> set[str]:
    """Return the set of per-task metric categories the scfm
    diagnostic engine actually consumes for the given config.

    Mapping
    -------
    * ``annotation`` — needed whenever the user names a ref / predicted
      / batch column, or whenever ``--scfm_batch_key`` is set (kBET,
      iLISI, PCR, graph connectivity all live in the annotation TSV).
    * ``cluster`` — needed whenever ``--scfm_ref_label`` is set
      (silhouette, ARI, dCor, rare-type gap all key off cluster
      labels in the cluster TSV).
    * ``dimred`` — needed whenever the user names an scFM embedding
      or any baseline embedding (dCor, Kruskal stress, Spearman rho
      are the cross-embedding comparisons problem 2 / 9 use).

    The function is pure — it depends only on the SCFMConfig fields,
    never on the AnnData. This makes it cheap to call early, before
    the AnnData is even loaded, when the CLI is still deciding
    whether the scfm invocation is well-formed.
    """
    needed: set[str] = set()
    if (
        scfm_config.ref_label
        or scfm_config.predicted_label
        or scfm_config.batch_key
    ):
        needed.add(_TASK_ANNOT)
    if scfm_config.ref_label:
        needed.add(_TASK_CLUSTER)
    if scfm_config.scfm_embedding or scfm_config.baseline_embeddings:
        needed.add(_TASK_DIMRED)
    return needed


def user_overrides(scfm_config: SCFMConfig) -> dict[str, Optional[str]]:
    """Return the user-specified keys as a plain dict, for inspection
    by the orchestrator. Values are ``None`` when the user did not
    pass the corresponding flag (meaning "let the column detector
    decide").

    The function does NOT mutate the config; it only reports which
    flags were set explicitly.
    """
    return {
        "ref_label": scfm_config.ref_label,
        "predicted_label": scfm_config.predicted_label,
        "batch_key": scfm_config.batch_key,
        "domain_key": scfm_config.domain_key,
        "patient_key": scfm_config.patient_key,
        "outcome_key": scfm_config.outcome_key,
    }


# ──────────────────────────────────────────────────────────────────
# Per-task gate checks
# ──────────────────────────────────────────────────────────────────


def _gate_cluster(
    adata: AnnData, detected: dict[str, Any], scfm_config: SCFMConfig
) -> tuple[bool, str]:
    """Return ``(ok, reason)``.

    ``ok`` is True iff the atlas has the keys required to run
    ``create_metric_cluster``. ``reason`` is a scientific explanation
    used in the orchestrator warning when ``ok`` is False.
    """
    has_cluster_label = any(
        col in adata.obs.columns
        for col, _ in detected["clustering"]["cluster_labels"]
    )
    if scfm_config.ref_label:
        if scfm_config.ref_label in adata.obs.columns:
            return True, ""
        return (
            False,
            f"user-specified --scfm_ref_label={scfm_config.ref_label!r} "
            f"is not present in adata.obs.columns",
        )
    if has_cluster_label:
        return True, ""
    return (
        False,
        "no cluster-label column detected in adata.obs "
        "(detector found no leiden / louvain / seurat_clusters / "
        "celltype / annotation-style column matching the cluster "
        "semantic patterns)",
    )


def _gate_annot(
    adata: AnnData, detected: dict[str, Any], scfm_config: SCFMConfig
) -> tuple[bool, str]:
    """Return ``(ok, reason)`` for the annotation task."""
    user_keys = [
        scfm_config.ref_label,
        scfm_config.predicted_label,
        scfm_config.batch_key,
    ]
    user_keys_present = [k for k in user_keys if k and k in adata.obs.columns]
    user_keys_missing = [k for k in user_keys if k and k not in adata.obs.columns]
    if user_keys_present:
        if user_keys_missing:
            return (
                False,
                f"user-specified keys are partially missing in adata.obs: "
                f"{user_keys_missing}",
            )
        return True, ""

    has_ref = bool(detected["annotation"]["reference"])
    has_pred = bool(detected["annotation"]["predicted"])
    has_batch = bool(detected["batch"])
    if has_ref or has_pred or has_batch:
        return True, ""
    return (
        False,
        "no reference, predicted, or batch column detected in adata.obs "
        "(detector found no celltype / annotation / leiden / donor / "
        "sample / batch-style column matching the annotation semantic "
        "patterns)",
    )


def _gate_dimred(
    adata: AnnData, detected: dict[str, Any], scfm_config: SCFMConfig
) -> tuple[bool, str]:
    """Return ``(ok, reason)`` for the dimred task."""
    user_embs = []
    if scfm_config.scfm_embedding:
        user_embs.append(scfm_config.scfm_embedding)
    user_embs.extend(scfm_config.baseline_embeddings)
    user_embs_missing = [e for e in user_embs if e not in adata.obsm]
    if user_embs:
        if user_embs_missing:
            return (
                False,
                f"user-specified embedding(s) are not present in adata.obsm: "
                f"{user_embs_missing}",
            )
        return True, ""

    detector_embs = [k for k, _ in detected["clustering"]["embeddings"]]
    if detector_embs:
        return True, ""
    return (
        False,
        "no embedding with > 3 components detected in adata.obsm "
        "(detector found no X_pca / X_scvi / X_umap-style key, and "
        "2-D UMAP / t-SNE visualisations are intentionally excluded "
        "because they cannot drive a kNN-based metric)",
    )


_GATE_FUNCS = {
    _TASK_CLUSTER: _gate_cluster,
    _TASK_ANNOT: _gate_annot,
    _TASK_DIMRED: _gate_dimred,
}


# ──────────────────────────────────────────────────────────────────
# Public orchestrator
# ──────────────────────────────────────────────────────────────────


def ensure_per_task_tsvs(
    adata: AnnData,
    atlas_info: dict,
    args,
    scfm_config: SCFMConfig,
) -> dict[str, Optional[Any]]:
    """Ensure the wide-format per-task TSVs needed by the scfm
    diagnostic engine exist on disk. If any are missing and the atlas
    has the required keys, auto-run the corresponding
    ``atlas.create_metric_*`` step. If the keys themselves are
    missing, emit a single scientific warning per missing task and
    move on.

    Parameters
    ----------
    adata : AnnData
        The preprocessed atlas.
    atlas_info : dict
        The same ``atlas_info`` dict that was used for
        ``preprocess_atlas``. Required so we can call
        ``atlas.create_metric_*`` with the right path / name.
    args : argparse.Namespace
        The full CLI args. Used for ``args.path``, ``args.n_jobs``,
        ``args.metric_cluster``, ``args.metric_annot``,
        ``args.metric_dimred`` (so the user can narrow the metric
        set even in the auto-run path).
    scfm_config : SCFMConfig

    Returns
    -------
    dict
        ``{"cluster": df|None, "annotation": df|None, "dimred": df|None}``,
        refreshed after the auto-runs. ``None`` for any task that
        could not be computed.
    """
    atlas_name = scfm_config.atlas_name
    base_path = getattr(args, "path", ".") if args is not None else "."

    # ── 0. Ensure the checkatlas folder tree exists ──────────────
    # The atlas.create_metric_* functions write their wide TSVs to
    # ``checkatlas_files/<task>/<atlas>.tsv`` and expect those
    # directories to exist. The CLI calls ``folders.checkatlas_folders``
    # once in main(), but notebook / programmatic callers (e.g. the
    # smoke test) may not. We materialise the tree here so the
    # auto-run path is self-contained.
    from ..utils import folders
    folders.checkatlas_folders(base_path)

    # ── 1. Column-detector summary (single INFO line) ──────────────
    detected = CheckAtlasColumnDetector(adata).detect_all_parameters()
    n_ref = len(detected["annotation"]["reference"])
    n_pred = len(detected["annotation"]["predicted"])
    n_clust = len(detected["clustering"]["cluster_labels"])
    n_emb = len(detected["clustering"]["embeddings"])
    n_batch = len(detected["batch"])
    logger.info(
        "scfm: column detector for %s found %d ref / %d predicted / "
        "%d cluster label / %d embedding / %d batch column(s)",
        atlas_name,
        n_ref,
        n_pred,
        n_clust,
        n_emb,
        n_batch,
    )

    # ── 2. Decide which tasks are needed ──────────────────────────
    needed = required_tasks(scfm_config)
    if not needed:
        logger.info(
            "scfm: no per-task inputs were requested via --scfm_* flags "
            "for %s; nothing to auto-run",
            atlas_name,
        )
        return scfm_io.load_per_task_tsvs(
            base_path, atlas_name, verbose=False
        )

    # ── 3. Inspect existing TSVs ──────────────────────────────────
    # We use verbose=False for the orchestrator's internal load to
    # avoid double-logging: the per-TSV "loaded N rows" line is
    # emitted by run_scfm_from_cache, not here.
    existing = scfm_io.load_per_task_tsvs(
        base_path, atlas_name, verbose=False
    )

    # ── 4. Per-task: skip / reuse / auto-run ──────────────────────
    _RUNNERS = {
        _TASK_CLUSTER: (
            atlas.create_metric_cluster,
            lambda: (
                "Running full clustering pipeline for %s" % atlas_name
            ),
        ),
        _TASK_ANNOT: (
            atlas.create_metric_annot,
            lambda: (
                "Running full annotation pipeline for %s" % atlas_name
            ),
        ),
        _TASK_DIMRED: (
            atlas.create_metric_dimred,
            lambda: (
                "Running full dimred pipeline for %s" % atlas_name
            ),
        ),
    }

    for task in _TASKS:
        if task not in needed:
            continue
        if existing.get(task) is not None:
            logger.info(
                "scfm: %s TSV present at %s — skipping auto-run",
                task,
                scfm_io._tsv_path(base_path, task, atlas_name),
            )
            continue
        ok, reason = _GATE_FUNCS[task](adata, detected, scfm_config)
        if not ok:
            logger.warning(
                "scfm: %s task cannot be computed for %s — %s. "
                "This task will be skipped; the scfm diagnostic will "
                "report n/a for the corresponding problems.",
                task,
                atlas_name,
                reason,
            )
            continue
        runner, _ = _RUNNERS[task]
        logger.info(
            "scfm: %s TSV missing at %s — auto-running "
            "checkatlas metric_%s for %s (column detector found %s)",
            task,
            scfm_io._tsv_path(base_path, task, atlas_name),
            task,
            atlas_name,
            _detected_summary_for(task, detected, adata, scfm_config),
        )
        t0 = time.time()
        try:
            # NOTE: we deliberately do NOT silence the engine's own
            # INFO / progress output. The user wants to see what
            # metrics are being computed at the backend.
            runner(adata, atlas_info, args)
        except Exception as exc:
            logger.warning(
                "scfm: metric_%s auto-run failed for %s: %s. "
                "Continuing with whatever TSVs are present on disk.",
                task,
                atlas_name,
                exc,
            )
            continue
        elapsed = time.time() - t0
        logger.info(
            "scfm: metric_%s complete for %s in %.1fs -> %s",
            task,
            atlas_name,
            elapsed,
            scfm_io._tsv_path(base_path, task, atlas_name),
        )

    # ── 5. Re-read so the caller gets the post-orchestration state ─
    # verbose=False: run_scfm_from_cache will print the "loaded N rows"
    # lines itself.
    return scfm_io.load_per_task_tsvs(
        base_path, atlas_name, verbose=False
    )


def _detected_summary_for(
    task: str,
    detected: dict[str, Any],
    adata: AnnData,
    scfm_config: SCFMConfig,
) -> str:
    """Return a short human-readable summary of what the column
    detector found for one task, for the orchestrator's auto-run
    announcement line.
    """
    if task == _TASK_CLUSTER:
        keys = [k for k, _ in detected["clustering"]["cluster_labels"]]
        return (
            f"{len(keys)} cluster label(s) in adata.obs"
            + (f": {keys[:3]}{'...' if len(keys) > 3 else ''}" if keys else "")
        )
    if task == _TASK_ANNOT:
        ref = [k for k, _ in detected["annotation"]["reference"]]
        pred = [k for k, _ in detected["annotation"]["predicted"]]
        batch = [k for k, _ in detected["batch"]]
        return (
            f"{len(ref)} ref / {len(pred)} predicted / {len(batch)} batch "
            f"column(s) in adata.obs"
        )
    if task == _TASK_DIMRED:
        embs = [k for k, _ in detected["clustering"]["embeddings"]]
        return (
            f"{len(embs)} embedding(s) in adata.obsm"
            + (f": {embs[:3]}{'...' if len(embs) > 3 else ''}" if embs else "")
        )
    return ""
