"""I/O helpers for the scFM QC layer.

This module bridges the gap between CheckAtlas' on-disk per-task
artefacts (the wide MultiQC-style TSVs written by
``checkatlas.atlas.create_metric_*``) and the long-format metric
table the scfm diagnostic engine expects.

    The convention in ``checkatlas_files/<task>/<atlas>.tsv``:

    * **cluster** — wide format, one row per ``(embedding, label_key)``,
      columns ``Clust_Sample, obs, <metric>, <metric>_running_time, ...``.
    * **annotation** — wide format, one row per
      ``(embedding, reference, predicted)``, columns
      ``Annot_Sample, Reference, obs, <metric>, <metric>_running_time, ...``.
    * **batch_correction** — wide format, one row per
      ``(embedding, batch_key)``, columns
      ``Batch_Sample, Embedding, Batch Key, <metric>,
      <metric>_running_time, ...``. The ``Metric Name`` in the wide
      header is the long-format metric name (e.g. ``iLISI``, ``cLISI``,
      ``kbet``, ``pcr``, ``graph_connectivity``).
    * **dimred** — wide format, one row per ``obsm`` key, columns
      ``Dimred_Sample, obsm, <metric>, <metric>_running_time, ...``.

    The function :func:`wide_to_long` converts any of these into the
    long format used by :func:`checkatlas.scfm.diagnostics.diagnose`::

        [Atlas Name, Task, Metric Name, Embedding,
         Reference/Input 1, Prediction/Input 2, Value, Time (s)]

NaN / empty cells are skipped — only computed metric values are
forwarded to the diagnostic engine. This keeps partial runs from
producing false-positive verdicts.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import pandas as pd

from ..utils import folders

logger = logging.getLogger("checkatlas")


def _tsv_path(path: str, task: str, atlas_name: str) -> str:
    """Return the canonical TSV path for one task under
    ``checkatlas_files/<task>/<atlas>.tsv``.
    """
    folder = folders.get_folder(path, task)
    return os.path.join(folder, f"{atlas_name}.tsv")


def load_per_task_tsvs(
    path: str,
    atlas_name: str,
    *,
    verbose: bool = True,
) -> dict[str, Optional[pd.DataFrame]]:
    """Load the three per-task TSVs from disk, returning ``None`` for
    any that are missing.

    Parameters
    ----------
    path : str
        The checkatlas working directory (the directory that contains
        ``checkatlas_files/``).
    atlas_name : str
    verbose : bool, default True
        If False, the per-TSV ``INFO`` "loaded N rows from ..." line
        is suppressed. Used by the orchestrator to avoid double-logging
        (it loads once internally and the caller loads again).

    Returns
    -------
    dict
        ``{"cluster": df|None, "annotation": df|None,
        "batch_correction": df|None, "dimred": df|None}``.
    """
    out: dict[str, Optional[pd.DataFrame]] = {}
    for task in ("cluster", "annotation", "batch_correction", "dimred"):
        p = _tsv_path(path, task, atlas_name)
        if not os.path.isfile(p):
            logger.debug("scfm: per-task TSV missing: %s", p)
            out[task] = None
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception as exc:
            logger.warning("scfm: failed to read %s: %s", p, exc)
            out[task] = None
            continue
        if df.empty:
            out[task] = None
            continue
        out[task] = df
        if verbose:
            logger.info(
                "scfm: loaded %d rows from %s",
                len(df),
                p,
            )
    return out


def _strip_atlas_prefix(sample: str, atlas_name: str) -> str:
    """Strip ``<atlas>_`` from a sample label, falling back to the
    full string if the prefix is not present.
    """
    if not isinstance(sample, str):
        return ""
    prefix = f"{atlas_name}_"
    if sample.startswith(prefix):
        return sample[len(prefix):]
    return sample


def _strip_label_suffix(rest: str, label: str) -> Optional[str]:
    """Strip a trailing ``_<label>`` suffix from ``rest``. Returns
    ``None`` if the suffix is not present (malformed row).
    """
    if not isinstance(label, str) or not label:
        return rest
    suffix = f"_{label}"
    if rest.endswith(suffix):
        return rest[: -len(suffix)]
    return None


def _emit_long_rows(
    atlas_name: str,
    task: str,
    metric_to_embedding: dict[str, str],
    embedding: str,
    reference: str,
    prediction: str,
    df: pd.DataFrame,
) -> list[dict]:
    """Emit one long-format row per non-NaN metric value in ``df``.

    ``df`` is a single-row dataframe (one wide-format row) holding
    alternating ``<metric>`` and ``<metric>_running_time`` columns.
    """
    rows: list[dict] = []
    for col in df.columns:
        if not col.endswith("_running_time"):
            continue
        metric = col[: -len("_running_time")]
        if metric not in df.columns:
            continue
        try:
            value = df[metric].iloc[0]
        except Exception:
            continue
        if pd.isna(value) or value == "":
            continue
        try:
            time_v = df[col].iloc[0]
        except Exception:
            time_v = 0.0
        try:
            time_val = float(time_v) if not pd.isna(time_v) else 0.0
        except (TypeError, ValueError):
            time_val = 0.0
        # Allow caller to override the embedding per metric
        emb = metric_to_embedding.get(metric, embedding)
        rows.append(
            {
                "Atlas Name": atlas_name,
                "Task": task,
                "Metric Name": metric,
                "Embedding": emb,
                "Reference/Input 1": reference,
                "Prediction/Input 2": prediction,
                "Value": float(value),
                "Time (s)": time_val,
            }
        )
    return rows


def wide_to_long(
    df: pd.DataFrame,
    kind: str,
    atlas_name: str,
) -> pd.DataFrame:
    """Unpivot a wide-format per-task TSV into the long format the
    scfm diagnostic engine consumes.

    Parameters
    ----------
    df : pd.DataFrame
        The wide-format TSV, as written by ``create_metric_*``.
    kind : str
        One of ``"cluster"``, ``"annotation"``,
        ``"batch_correction"``, ``"dimred"``.
    atlas_name : str

    Returns
    -------
    pd.DataFrame
        Long-format table with columns
        ``[Atlas Name, Task, Metric Name, Embedding,
        Reference/Input 1, Prediction/Input 2, Value, Time (s)]``.
    """
    if df is None or df.empty:
        return pd.DataFrame(
            columns=[
                "Atlas Name",
                "Task",
                "Metric Name",
                "Embedding",
                "Reference/Input 1",
                "Prediction/Input 2",
                "Value",
                "Time (s)",
            ]
        )

    columns = [
        "Atlas Name",
        "Task",
        "Metric Name",
        "Embedding",
        "Reference/Input 1",
        "Prediction/Input 2",
        "Value",
        "Time (s)",
    ]

    if kind == "annotation":
        if "Reference" not in df.columns or "obs" not in df.columns:
            logger.warning("scfm: annotation TSV missing Reference/obs columns")
            return pd.DataFrame(columns=columns)
        rows: list[dict] = []
        for _, r in df.iterrows():
            reference = str(r.get("Reference", "") or "")
            prediction = str(r.get("obs", "") or "")
            # ``Reference`` is the embedding key for the metric
            # (e.g. ``X_pca``, ``X_scvi``). Some TSV rows store the
            # actual reference column name here (e.g. ``cell_type``)
            # because the annotator ran the metric against a
            # reference column rather than an embedding. The
            # diagnostic engine looks up rows by the ``Embedding``
            # column, so we forward ``Reference`` as-is. Rows whose
            # ``Reference`` is a column name rather than an
            # ``obsm`` key are simply not matched by
            # ``diagnostics.diagnose`` and have no effect.
            embedding = reference
            rows.extend(
                _emit_long_rows(
                    atlas_name,
                    "annot",
                    {},
                    embedding,
                    reference,
                    prediction,
                    pd.DataFrame([r]).reset_index(drop=True),
                )
            )
        return pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)

    if kind == "cluster":
        if "obs" not in df.columns:
            logger.warning("scfm: cluster TSV missing obs column")
            return pd.DataFrame(columns=columns)
        rows = []
        for _, r in df.iterrows():
            label = str(r.get("obs", "") or "")
            sample = str(r.get("Clust_Sample", "") or "")
            rest = _strip_atlas_prefix(sample, atlas_name)
            embedding = _strip_label_suffix(rest, label) or rest
            rows.extend(
                _emit_long_rows(
                    atlas_name,
                    "cluster",
                    {},
                    embedding,
                    "",
                    label,
                    pd.DataFrame([r]).reset_index(drop=True),
                )
            )
        return pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)

    if kind == "batch_correction":
        if "Embedding" not in df.columns or "Batch Key" not in df.columns:
            logger.warning(
                "scfm: batch_correction TSV missing Embedding / Batch Key columns"
            )
            return pd.DataFrame(columns=columns)
        rows = []
        for _, r in df.iterrows():
            embedding = str(r.get("Embedding", "") or "")
            batch_key = str(r.get("Batch Key", "") or "")
            rows.extend(
                _emit_long_rows(
                    atlas_name,
                    "batch_correction",
                    {},
                    embedding,
                    embedding,
                    batch_key,
                    pd.DataFrame([r]).reset_index(drop=True),
                )
            )
        return (
            pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)
        )

    if kind == "dimred":
        if "obsm" not in df.columns:
            logger.warning("scfm: dimred TSV missing obsm column")
            return pd.DataFrame(columns=columns)
        rows = []
        for _, r in df.iterrows():
            emb = str(r.get("obsm", "") or "")
            rows.extend(
                _emit_long_rows(
                    atlas_name,
                    "dimred",
                    {},
                    emb,
                    "X",
                    "",
                    pd.DataFrame([r]).reset_index(drop=True),
                )
            )
        return pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)

    raise ValueError(f"scfm: unknown kind {kind!r}")
