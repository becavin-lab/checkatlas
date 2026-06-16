"""Rare-type recall metric (Problem 5).

Reproduces the Atti & Subramaniam (2025) and DenAdel et al. (2025)
findings that scFM tokenization / over-parameterisation leads to
*overprediction of common cell types* and *underprediction of rare
types* (the scGPT CD14 overprediction finding).

The metric computes per-class F1 on rare-vs-common cells, and
returns the (rare_F1, common_F1, gap) triple plus a long-format row
suitable for the existing CheckAtlas pipeline.
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import f1_score

logger = logging.getLogger("checkatlas")


def _safe_str(arr) -> np.ndarray:
    try:
        return np.asarray(arr).astype(str)
    except Exception:
        return np.asarray(arr, dtype=object).astype(str)


def run(
    adata: AnnData,
    embeddings: Sequence[str],
    *,
    pred_label: str,
    ref_label: str,
    rare_quantile: float = 0.10,
) -> pd.DataFrame:
    """Per-class F1 on rare and common cell types.

    Parameters
    ----------
    adata : AnnData
    embeddings : sequence of str
        Embedding keys. (Not used directly for the F1 itself — the
        score is label-based — but kept for symmetry with the rest of
        the metric pipeline so it can be pivoted with embedding on
        the row axis in MultiQC.)
    pred_label, ref_label : str
        Predicted and reference cell-type columns.
    rare_quantile : float
        Cells whose class frequency is in the bottom ``rare_quantile``
        are "rare"; the rest are "common". Default 0.10.

    Returns
    -------
    pd.DataFrame
        Long format with columns
        ``[Embedding, Metric, Subset, Value]`` where
        ``Metric in {"rare_f1", "common_f1", "rare_common_gap"}``.
    """
    if pred_label not in adata.obs.columns or ref_label not in adata.obs.columns:
        return pd.DataFrame(
            columns=["Embedding", "Metric", "Subset", "Value"]
        )

    y_true = _safe_str(adata.obs[ref_label].values)
    y_pred = _safe_str(adata.obs[pred_label].values)

    classes, counts = np.unique(y_true, return_counts=True)
    if classes.size < 2:
        return pd.DataFrame(
            columns=["Embedding", "Metric", "Subset", "Value"]
        )

    # Rare = the bottom ``rare_quantile`` *fraction* of classes when
    # sorted by frequency (ascending). If ``rare_quantile=0.34`` and
    # there are 3 classes, the rarest 1 class is "rare" and the rest
    # are "common". This matches the Souza & Mehta / Atti 2025
    # definition of "rare" in the scFM literature.
    order = np.argsort(counts)
    n_rare = max(1, int(np.floor(len(classes) * rare_quantile)))
    rare_classes = classes[order[:n_rare]]
    common_classes = classes[order[n_rare:]]

    rare_mask = np.isin(y_true, rare_classes)
    common_mask = np.isin(y_true, common_classes)

    rows: list[dict] = []
    for emb in embeddings:
        if emb not in adata.obsm:
            continue
        try:
            rare_f1 = float(
                f1_score(
                    y_true[rare_mask],
                    y_pred[rare_mask],
                    average="weighted",
                    zero_division=0,
                )
            )
        except Exception:
            rare_f1 = float("nan")
        try:
            common_f1 = float(
                f1_score(
                    y_true[common_mask],
                    y_pred[common_mask],
                    average="weighted",
                    zero_division=0,
                )
            )
        except Exception:
            common_f1 = float("nan")
        gap = (
            float(common_f1 - rare_f1)
            if (not np.isnan(rare_f1) and not np.isnan(common_f1))
            else float("nan")
        )
        for stat, val in (
            ("rare_f1", rare_f1),
            ("common_f1", common_f1),
            ("rare_common_gap", gap),
        ):
            rows.append(
                {
                    "Embedding": emb,
                    "Metric": stat,
                    "Subset": "all",
                    "Value": val,
                }
            )
    return pd.DataFrame(rows)
