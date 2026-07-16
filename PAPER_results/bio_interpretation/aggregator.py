"""
Data aggregation for CheckAtlas biological interpretation.

Reads all TSV files and log files, merges into unified DataFrames,
applies min-max scaling following the scIB protocol, and computes
aggregated per-atlas, per-embedding scores.
"""

import os
import re
from pathlib import Path
import pandas as pd
import numpy as np


# ── Paths ────────────────────────────────────────────────────────────────
RESULT_DIR = Path("/data/analysis/data_ganguly/checkatlas_files")
LOG_DIR = Path("/home/ganguly/checkatlas/PAPER_results/atlas_logs_8_core")
OUT_DIR = Path("/home/ganguly/checkatlas/PAPER_results/bio_interpretation")

CLUSTER_DIR = RESULT_DIR / "cluster"
ANNOT_DIR = RESULT_DIR / "annotation"
BATCH_DIR = RESULT_DIR / "batch_correction"
DIMRED_DIR = RESULT_DIR / "dimred"

# ── Atlas metadata (from log analysis) ───────────────────────────────────
ATLAS_META = {
    "bone_marrow": {"cells": 27112, "genes": 60606, "total_time_s": 1239, "has_cluster": True},
    "TS_neural": {"cells": 2685, "genes": 60606, "total_time_s": 67, "has_cluster": True},
    "liver": {"cells": 22214, "genes": 60606, "total_time_s": 879, "has_cluster": True},
    "blood": {"cells": 85233, "genes": 60606, "total_time_s": 5031, "has_cluster": True},
    "blood_pbmc": {"cells": 43512, "genes": 21649, "total_time_s": 734, "has_cluster": False},
    "lung": {"cells": 347970, "genes": 29824, "total_time_s": 1273, "has_cluster": False},
    "TS_skin": {"cells": 17786, "genes": 60606, "total_time_s": 439, "has_cluster": True},
    "TS_immune": {"cells": 592317, "genes": 60606, "total_time_s": 5259, "has_cluster": True},
    "lung_2M": {"cells": 2282447, "genes": 55329, "total_time_s": 7219, "has_cluster": False},
}

ATLAS_ORDER = ["TS_neural", "TS_skin", "bone_marrow", "liver", "blood_pbmc", "blood", "lung", "TS_immune", "lung_2M"]

# ── Metric definitions ───────────────────────────────────────────────────
CLUSTER_METRICS = ["silhouette", "calinski_harabasz", "davies_bouldin", "dbcvi", "kolmogorov_smirnov"]
CLUSTER_METRICS_DIRECTION = {"silhouette": 1, "calinski_harabasz": 1, "davies_bouldin": -1, "dbcvi": 1, "kolmogorov_smirnov": 1}

BATCH_METRICS = ["iLISI", "cLISI", "kbet", "pcr", "graph_connectivity"]
BATCH_METRICS_DIRECTION = {"iLISI": 1, "cLISI": -1, "kbet": -1, "pcr": -1, "graph_connectivity": 1}

ANNOT_METRICS = ["adj_rand_index", "adj_mutual_info", "normalized_mutual_info", "fowlkes_mallows", "vmeasure", "isolated_f1_score"]
ANNOT_METRICS_DIRECTION = {"adj_rand_index": 1, "adj_mutual_info": 1, "normalized_mutual_info": 1, "fowlkes_mallows": 1, "vmeasure": 1, "isolated_f1_score": 1}

DIMRED_METRICS = ["kruskal_stress", "spearman_rho", "dCor", "coknn", "entourage", "lcmc", "avg_jaccard_dis", "ged", "den_pre"]
DIMRED_METRICS_DIRECTION = {"kruskal_stress": -1, "spearman_rho": 1, "dCor": 1, "coknn": 1, "entourage": 1, "lcmc": 1, "avg_jaccard_dis": -1, "ged": -1, "den_pre": 1}


def _extract_embedding(sample_str):
    """Extract embedding name from Clust_Sample/Annot_Sample/etc."""
    parts = sample_str.rsplit("_", 1)
    atlas_embed = parts[0]
    atlas_parts = atlas_embed.split("_", 1)
    if len(atlas_parts) > 1:
        return atlas_parts[1]
    return atlas_embed


def load_cluster_metrics() -> pd.DataFrame:
    """Load all cluster TSVs into long-format DataFrame."""
    rows = []
    for tsv in CLUSTER_DIR.glob("*.tsv"):
        atlas = tsv.stem
        if atlas == "lung":
            continue  # lung has no cluster metrics
        df = pd.read_csv(tsv, sep="\t")
        for _, row_df in df.iterrows():
            sample = row_df["Clust_Sample"]
            embedding = _extract_embedding(sample)
            label_key = row_df["obs"]
            for metric in CLUSTER_METRICS:
                val = row_df.get(metric, np.nan)
                time_val = row_df.get(f"{metric}_running_time", np.nan)
                if pd.notna(val):
                    rows.append({
                        "atlas": atlas, "task": "cluster", "metric": metric,
                        "embedding": embedding, "label_key": label_key,
                        "value": val, "time_s": time_val,
                    })
    return pd.DataFrame(rows)


def load_annot_metrics() -> pd.DataFrame:
    """Load all annotation TSVs into long-format."""
    rows = []
    for tsv in ANNOT_DIR.glob("*.tsv"):
        atlas = tsv.stem
        df = pd.read_csv(tsv, sep="\t")
        for _, row_df in df.iterrows():
            sample = row_df.get("Annot_Sample", "")
            ref = row_df.get("Reference", "")
            pred = row_df.get("obs", "")

            for metric in ANNOT_METRICS:
                time_col = f"{metric}_running_time"
                val = row_df.get(metric, np.nan)
                if pd.notna(val) and time_col in df.columns:
                    time_val = row_df.get(time_col, np.nan)
                    if ref and not ref.startswith("X_"):
                        rows.append({
                            "atlas": atlas, "task": "annotation", "metric": metric,
                            "ref": ref, "pred": pred,
                            "value": val, "time_s": time_val,
                        })
    return pd.DataFrame(rows)


def load_batch_metrics() -> pd.DataFrame:
    """Load all batch correction TSVs into long-format."""
    rows = []
    for tsv in BATCH_DIR.glob("*.tsv"):
        atlas = tsv.stem
        df = pd.read_csv(tsv, sep="\t")
        for _, row_df in df.iterrows():
            embedding = str(row_df.get("Embedding", "")).replace(f"{atlas}_", "")
            batch_key = str(row_df.get("Batch Key", ""))
            for metric in BATCH_METRICS:
                val = row_df.get(metric, np.nan)
                if pd.notna(val):
                    rows.append({
                        "atlas": atlas, "task": "batch_correction", "metric": metric,
                        "embedding": embedding, "batch_key": batch_key,
                        "value": val,
                    })
    return pd.DataFrame(rows)


def load_dimred_metrics() -> pd.DataFrame:
    """Load all dimred TSVs into long-format."""
    rows = []
    for tsv in DIMRED_DIR.glob("*.tsv"):
        atlas = tsv.stem
        if atlas == "pancreas_scib":
            continue  # exclude scIB test atlas
        df = pd.read_csv(tsv, sep="\t")
        for _, row_df in df.iterrows():
            sample = row_df["Dimred_Sample"]
            embedding = _extract_embedding(sample)
            for metric in DIMRED_METRICS:
                val = row_df.get(metric, np.nan)
                if pd.notna(val):
                    rows.append({
                        "atlas": atlas, "task": "dimred", "metric": metric,
                        "embedding": embedding,
                        "value": val,
                    })
    return pd.DataFrame(rows)


def _is_bio_batch_key(key: str) -> bool:
    """Heuristic: does this batch key represent a genuine technical batch effect?"""
    key_lower = key.lower()
    biological = {"cell_type", "free_annotation", "scvi_leiden", "scanvi_label"}
    for b in biological:
        if b in key_lower:
            return False
    return True


def minmax_scale_per_metric(df: pd.DataFrame, metrics: list, direction: dict) -> pd.DataFrame:
    """Apply min-max scaling per metric across all atlases, flipping direction so higher=better."""
    df = df.copy()
    for metric in metrics:
        mask = df["metric"] == metric
        vals = df.loc[mask, "value"].copy()
        if len(vals.dropna()) < 2:
            df.loc[mask, "scaled"] = 0.5
            continue
        vmin, vmax = vals.min(), vals.max()
        if vmax == vmin:
            df.loc[mask, "scaled"] = 0.5
        else:
            scaled = (vals - vmin) / (vmax - vmin)
            if direction.get(metric, 1) == -1:
                scaled = 1.0 - scaled
            df.loc[mask, "scaled"] = scaled
    return df


def aggregate_all() -> dict:
    """Load, merge, scale all metrics. Return dict of DataFrames."""
    cluster = load_cluster_metrics()
    annot = load_annot_metrics()
    batch = load_batch_metrics()
    dimred = load_dimred_metrics()

    # Fix embedding names: strip atlas prefix
    for df_part in [cluster, batch, dimred]:
        df_part["embedding_clean"] = df_part["embedding"].apply(
            lambda x: re.sub(r'^(blood_pbmc_|lung_2M_|lung_|TS_immune_|TS_skin_|TS_neural_|bone_marrow_|liver_|blood_)', '', str(x))
        )

    df_all = pd.concat([cluster, annot, batch, dimred], ignore_index=True)

    # Scale per task
    cluster = minmax_scale_per_metric(cluster, CLUSTER_METRICS, CLUSTER_METRICS_DIRECTION)
    annot = minmax_scale_per_metric(annot, ANNOT_METRICS, ANNOT_METRICS_DIRECTION)
    batch = minmax_scale_per_metric(batch, BATCH_METRICS, BATCH_METRICS_DIRECTION)
    dimred = minmax_scale_per_metric(dimred, DIMRED_METRICS, DIMRED_METRICS_DIRECTION)

    # Filter batch to technical batch keys only
    batch_tech = batch[batch["batch_key"].apply(_is_bio_batch_key)].copy()
    batch_bio = batch[~batch["batch_key"].apply(_is_bio_batch_key)].copy()

    # Per-embedding aggregate scores (for batch-vs-bio)
    batch_embed_scores = batch_tech.groupby(["atlas", "embedding", "metric"]).agg(
        value_mean=("value", "mean"), scaled_mean=("scaled", "mean")
    ).reset_index()

    # cLISI per bio-key
    if "cLISI" in batch_bio["metric"].values:
        clisi_bio = batch_bio[batch_bio["metric"] == "cLISI"].groupby(
            ["atlas", "embedding"]
        ).agg(cLISI_mean=("value", "mean")).reset_index()
    else:
        clisi_bio = pd.DataFrame(columns=["atlas", "embedding", "cLISI_mean"])

    # iLISI per tech key
    ilisi_tech = batch_tech[batch_tech["metric"] == "iLISI"].groupby(
        ["atlas", "embedding"]
    ).agg(iLISI_mean=("value", "mean")).reset_index()

    # Aggregate per-task scores
    def _task_score(df: pd.DataFrame, task_name: str) -> pd.DataFrame:
        return df.groupby("atlas").agg(
            **{f"{task_name}_score": ("scaled", "mean")}
        ).reset_index()

    cluster_score = _task_score(cluster, "cluster")
    annot_score = _task_score(annot, "annotation")
    batch_score = _task_score(batch_tech, "batch")
    dimred_score = _task_score(dimred, "dimred")

    # Per-atlas summary
    atlas_summary = pd.DataFrame({"atlas": ATLAS_ORDER})
    atlas_summary["cells"] = atlas_summary["atlas"].map(lambda a: ATLAS_META.get(a, {}).get("cells", 0))
    atlas_summary["genes"] = atlas_summary["atlas"].map(lambda a: ATLAS_META.get(a, {}).get("genes", 0))
    atlas_summary["total_time_s"] = atlas_summary["atlas"].map(lambda a: ATLAS_META.get(a, {}).get("total_time_s", 0))
    atlas_summary["has_cluster"] = atlas_summary["atlas"].map(lambda a: ATLAS_META.get(a, {}).get("has_cluster", False))

    for score_df in [cluster_score, annot_score, batch_score, dimred_score]:
        atlas_summary = atlas_summary.merge(score_df, on="atlas", how="left")

    # Overall score: 0.25 per task for atlases missing cluster/annot
    for col in ["cluster_score", "annotation_score", "batch_score", "dimred_score"]:
        atlas_summary[col] = atlas_summary[col].fillna(0)

    atlas_summary["overall_score"] = 0.25 * (
        atlas_summary["cluster_score"] +
        atlas_summary["annotation_score"] +
        atlas_summary["batch_score"] +
        atlas_summary["dimred_score"]
    )

    return {
        "cluster": cluster,
        "annot": annot,
        "batch": batch,
        "batch_tech": batch_tech,
        "batch_bio": batch_bio,
        "dimred": dimred,
        "df_all": df_all,
        "atlas_summary": atlas_summary,
        "batch_embed_scores": batch_embed_scores,
        "clisi_bio": clisi_bio,
        "ilisi_tech": ilisi_tech,
        "atlas_meta": ATLAS_META,
    }


if __name__ == "__main__":
    data = aggregate_all()
    for k, v in data.items():
        if isinstance(v, pd.DataFrame):
            print(f"{k}: {v.shape}")
    print("\nAtlas Summary:")
    print(data["atlas_summary"].to_string())
