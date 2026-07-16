"""
Table generation for CheckAtlas biological interpretation.
Produces LaTeX tables for the manuscript.
"""

import pandas as pd
import numpy as np
from pathlib import Path

OUT_DIR = Path("/home/ganguly/checkatlas/PAPER_results/bio_interpretation")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _fmt3(x):
    """Format to 3 decimal places."""
    if pd.isna(x):
        return "---"
    return f"{x:.3f}"


def _fmt_cells(n):
    """Human-readable cell count."""
    if n >= 1_000_000:
        return f"{n/1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n/1_000:.0f}k"
    return str(n)


def _fmt_s(seconds):
    """Human-readable runtime."""
    if seconds < 60:
        return f"{int(seconds)}s"
    if seconds < 3600:
        return f"{int(seconds // 60)}m {int(seconds % 60)}s"
    return f"{int(seconds // 3600)}h {int((seconds % 3600) // 60)}m"


def table1_atlas_summary(data: dict) -> str:
    """Atlas-level summary table (Table 1)."""
    df = data["atlas_summary"].copy()
    atlas_meta = data["atlas_meta"]

    rows = []
    for _, row in df.iterrows():
        atlas = row["atlas"]
        meta = atlas_meta.get(atlas, {})
        rows.append({
            "Atlas": atlas.replace("_", "\\_"),
            "Cells": _fmt_cells(meta.get("cells", 0)),
            "Genes": f"{meta.get('genes', 0):,}",
            "Cluster": _fmt3(row["cluster_score"]),
            "Annotation": _fmt3(row["annotation_score"]),
            "Batch": _fmt3(row["batch_score"]),
            "Dimred": _fmt3(row["dimred_score"]),
            "Overall": _fmt3(row["overall_score"]),
            "Runtime": _fmt_s(int(row["total_time_s"])),
        })

    result_df = pd.DataFrame(rows)

    latex = r"""
\begin{table}[htbp]
\centering
\caption{\textbf{Atlas-level summary of CheckAtlas benchmark results.}
Scores are min-max scaled per metric within each task (0=worst, 1=best).
Overall score is the unweighted mean across the four tasks.
Atlas names are abbreviated; \textit{TS\_} prefix indicates Tabula Sapiens atlas.
\emph{lung\_2M} (2.3 million cells) is the largest atlas with partial metrics.}
\label{tab:atlas_summary}
\small
\begin{tabular}{l r r c c c c c c}
\toprule
\textbf{Atlas} & \textbf{Cells} & \textbf{Genes} & \textbf{Cluster} & \textbf{Annotation} & \textbf{Batch} & \textbf{Dimred} & \textbf{Overall} & \textbf{Runtime} \\
\midrule
"""

    for _, row in result_df.iterrows():
        latex += f"{row['Atlas']} & {row['Cells']} & {row['Genes']} & {row['Cluster']} & {row['Annotation']} & {row['Batch']} & {row['Dimred']} & {row['Overall']} & {row['Runtime']} \\\\\n"

    latex += r"""\bottomrule
\end{tabular}
\end{table}
"""
    return latex


def table2_batch_vs_bio(data: dict) -> str:
    """Per-embedding batch-vs-biology table (Table 2)."""
    clisi = data["clisi_bio"].copy()
    ilisi = data["ilisi_tech"].copy()
    df = data["batch_tech"].copy()

    if clisi.empty or ilisi.empty:
        return "No batch-vs-biology data available."

    merged = ilisi.merge(clisi, on=["atlas", "embedding"], how="inner")

    # Add mean kbet per embedding
    kbet_mean = df[df["metric"] == "kbet"].groupby(["atlas", "embedding"])["value"].mean().reset_index()
    kbet_mean.columns = ["atlas", "embedding", "kbet_mean"]

    # Add mean PCR per embedding
    pcr_mean = df[df["metric"] == "pcr"].groupby(["atlas", "embedding"])["value"].mean().reset_index()
    pcr_mean.columns = ["atlas", "embedding", "pcr_mean"]

    merged = merged.merge(kbet_mean, on=["atlas", "embedding"], how="left")
    merged = merged.merge(pcr_mean, on=["atlas", "embedding"], how="left")

    # Shorten embedding names
    merged["embedding_short"] = merged["embedding"].apply(
        lambda x: x[:30] + "..." if len(str(x)) > 30 else str(x)
    )

    rows = []
    for _, row in merged.iterrows():
        batch_score = (row["iLISI_mean"] + (1 - row["kbet_mean"]) + (1 - row["pcr_mean"])) / 3
        bio_score = 1 - row["cLISI_mean"]
        overall = 0.5 * batch_score + 0.5 * bio_score
        rows.append({
            "Atlas": row["atlas"].replace("_", "\\_"),
            "Embedding": row["embedding_short"],
            "iLISI": _fmt3(row["iLISI_mean"]),
            "cLISI": _fmt3(row["cLISI_mean"]),
            "kBET": _fmt3(row["kbet_mean"]),
            "PCR": _fmt3(row["pcr_mean"]),
            "Batch": _fmt3(batch_score),
            "Bio": _fmt3(bio_score),
            "Overall": _fmt3(overall),
        })

    result_df = pd.DataFrame(rows)

    latex = r"""
\begin{table}[htbp]
\centering
\caption{\textbf{Batch-vs-biology tradeoff across atlas embeddings.}
iLISI measures batch mixing (higher=better), cLISI measures cell-type purity (lower=better).
kBET and PCR measure specific aspects of batch effect removal.
Batch score = mean(iLISI, 1-kBET, 1-PCR). Bio score = 1-cLISI. Overall = mean(Batch, Bio).}
\label{tab:batch_vs_bio}
\small
\begin{tabular}{l l c c c c c c c}
\toprule
\textbf{Atlas} & \textbf{Embedding} & \textbf{iLISI$\uparrow$} & \textbf{cLISI$\downarrow$} & \textbf{kBET$\downarrow$} & \textbf{PCR$\downarrow$} & \textbf{Batch} & \textbf{Bio} & \textbf{Overall} \\
\midrule
"""

    for _, row in result_df.iterrows():
        latex += f"{row['Atlas']} & {row['Embedding']} & {row['iLISI']} & {row['cLISI']} & {row['kBET']} & {row['PCR']} & {row['Batch']} & {row['Bio']} & {row['Overall']} \\\\\n"

    latex += r"""\bottomrule
\end{tabular}
\end{table}
"""
    return latex


def table3_dimred_summary(data: dict) -> str:
    """Dimred per-atlas summary (Table 3)."""
    df = data["dimred"].copy()

    metrics_global = ["spearman_rho", "dCor"]
    metrics_local = ["coknn", "entourage", "lcmc", "den_pre", "kruskal_stress", "avg_jaccard_dis", "ged"]

    rows = []
    for atlas in sorted(df["atlas"].unique()):
        atlas_data = df[df["atlas"] == atlas]
        global_vals = atlas_data[atlas_data["metric"].isin(metrics_global)]
        local_vals = atlas_data[atlas_data["metric"].isin(metrics_local)]

        global_score = global_vals["scaled"].mean() if "scaled" in global_vals.columns and len(global_vals) > 0 else np.nan
        local_score = local_vals["scaled"].mean() if "scaled" in local_vals.columns and len(local_vals) > 0 else np.nan

        n_emb = atlas_data["embedding"].nunique()

        rows.append({
            "Atlas": atlas.replace("_", "\\_"),
            "Embeddings": n_emb,
            "Global": _fmt3(global_score),
            "Local": _fmt3(local_score),
            "Spearman": _fmt3(atlas_data[atlas_data["metric"] == "spearman_rho"]["value"].mean()),
            "dCor": _fmt3(atlas_data[atlas_data["metric"] == "dCor"]["value"].mean()),
            "CoKNN": _fmt3(atlas_data[atlas_data["metric"] == "coknn"]["value"].mean()),
            "LCMC": _fmt3(atlas_data[atlas_data["metric"] == "lcmc"]["value"].mean()),
        })

    result_df = pd.DataFrame(rows)

    latex = r"""
\begin{table}[htbp]
\centering
\caption{\textbf{Dimensionality reduction quality per atlas.}
Global = mean(scaled spearman\_rho, dCor). Local = mean(scaled coKNN, entourage, lcmc, den\_pre, kruskal\_stress, avg\_jaccard\_dis, ged).
Higher is always better after direction-aware min-max scaling.}
\label{tab:dimred_summary}
\small
\begin{tabular}{l c c c c c c c}
\toprule
\textbf{Atlas} & \textbf{Emb.} & \textbf{Global} & \textbf{Local} & \textbf{Spearman$\rho$} & \textbf{dCor} & \textbf{CoKNN} & \textbf{LCMC} \\
\midrule
"""

    for _, row in result_df.iterrows():
        latex += f"{row['Atlas']} & {row['Embeddings']} & {row['Global']} & {row['Local']} & {row['Spearman']} & {row['dCor']} & {row['CoKNN']} & {row['LCMC']} \\\\\n"

    latex += r"""\bottomrule
\end{tabular}
\end{table}
"""
    return latex


def generate_all_tables(data: dict):
    """Generate all tables."""
    with open(OUT_DIR / "Table1_atlas_summary.tex", "w") as f:
        f.write(table1_atlas_summary(data))
    print("Table 1 written to Table1_atlas_summary.tex")

    with open(OUT_DIR / "Table2_batch_vs_bio.tex", "w") as f:
        f.write(table2_batch_vs_bio(data))
    print("Table 2 written to Table2_batch_vs_bio.tex")

    with open(OUT_DIR / "Table3_dimred_summary.tex", "w") as f:
        f.write(table3_dimred_summary(data))
    print("Table 3 written to Table3_dimred_summary.tex")

    # Also write CSV for easy inspection
    data["atlas_summary"].to_csv(OUT_DIR / "atlas_summary.csv", index=False)
    print("CSV: atlas_summary.csv")


if __name__ == "__main__":
    from aggregator import aggregate_all
    data = aggregate_all()
    generate_all_tables(data)
