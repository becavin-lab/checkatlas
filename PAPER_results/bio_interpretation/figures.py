"""
Publication-quality figure generation for CheckAtlas benchmarking.

Follows scIB figure conventions (Luecken et al., Nature Methods, 2022)
and Nature Methods formatting guidelines.

Colour scheme:
  - Batch removal: blues (#2166ac, #4393c3, #92c5de)
  - Bio-conservation: pinks/reds (#b2182b, #d6604d, #f4a582)
  - Dimred: greens (#1b7837, #5aae61, #a6dba0)
  - Cluster: purples (#762a83, #9970ab, #c2a5cf)
  - Overall: dark grey (#404040)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

OUT_DIR = Path("/home/ganguly/checkatlas/PAPER_results/bio_interpretation/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Nature Methods style ─────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "lines.linewidth": 0.8,
})

# ── Consistent colour palette ────────────────────────────────────────────
C_BATCH = "#2166ac"       # Batch removal - blue
C_BIO = "#b2182b"         # Bio-conservation - red
C_DIMRED = "#1b7837"      # Dimred - green
C_CLUSTER = "#762a83"     # Cluster - purple
C_OVERALL = "#404040"     # Overall - grey
C_SCALIBILITY = "#e08214" # Scalability - orange

COLORS_TASK = {
    "cluster": C_CLUSTER,
    "annotation": C_BIO,
    "batch_correction": C_BATCH,
    "dimred": C_DIMRED,
}

EMBEDDING_COLORS = {
    "X_pca": "#2166ac",
    "X_scvi": "#b2182b",
    "X_umap": "#1b7837",
    "X_tsne": "#762a83",
    "X_umap_scvi_full_donorassay": "#e08214",
    "X_umap_tissue_scvi_donorassay": "#5aae61",
    "X_tissue_uncorrected_umap": "#9970ab",
    "X_uncorrected_umap": "#d6604d",
    "X_compartment_uncorrected_umap": "#4393c3",
    "X_umap_compartment_scvi_donorassay": "#f4a582",
    "X_scanvi_emb": "#c51b7d",
}

ATLAS_LABELS = {
    "TS_neural": "Neural\n(2.7k)",
    "TS_skin": "Skin\n(18k)",
    "bone_marrow": "Bone Marrow\n(27k)",
    "liver": "Liver\n(22k)",
    "blood_pbmc": "PBMC\n(44k)",
    "blood": "Blood\n(85k)",
    "lung": "Lung\n(348k)",
    "TS_immune": "Immune\n(592k)",
    "lung_2M": "Lung 2M\n(2.3M)",
}


def _fmt_atlas(label: str) -> str:
    return ATLAS_LABELS.get(label, label)


def _save(fig, name: str, dpi: int = 300):
    """Save figure as PNG and PDF."""
    for fmt, ext in [("png", "png"), ("pdf", "pdf")]:
        path = OUT_DIR / f"{name}.{ext}"
        fig.savefig(path, dpi=dpi, format=ext, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Overview heatmap: all metrics per atlas (scIB Fig 2a style)
# ═══════════════════════════════════════════════════════════════════════════

def fig1_overview_heatmap(data: dict):
    """
    Heatmap showing scaled metric scores per atlas per task.
    Rows = Atlases, Columns grouped by Task with per-metric subcolumns.
    """
    atlas_summary = data["atlas_summary"].copy()
    df_cluster = data["cluster"]
    df_annot = data["annot"]
    df_batch = data["batch_tech"]
    df_dimred = data["dimred"]

    atlases = [a for a in ATLAS_LABELS if a in atlas_summary["atlas"].values]

    # Build matrix: rows=atlas, columns=metric
    all_metrics = []
    metric_tasks = {}
    for task, metrics, df in [
        ("Cluster", ["silhouette", "calinski_harabasz", "davies_bouldin", "dbcvi", "kolmogorov_smirnov"], df_cluster),
        ("Annotation", ["adj_rand_index", "adj_mutual_info", "normalized_mutual_info", "fowlkes_mallows", "vmeasure", "isolated_f1_score", "average_silhouette_width", "dunn_index"], df_annot),
        ("Batch", ["iLISI", "cLISI", "kbet", "pcr", "graph_connectivity"], df_batch),
        ("Dimred", ["spearman_rho", "dCor", "coknn", "entourage", "lcmc", "den_pre", "kruskal_stress", "avg_jaccard_dis", "ged"], df_dimred),
    ]:
        for m in metrics:
            metric_id = f"{task[:4]}_{m}"
            all_metrics.append(metric_id)
            metric_tasks[metric_id] = task

    n_atlases = len(atlases)
    n_metrics = len(all_metrics)

    fig, ax = plt.subplots(figsize=(16, 5))

    heatmap_data = np.full((n_atlases, n_metrics), np.nan)
    for i, atlas in enumerate(atlases):
        for j, mid in enumerate(all_metrics):
            task = metric_tasks[mid]
            metric = mid[5:]  # Remove prefix
            if task == "Cluster":
                subset = df_cluster[df_cluster["atlas"] == atlas]
                subset = subset[subset["metric"] == metric]
            elif task == "Annotation":
                subset = df_annot[df_annot["atlas"] == atlas]
                subset = subset[subset["metric"] == metric]
            elif task == "Batch":
                subset = df_batch[df_batch["atlas"] == atlas]
                subset = subset[subset["metric"] == metric]
            else:
                subset = df_dimred[df_dimred["atlas"] == atlas]
                subset = subset[subset["metric"] == metric]

            if len(subset) > 0 and "scaled" in subset.columns:
                heatmap_data[i, j] = subset["scaled"].mean()

    mask = np.isnan(heatmap_data)

    cmap = sns.diverging_palette(10, 250, s=70, l=50, center="light", as_cmap=True)

    im = ax.imshow(heatmap_data, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")

    # Overlay NaN with hatched grey
    for i in range(n_atlases):
        for j in range(n_metrics):
            if mask[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=True, facecolor="#e0e0e0",
                                           edgecolor="#cccccc", linewidth=0.3, hatch="///"))

    # Tick labels
    short_metrics = [m.replace("calinski_harabasz", "CH").replace("davies_bouldin", "DB")
                      .replace("kolmogorov_smirnov", "KS").replace("adj_rand_index", "ARI")
                      .replace("adj_mutual_info", "AMI").replace("normalized_mutual_info", "NMI")
                      .replace("fowlkes_mallows", "FMI").replace("isolated_f1_score", "IsolF1")
                      .replace("average_silhouette_width", "ASW").replace("dunn_index", "Dunn")
                      .replace("graph_connectivity", "Graph").replace("spearman_rho", "Spearman")
                      .replace("kruskal_stress", "Stress").replace("avg_jaccard_dis", "Jaccard")
                      .replace("entourage", "Entour").replace("den_pre", "DenPre")
                      .replace("batch_correction", "Batch").replace("_", "\n")
                      for m in all_metrics]

    ax.set_xticks(range(n_metrics))
    ax.set_xticklabels(short_metrics, fontsize=5.5, rotation=90, ha="center")
    ax.set_yticks(range(n_atlases))
    ax.set_yticklabels([_fmt_atlas(a) for a in atlases], fontsize=6.5)

    # Task grouping lines
    task_boundaries = []
    prev = None
    for j, mid in enumerate(all_metrics):
        t = metric_tasks[mid]
        if t != prev:
            task_boundaries.append(j)
        prev = t
    task_boundaries.append(n_metrics)

    task_names = []
    for j in range(len(task_boundaries) - 1):
        mid = (task_boundaries[j] + task_boundaries[j + 1]) / 2
        t = metric_tasks[all_metrics[task_boundaries[j]]]
        ax.text(mid, -1.5, t, fontsize=6.5, fontweight="bold", ha="center",
                va="bottom", color=COLORS_TASK.get(t.lower(), "#333"))

    for b in task_boundaries[1:-1]:
        ax.axvline(x=b - 0.5, color="white", linewidth=1.5)

    # Colour bar
    cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.02, shrink=0.75)
    cbar.set_label("Scaled score (0–1)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.5)

    ax.set_title("CheckAtlas benchmarking overview — scaled metric scores per atlas", fontsize=8, fontweight="bold", pad=10)
    ax.set_xlabel("Metrics (grouped by task)", fontsize=7)

    _save(fig, "Fig1_overview_heatmap")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — iLISI vs cLISI scatter: batch-vs-biology frontier (scIB Fig 3a)
# ═══════════════════════════════════════════════════════════════════════════

def fig2_ilisi_vs_clisi(data: dict):
    """
    Scatter plot: iLISI (batch mixing) vs 1-cLISI (cell-type purity) per atlas per embedding.
    Shows the batch-vs-biology Pareto frontier. Both axes: higher=better.
    """
    clisi = data["clisi_bio"].copy()
    ilisi = data["ilisi_tech"].copy()
    atlas_meta = data["atlas_meta"]

    merged = ilisi.merge(clisi, on=["atlas", "embedding"], how="inner")
    if merged.empty:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.text(0.5, 0.5, "No data available for iLISI/cLISI comparison", ha="center", va="center", transform=ax.transAxes)
        _save(fig, "Fig2_ilisi_clisi_scatter")
        return fig

    merged["cells"] = merged["atlas"].map(lambda a: atlas_meta.get(a, {}).get("cells", 0))
    # Convert cLISI to purity: 1 - cLISI (higher = better cell-type preservation)
    merged["purity"] = 1.0 - merged["cLISI_mean"]

    merged["embedding_short"] = merged["embedding"].apply(
        lambda x: x.replace("_scvi_full_donorassay", "_scvi_full")
        .replace("_tissue_scvi_donorassay", "_tissue_scvi")
        .replace("_compartment_scvi_donorassay", "_comp_scvi")
        .replace("_uncorrected_umap", "_unco").replace("tissue_uncorrected_umap", "X_tissue_unco")
        .replace("_umap_compartment_scvi_donorassay", "_comp_scvi")
        .replace("_umap_scvi_full_donorassay", "_scvi_full")
    )

    n_col = min(3, merged["atlas"].nunique())
    n_row = int(np.ceil(merged["atlas"].nunique() / n_col))

    fig, axes = plt.subplots(n_row, n_col, figsize=(3.3 * n_col, 2.8 * n_row), squeeze=False,
                            constrained_layout=True)
    axes = axes.flatten()

    atlases_sorted = sorted(merged["atlas"].unique(), key=lambda a: atlas_meta.get(a, {}).get("cells", 0))

    for idx, atlas in enumerate(atlases_sorted):
        ax = axes[idx]
        atlas_data = merged[merged["atlas"] == atlas]

        for _, row in atlas_data.iterrows():
            emb = row["embedding_short"]
            color = EMBEDDING_COLORS.get(row["embedding"], "#999999")
            size = max(30, min(140, np.sqrt(float(row["cells"])) / 6))
            ax.scatter(row["iLISI_mean"], row["purity"], c=color, s=size,
                      edgecolors="white", linewidth=0.4, zorder=3, alpha=0.9,
                      label=emb)

        # Ideal region = top-right (high iLISI, high purity)
        ax.fill_between([0.2, 1.0], 0.95, 1.0, alpha=0.08, color="#2166ac")

        # Set limits
        x_max = max(atlas_data["iLISI_mean"].max() * 1.15, 0.3)
        ax.set_xlim(-0.02, min(x_max, 0.85))
        ax.set_ylim(0.94, 1.002)

        cells = atlas_meta.get(atlas, {}).get("cells", 0)
        ax.set_title(f"{atlas} ({cells:,} cells)", fontsize=7, fontweight="bold")
        ax.set_xlabel("iLISI (batch mixing)", fontsize=6.5)
        ax.set_ylabel("1 − cLISI (cell-type purity)", fontsize=6.5)

        # Add diagonal reference
        ax.plot([0, 0.85], [0.95, 0.95], "k-", linewidth=0.3, alpha=0.15)

    for idx in range(len(atlases_sorted), len(axes)):
        axes[idx].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    unique = {}
    for h, l in zip(handles, labels):
        if l not in unique:
            unique[l] = h
    if len(unique) <= 10:
        legend = axes[0].legend(unique.values(), unique.keys(),
                                fontsize=5, loc="upper left",
                                bbox_to_anchor=(1.02, 1.0),
                                frameon=True, fancybox=True, framealpha=0.9,
                                title="Embedding", title_fontsize=5.5)
        legend.get_frame().set_linewidth(0.3)

    fig.suptitle("Batch mixing vs. cell-type purity across atlases and embeddings",
                 fontsize=8.5, fontweight="bold")
    _save(fig, "Fig2_ilisi_clisi_scatter")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Runtime vs cells: scalability (log-log)
# ═══════════════════════════════════════════════════════════════════════════

def fig3_scalability(data: dict):
    """Log-log plot: total runtime vs number of cells."""
    atlas_summary = data["atlas_summary"].copy()
    atlas_meta = data["atlas_meta"]

    fig, ax = plt.subplots(figsize=(5.5, 4.2))

    for _, row in atlas_summary.iterrows():
        atlas = row["atlas"]
        cells = row["cells"]
        time_s = row["total_time_s"]
        has_cluster = row["has_cluster"]
        is_pca = cells > 300000

        if time_s <= 0:
            continue

        marker = "s" if is_pca else "o"
        edgecolor = "#d6604d" if not has_cluster else "#2166ac"
        facecolor = edgecolor
        size = 60

        ax.scatter(cells, time_s, c=facecolor, s=size, marker=marker,
                  edgecolors="white", linewidth=0.5, zorder=3, alpha=0.9)

        offset_x = 1.15
        offset_y = 1.05
        ax.annotate(atlas, (cells, time_s),
                    xytext=(cells * offset_x, time_s * offset_y),
                    fontsize=5.5, color="#333333",
                    ha="left", va="bottom",
                    arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.5))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of cells", fontsize=7.5)
    ax.set_ylabel("Total runtime (seconds)", fontsize=7.5)

    # Reference lines
    x_ref = np.logspace(3, 7, 50)
    ax.plot(x_ref, x_ref * 0.002, "k--", linewidth=0.5, alpha=0.4, label="O(N)")
    ax.plot(x_ref, x_ref * np.log(x_ref) * 0.0001, "k:", linewidth=0.5, alpha=0.4, label="O(N log N)")

    ax.axvline(x=300000, color="#e08214", linewidth=0.6, linestyle="--", alpha=0.6)
    ax.text(310000, 0.5, "PCA subsampling\nthreshold", fontsize=5.5, color="#e08214",
            rotation=90, va="bottom")

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="#2166ac", edgecolor="white", label="Full-feature atlas (with cluster)"),
        mpatches.Patch(facecolor="#d6604d", edgecolor="white", label="No cluster labels detected"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=6, label="Full gene space"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="grey", markersize=6, label="PCA reference (≥300k cells)"),
    ]
    ax.legend(handles=legend_elements, fontsize=5.5, loc="upper left", framealpha=0.9)

    ax.set_title("CheckAtlas scalability: runtime vs. atlas size", fontsize=8.5, fontweight="bold")
    ax.grid(True, alpha=0.2, linewidth=0.3)

    _save(fig, "Fig3_scalability")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Dimred metrics overview heatmap
# ═══════════════════════════════════════════════════════════════════════════

def fig4_dimred_heatmap(data: dict):
    """Heatmap of dimred metrics: rows=atlas×embedding, columns=metric."""
    df = data["dimred"].copy()

    df["atlas_embedding"] = df["atlas"] + " | " + df["embedding"].apply(
        lambda x: x.replace("_scvi_full_donorassay", "_scvi_full")
        .replace("_tissue_scvi_donorassay", "_tissue_scvi")
        .replace("_compartment_scvi_donorassay", "_comp_scvi")
        .replace("_tissue_uncorrected_umap", "_unco_umap")
        .replace("_compartment_uncorrected_umap", "_comp_unco")
        .replace("_uncorrected_umap", "_unco")
        .replace("tissue_uncorrected_umap", "X_unco_umap")
        .replace("_umap_scvi_full_donorassay", "_scvi_full")
        .replace("_umap_compartment_scvi_donorassay", "_comp_scvi")
    )

    metrics_ordered = ["spearman_rho", "dCor", "coknn", "entourage", "lcmc", "den_pre", "kruskal_stress", "avg_jaccard_dis", "ged"]
    metric_labels = ["Spearman ρ", "dCor", "coKNN", "Entourage", "LCMC", "DenPre", "Kruskal\nStress", "Jaccard\nDist", "GED"]

    rows = sorted(df["atlas_embedding"].unique())

    heatmap_data = np.full((len(rows), len(metrics_ordered)), np.nan)
    for i, ae in enumerate(rows):
        atlas = ae.split(" | ")[0]
        emb_part = " | ".join(ae.split(" | ")[1:])
        for j, m in enumerate(metrics_ordered):
            subset = df[(df["atlas"] == atlas.split(" | ")[0]) & (df["atlas_embedding"] == ae) & (df["metric"] == m)]
            if len(subset) > 0 and "scaled" in subset.columns:
                heatmap_data[i, j] = subset["scaled"].mean()

    fig, ax = plt.subplots(figsize=(9, max(5, len(rows) * 0.25)), constrained_layout=True)

    cmap = sns.color_palette("RdYlGn", as_cmap=True)
    im = ax.imshow(heatmap_data, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")

    # Overlay NaN
    for i in range(len(rows)):
        for j in range(len(metrics_ordered)):
            if np.isnan(heatmap_data[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=True, facecolor="#e0e0e0", edgecolor="#cccccc",
                                           linewidth=0.3, hatch="///"))

    ax.set_xticks(range(len(metrics_ordered)))
    ax.set_xticklabels(metric_labels, fontsize=5.5, rotation=45, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=5.5)

    # Group by atlas
    prev_atlas = None
    for i, ae in enumerate(rows):
        curr_atlas = ae.split(" | ")[0]
        if curr_atlas != prev_atlas:
            ax.axhline(y=i - 0.5, color="white", linewidth=1.2)
        prev_atlas = curr_atlas

    # Color bar
    cbar = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.02, shrink=0.8)
    cbar.set_label("Scaled score (0=worst, 1=best)", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # Annotate global vs local boundary
    ax.axvline(x=1.5, color="white", linewidth=1.5)
    ax.text(0.75, -1.3, "Global structure", fontsize=6, ha="center", fontweight="bold", color="#2166ac")
    ax.text(5.5, -1.3, "Local neighborhood", fontsize=6, ha="center", fontweight="bold", color="#b2182b")

    ax.set_title("Dimensionality reduction quality: embedding fidelity across atlases",
                 fontsize=8, fontweight="bold", pad=10)

    _save(fig, "Fig4_dimred_heatmap")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Per-task score distributions across atlases
# ═══════════════════════════════════════════════════════════════════════════

def fig5_task_scores(data: dict):
    """Bar chart: per-task scores per atlas with overall score."""
    atlas_summary = data["atlas_summary"].copy()

    atlases = [a for a in ATLAS_LABELS if a in atlas_summary["atlas"].values]
    tasks = ["cluster_score", "annotation_score", "batch_score", "dimred_score"]
    task_labels = ["Cluster", "Annotation", "Batch Corr.", "Dimred"]
    task_colors = [C_CLUSTER, C_BIO, C_BATCH, C_DIMRED]

    x = np.arange(len(atlases))
    width = 0.18
    fig, ax = plt.subplots(figsize=(11, 4.5))

    for i, (task, label, color) in enumerate(zip(tasks, task_labels, task_colors)):
        scores = []
        for atlas in atlases:
            row = atlas_summary[atlas_summary["atlas"] == atlas]
            scores.append(row[task].values[0] if len(row) > 0 else 0)
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, scores, width, label=label, color=color, alpha=0.85,
                      edgecolor="white", linewidth=0.3)

    # Overall score as black dots
    overall_scores = []
    for atlas in atlases:
        row = atlas_summary[atlas_summary["atlas"] == atlas]
        overall_scores.append(row["overall_score"].values[0] if len(row) > 0 else 0)

    ax.scatter(x, overall_scores, c=C_OVERALL, s=50, zorder=5, marker="D",
              edgecolors="white", linewidth=0.5, label="Overall")

    ax.set_xticks(x)
    ax.set_xticklabels([_fmt_atlas(a) for a in atlases], fontsize=6.5)
    ax.set_ylabel("Scaled score (0–1)", fontsize=7)
    ax.set_ylim(0, 0.65)
    ax.set_title("Per-task and overall atlas quality scores", fontsize=8.5, fontweight="bold")

    ax.legend(fontsize=5.5, loc="upper right", framealpha=0.9, ncol=5,
             columnspacing=0.8, handletextpad=0.5)
    ax.grid(axis="y", alpha=0.2, linewidth=0.3)

    _save(fig, "Fig5_task_scores")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE S1 — Silhouette and Davies-Bouldin per embedding per atlas
# ═══════════════════════════════════════════════════════════════════════════

def fig_s1_cluster_details(data: dict):
    """Multi-panel: Silhouette width per embedding per atlas (cluster quality)."""
    df = data["cluster"].copy()
    atlases_with_cluster = sorted(df["atlas"].unique())
    n_col = 3
    n_row = int(np.ceil(len(atlases_with_cluster) / n_col))

    fig, axes = plt.subplots(n_row, n_col, figsize=(3.5 * n_col, 2.5 * n_row),
                            squeeze=False, constrained_layout=True)
    axes = axes.flatten()

    for idx, atlas in enumerate(atlases_with_cluster):
        ax = axes[idx]
        sil_data = df[(df["atlas"] == atlas) & (df["metric"] == "silhouette")].dropna(subset=["value"]).copy()

        if sil_data.empty:
            ax.set_visible(False)
            continue

        if "embedding_clean" not in sil_data.columns:
            sil_data["embedding_clean"] = sil_data["embedding"].apply(
                lambda x: x.replace("_scvi_full_donorassay", "_scvi_f").replace("_tissue_scvi_donorassay", "_t_scvi")
                .replace("_compartment_scvi_donorassay", "_c_scvi").replace("_uncorrected_umap", "_unco")
                .replace("tissue_uncorrected_umap", "X_tissue_unco").replace("_umap_compartment_scvi_donorassay", "_c_scvi")
            )

        emb_avg = sil_data.groupby("embedding_clean")["value"].mean().sort_values()

        colors = []
        for e in emb_avg.index:
            for orig, short in zip(sil_data["embedding"], sil_data["embedding_clean"]):
                if short == e:
                    colors.append(EMBEDDING_COLORS.get(orig, "#999999"))
                    break
            else:
                colors.append("#999999")

        bars = ax.barh(range(len(emb_avg)), emb_avg.values, color=colors,
                       edgecolor="white", linewidth=0.4, height=0.65)
        ax.axvline(x=0, color="black", linewidth=0.5)

        ax.set_yticks(range(len(emb_avg)))
        ax.set_yticklabels(emb_avg.index, fontsize=4.8)
        ax.set_xlabel("Silhouette width", fontsize=6)
        ax.set_title(f"{_fmt_atlas(atlas)}", fontsize=7.5, fontweight="bold")
        ax.set_xlim(max(-0.5, emb_avg.min() - 0.1), min(1.0, emb_avg.max() + 0.15))

        for i, (v, e) in enumerate(zip(emb_avg.values, emb_avg.index)):
            if v > 0:
                ax.text(v + 0.01, i, f"{v:.2f}", fontsize=4.5, va="center")
            else:
                ax.text(v - 0.08, i, f"{v:.2f}", fontsize=4.5, va="center", ha="right")

    for idx in range(len(atlases_with_cluster), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Cluster quality: mean silhouette width per embedding per atlas",
                fontsize=9, fontweight="bold")
    _save(fig, "FigS1_cluster_details")
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE S2 — Runtime breakdown stacked bar
# ═══════════════════════════════════════════════════════════════════════════

def fig_s2_runtime_breakdown(data: dict):
    """Stacked bar: runtime breakdown per atlas from log data."""
    # Hardcoded from log analysis (preprocess times were provided per atlas)
    runtime_data = {
        "TS_neural": {"Dimred prep": 22, "Other prep": 7, "Cluster": 14, "Annot": 0, "Batch": 11, "Dimred": 13},
        "TS_skin": {"Dimred prep": 189, "Other prep": 30, "Cluster": 121, "Annot": 44, "Batch": 22, "Dimred": 33},
        "bone_marrow": {"Dimred prep": 490, "Other prep": 48, "Cluster": 340, "Annot": 253, "Batch": 44, "Dimred": 64},
        "liver": {"Dimred prep": 346, "Other prep": 42, "Cluster": 251, "Annot": 148, "Batch": 34, "Dimred": 57},
        "blood_pbmc": {"Dimred prep": 649, "Other prep": 32, "Cluster": 0, "Annot": 7, "Batch": 20, "Dimred": 26},
        "blood": {"Dimred prep": 3760, "Other prep": 51, "Cluster": 994, "Annot": 76, "Batch": 68, "Dimred": 81},
        "lung": {"Dimred prep": 677, "Other prep": 448, "Cluster": 0, "Annot": 13, "Batch": 112, "Dimred": 22},
        "TS_immune": {"Dimred prep": 2285, "Other prep": 1353, "Cluster": 741, "Annot": 89, "Batch": 649, "Dimred": 142},
        "lung_2M": {"Dimred prep": 3213, "Other prep": 2412, "Cluster": 0, "Annot": 30, "Batch": 1413, "Dimred": 149},
    }

    atlases = [a for a in ATLAS_LABELS if a in runtime_data]
    categories = ["Dimred prep", "Other prep", "Cluster", "Annot", "Batch", "Dimred"]
    cat_colors = ["#4393c3", "#92c5de", C_CLUSTER, C_BIO, C_BATCH, C_DIMRED]

    fig, ax = plt.subplots(figsize=(10, 4.5))

    x = np.arange(len(atlases))
    bottom = np.zeros(len(atlases))

    for i, (cat, color) in enumerate(zip(categories, cat_colors)):
        vals = [runtime_data[a].get(cat, 0) for a in atlases]
        ax.bar(x, vals, 0.65, bottom=bottom, label=cat, color=color,
               edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([_fmt_atlas(a) for a in atlases], fontsize=6.5)
    ax.set_ylabel("Runtime (seconds)", fontsize=7)
    ax.set_title("Runtime breakdown per atlas per task", fontsize=8.5, fontweight="bold")
    ax.legend(fontsize=5.5, loc="upper left", framealpha=0.9, ncol=3)

    # Add total time annotations
    for i, atlas in enumerate(atlases):
        total = sum(runtime_data[atlas].values())
        ax.text(i, bottom[i] + max(bottom) * 0.02, f"{total:,}s",
                fontsize=5, ha="center", va="bottom", fontweight="bold")

    ax.grid(axis="y", alpha=0.15, linewidth=0.3)

    _save(fig, "FigS2_runtime_breakdown")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE S3 — Metric correlation matrix
# ═══════════════════════════════════════════════════════════════════════════

def fig_s3_metric_correlations(data: dict):
    """Inter-metric correlation heatmap within each task."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)

    task_configs = [
        ("cluster", data["cluster"], ["silhouette", "calinski_harabasz", "davies_bouldin", "dbcvi", "kolmogorov_smirnov"]),
        ("dimred", data["dimred"], ["spearman_rho", "dCor", "coknn", "entourage", "lcmc", "den_pre", "kruskal_stress", "avg_jaccard_dis", "ged"]),
        ("batch_correction", data["batch_tech"], ["iLISI", "cLISI", "kbet", "pcr", "graph_connectivity"]),
    ]

    for ax_idx, (task_name, df, metrics) in enumerate(task_configs):
        ax = axes[ax_idx]

        # Build correlation matrix
        df_pivot = df.pivot_table(values="value", index=["atlas", "embedding"], columns="metric", aggfunc="mean").dropna(axis=1)
        available = [m for m in metrics if m in df_pivot.columns]
        if len(available) < 2:
            ax.text(0.5, 0.5, f"No {task_name} data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(task_name, fontsize=8, fontweight="bold")
            continue

        corr = df_pivot[available].corr()

        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        cmap = sns.diverging_palette(240, 10, s=75, l=45, center="light", as_cmap=True)

        im = ax.imshow(corr.values, cmap=cmap, vmin=-1, vmax=1, aspect="auto")

        # Annotate values
        for i in range(len(available)):
            for j in range(len(available)):
                if i != j:
                    val = corr.values[i, j]
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                           fontsize=5.5, color="black" if abs(val) < 0.6 else "white",
                           fontweight="bold" if abs(val) > 0.7 else "normal")

        ax.set_xticks(range(len(available)))
        ax.set_xticklabels(available, fontsize=5.5, rotation=45, ha="right")
        ax.set_yticks(range(len(available)))
        ax.set_yticklabels(available, fontsize=5.5)
        ax.set_title(f"{task_name} metrics", fontsize=7.5, fontweight="bold")

    cbar = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, shrink=0.8)
    cbar.set_label("Pearson r", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.5)

    fig.suptitle("Inter-metric correlations within each task", fontsize=9, fontweight="bold")
    _save(fig, "FigS3_metric_correlations")

    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE S4 — Isolated F1 score analysis
# ═══════════════════════════════════════════════════════════════════════════

def fig_s4_isolated_f1(data: dict):
    """Box/strip plot of isolated F1 scores across atlases."""
    df = data["annot"].copy()
    df_f1 = df[df["metric"] == "isolated_f1_score"].dropna(subset=["value"])

    if df_f1.empty:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "No isolated F1 data available", ha="center", va="center", transform=ax.transAxes)
        _save(fig, "FigS4_isolated_f1")
        return fig

    atlases = sorted(df_f1["atlas"].unique())

    fig, ax = plt.subplots(figsize=(max(6, len(atlases) * 0.9), 3.5))

    positions = range(len(atlases))
    for i, atlas in enumerate(atlases):
        vals = df_f1[df_f1["atlas"] == atlas]["value"].dropna()
        if len(vals) > 0:
            parts = ax.violinplot(vals, positions=[i], showmeans=True,
                                  showmedians=True, widths=0.6)
            for body in parts["bodies"]:
                body.set_facecolor(C_BIO)
                body.set_alpha(0.5)

    ax.set_xticks(positions)
    ax.set_xticklabels([_fmt_atlas(a) for a in atlases], fontsize=6.5)
    ax.set_ylabel("Isolated F1 score", fontsize=7)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=0.5, color="#999999", linewidth=0.4, linestyle="--")
    ax.text(len(atlases) - 0.5, 0.51, "scIB threshold (0.50)", fontsize=5.5,
           color="#666666", ha="right", va="bottom")
    ax.set_title("Rare cell-type preservation: Isolated F1 scores per atlas",
                fontsize=8, fontweight="bold")
    ax.grid(axis="y", alpha=0.15, linewidth=0.3)

    _save(fig, "FigS4_isolated_f1")
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE S5 — PCR analysis: batch predictability
# ═══════════════════════════════════════════════════════════════════════════

def fig_s5_pcr_analysis(data: dict):
    """Per-atlas PCR scores: how well can batch be predicted from embedding?"""
    df = data["batch_tech"].copy()
    df_pcr = df[df["metric"] == "pcr"].dropna(subset=["value"])

    if df_pcr.empty:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "No PCR data available", ha="center", va="center", transform=ax.transAxes)
        _save(fig, "FigS5_pcr_analysis")
        return fig

    atlases = sorted(df_pcr["atlas"].unique(), key=lambda a: ATLAS_LABELS.get(a, a))

    top_keys = df_pcr.groupby("batch_key")["value"].mean().sort_values(ascending=False).head(8).index.tolist()

    n_rows = int(np.ceil(len(atlases) / 3))
    fig, axes = plt.subplots(n_rows, 3, figsize=(12, 3 * n_rows),
                            squeeze=False, constrained_layout=True)
    axes = axes.flatten()

    for idx, atlas in enumerate(atlases):
        ax = axes[idx]
        atlas_data = df_pcr[df_pcr["atlas"] == atlas].copy()

        if "embedding_clean" not in atlas_data.columns:
            atlas_data["embedding_clean"] = atlas_data["embedding"]

        pca_data = atlas_data[atlas_data["embedding_clean"].str.contains("pca", case=False)]
        scvi_data = atlas_data[atlas_data["embedding_clean"].str.contains("scvi", case=False)]

        common_keys = set(pca_data["batch_key"]) & set(scvi_data["batch_key"])
        if not common_keys:
            ax.text(0.5, 0.5, f"No comparable\nbatch keys", ha="center", va="center",
                   fontsize=8, transform=ax.transAxes)
            ax.set_title(_fmt_atlas(atlas), fontsize=7.5, fontweight="bold")
            continue

        common_keys = sorted(common_keys)[:5]
        x = np.arange(len(common_keys))
        width = 0.35

        pca_vals = [pca_data[pca_data["batch_key"] == k]["value"].mean() for k in common_keys]
        scvi_vals = [scvi_data[scvi_data["batch_key"] == k]["value"].mean() for k in common_keys]

        ax.bar(x - width / 2, pca_vals, width, label="PCA", color="#2166ac", alpha=0.8,
              edgecolor="white", linewidth=0.3)
        ax.bar(x + width / 2, scvi_vals, width, label="scVI", color="#b2182b", alpha=0.8,
              edgecolor="white", linewidth=0.3)

        ax.set_xticks(x)
        ax.set_xticklabels(common_keys, fontsize=4.5, rotation=45, ha="right")
        ax.set_ylabel("PCR (lower=better)", fontsize=6)
        ax.set_ylim(0, 1.05)
        ax.set_title(_fmt_atlas(atlas), fontsize=7, fontweight="bold")

        if idx == 0:
            ax.legend(fontsize=6, loc="upper right")

    for idx in range(len(atlases), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Principal Component Regression: batch predictability before and after scVI integration",
                fontsize=9, fontweight="bold", y=1.01)
    fig.tight_layout()
    _save(fig, "FigS5_pcr_analysis")
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Master: generate all figures
# ═══════════════════════════════════════════════════════════════════════════

def generate_all_figures(data: dict):
    """Generate all main and supplementary figures."""
    print("Generating Fig 1: Overview heatmap...")
    fig1_overview_heatmap(data)

    print("Generating Fig 2: iLISI vs cLISI scatter...")
    fig2_ilisi_vs_clisi(data)

    print("Generating Fig 3: Scalability...")
    fig3_scalability(data)

    print("Generating Fig 4: Dimred heatmap...")
    fig4_dimred_heatmap(data)

    print("Generating Fig 5: Task scores...")
    fig5_task_scores(data)

    print("Generating Fig S1: Cluster details...")
    fig_s1_cluster_details(data)

    print("Generating Fig S2: Runtime breakdown...")
    fig_s2_runtime_breakdown(data)

    print("Generating Fig S3: Metric correlations...")
    fig_s3_metric_correlations(data)

    print("Generating Fig S4: Isolated F1...")
    fig_s4_isolated_f1(data)

    print("Generating Fig S5: PCR analysis...")
    fig_s5_pcr_analysis(data)

    print(f"\nAll figures saved to: {OUT_DIR}")
