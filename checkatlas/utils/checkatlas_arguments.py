import argparse
from importlib.resources import files

try:
    from .. import atlas
    from ..metrics import annot, cluster, dimred
except ImportError:
    from checkatlas import atlas
    from checkatlas.metrics import annot, cluster, dimred


MAX_N_JOBS = 48
"""Upper limit on the number of CPU threads any CheckAtlas process may
consume.  Capped at 48 so that on a 80-core workstation at least 32
threads remain free for other users / pipelines.  See
:func:`cap_n_jobs` for the enforcement helper."""


def cap_n_jobs(n_jobs):
    """Cap *n_jobs* at :data:`MAX_N_JOBS` and treat ``None`` / ``-1``
    as "use the system default" by resolving to the cap.

    Rules:
      * ``None``  → :data:`MAX_N_JOBS` (covers the legacy path where
        ``getattr(args, "n_jobs", -1)`` returned the fallback).
      * ``-1``    → :data:`MAX_N_JOBS` (sklearn / joblib convention:
        "all available cores", capped for politeness).
      * values > :data:`MAX_N_JOBS` → :data:`MAX_N_JOBS`.
      * any positive integer ≤ :data:`MAX_N_JOBS` → unchanged.
      * anything else (non-integer, zero, negative) → :data:`MAX_N_JOBS`.
    """
    if n_jobs is None or n_jobs == -1:
        return MAX_N_JOBS
    try:
        n_jobs = int(n_jobs)
    except (TypeError, ValueError):
        return MAX_N_JOBS
    if n_jobs > MAX_N_JOBS:
        return MAX_N_JOBS
    if n_jobs < 1:
        return MAX_N_JOBS
    return n_jobs


def create_parser():
    parser = argparse.ArgumentParser(
        prog="checkatlas",
        usage="checkatlas [OPTIONS] process path/ [--atlas_name ATLAS_NAME]",
        formatter_class=argparse.RawTextHelpFormatter,
        description="CheckAtlas is a one liner tool to check the "
        "quality of your single-cell atlases. For "
        "every atlas, it produces the quality control "
        "tables and figures which can be then processed "
        "by multiqc. CheckAtlas is able to load Scanpy, "
        "Seurat, and CellRanger files.",
        epilog="Enjoy the checkatlas functionality!",
    )

    # All Program arguments
    # All Program arguments
    # main_options = parser.add_argument_group("Main arguments")
    parser.add_argument(
        "process",
        type=str,
        help=(
            "Process to run. Choose from:\n"
            "  run            — launch the full Nextflow pipeline on a folder\n"
            "  preprocess     — load and preprocess a single atlas\n"
            "  summary        — generate summary tables and UMAP/t-SNE figures\n"
            "  qc             — compute QC tables and plots\n"
            "  metric_cluster — compute clustering metrics\n"
            "  metric_annot   — compute annotation metrics\n"
            "  metric_dimred  — compute dimensionality-reduction metrics\n"
            "  metric         — run all metric processes on a single atlas\n"
            "  analyse        — run all analysis processes on a single atlas\n"
            "  scfm           — run scFM quality-control (9 problems, FMF, BF, PR)\n"
        ),
        default="",
    )

    parser.add_argument(
        "path",
        type=str,
        help="Required argument: Your folder containing "
        "Scanpy, CellRanger and Seurat atlases",
        default=".",
    )

    parser.add_argument(
        "--atlas_name",
        type=str,
        help="Required for all processes except 'run': name of the atlas to "
        "process. Atlas_name should be found in one of the samplesheet "
        "provided to nf-checkatlas, or directly created by "
        "checkatlas.list_all_atlases() function",
        default=None,
    )

    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print out all debug messages.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Checkatlas, version {get_version()}",
        help="Display checkatlas version.",
    )

    # Arguments linked to QC
    qc_options = parser.add_argument_group("QC options")
    qc_options.add_argument(
        "--qc_display",
        nargs="+",
        type=str,
        default=[
            "violin_plot",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ],
        help="List of QC to display. "
        "Available qc = violin_plot, total_counts, "
        "n_genes_by_counts, pct_counts_mt. "
        "Default: --qc_display violin_plot total_counts "
        "n_genes_by_counts pct_counts_mt",
    )
    qc_options.add_argument(
        "--plot_celllimit",
        type=int,
        default=0,
        help="Set the maximum number of cells"
        "to plot in QC, UMAP, t-SNE, etc...."
        "If plot_celllimit=0, no limit will"
        "be applied.",
    )
    qc_options.add_argument(
        "--n_jobs",
        type=int,
        default=MAX_N_JOBS,
        help="Maximum number of CPU threads to use for parallel jobs "
        "(sklearn pairwise_distances, NearestNeighbors, kNN graphs, "
        "metric dispatch). Capped at 48 so other users / pipelines "
        "keep some cores free on a 80-core host. Pass any positive "
        "integer <= 48 to lower the cap further; higher values are "
        f"silently clamped to {MAX_N_JOBS}.",
    )

    # Arguments linked to metric
    metric_options = parser.add_argument_group("Metric options")
    metric_options.add_argument(
        "--obs_cluster",
        nargs="+",
        type=str,
        default=atlas.OBS_CLUSTERS,
        help="List of obs from the adata file to "
        "use in the clustering metric calculus."
        "Example: --obs_cluster celltype leuven seurat_clusters",
    )
    metric_options.add_argument(
        "--metric_cluster",
        nargs="+",
        type=str,
        default=cluster.__all__,
        help="List of clustering metrics to calculate.\n"
        "   By default all available clustering metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_cluster silhouette\n"
        "   Example: --metric_cluster silhouette davies_bouldin\n"
        "   Example: --metric_cluster none\n"
        f"   Available: {cluster.__all__}",
    )
    metric_options.add_argument(
        "--metric_annot",
        nargs="+",
        type=str,
        default=annot.__all__,
        help="List of annotation metrics to calculate.\n"
        "   By default all available annotation metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_annot rand_index\n"
        "   Example: --metric_annot rand_index fowlkes_mallows\n"
        "   Example: --metric_annot none\n"
        f"   Available: {annot.__all__}",
    )
    metric_options.add_argument(
        "--metric_dimred",
        nargs="+",
        type=str,
        default=dimred.__all__,
        help="List of dimensionality reduction metrics to calculate.\n"
        "   By default all available dimensionality reduction metrics are run.\n"
        "   Specify one or more metrics, or 'none' to skip this category.\n"
        "   Example: --metric_dimred kruskal_stress\n"
        "   Example: --metric_dimred kruskal_stress spearman_rho\n"
        "   Example: --metric_dimred none\n"
        f"   Available: {dimred.__all__}",
    )

    # Nextflow reporting options (only used with 'run' process)
    nf_options = parser.add_argument_group(
        "Nextflow reporting options (only used with 'run' process)"
    )
    nf_options.add_argument(
        "--with_report",
        nargs="?",
        const="Checkatlas_report.html",
        default="Checkatlas_report.html",
        metavar="FILE",
        help="Generate a Nextflow execution report (enabled by default). "
        "Optionally specify the output filename (default: Checkatlas_report.html).",
    )
    nf_options.add_argument(
        "--with_dag",
        nargs="?",
        const="dag.html",
        default=None,
        metavar="FILE",
        help="Generate a Nextflow pipeline DAG. "
        "Optionally specify the output filename (default: dag.html). "
        "Use .dot extension for Graphviz format.",
    )
    nf_options.add_argument(
        "--with_timeline",
        nargs="?",
        const="Checkatlas_timeline.html",
        default="Checkatlas_timeline.html",
        metavar="FILE",
        help="Generate a Nextflow execution timeline (enabled by default). "
        "Optionally specify the output filename (default: Checkatlas_timeline.html).",
    )

    # scFM QC options (only used with 'scfm' process)
    scfm_options = parser.add_argument_group(
        "scFM QC options (only used with 'scfm' process). "
        "Default behaviour with no --scfm_* flag: 'all combinations' "
        "mode — the orchestrator runs all 3 per-task engines and the "
        "diagnostic engine evaluates every (ref, pred, batch, scfm, "
        "baseline) combination the column detector found. Use "
        "--scfm_fast to opt out and run only the first combo."
    )
    scfm_options.add_argument(
        "--scfm_embedding",
        type=str,
        default="",
        help="Embedding key for the scFM under test "
        "(e.g. X_geneformer). Optional: scFM embeddings are also "
        "auto-detected from adata.obsm via the column detector.",
    )
    scfm_options.add_argument(
        "--baseline_embeddings",
        nargs="+",
        type=str,
        default=[],
        help="Baseline embedding keys for scFM comparison "
        "(e.g. --baseline_embeddings X_pca X_scvi).",
    )
    scfm_options.add_argument(
        "--scfm_predicted_label",
        type=str,
        default=None,
        help="obs column with the scFM's predicted cell types.",
    )
    scfm_options.add_argument(
        "--scfm_ref_label",
        type=str,
        default=None,
        help="obs column with the author/curated reference labels.",
    )
    scfm_options.add_argument(
        "--scfm_batch_key",
        type=str,
        default=None,
        help="obs column with batch / donor labels (for iLISI, kBET, PCR).",
    )
    scfm_options.add_argument(
        "--scfm_domain_key",
        type=str,
        default=None,
        help="obs column with domain labels (species, tissue, assay) "
        "for cross-domain generalisation metric (Problem 6).",
    )
    scfm_options.add_argument(
        "--scfm_patient_key",
        type=str,
        default=None,
        help="obs column with patient IDs (for Problem 7).",
    )
    scfm_options.add_argument(
        "--scfm_outcome_key",
        type=str,
        default=None,
        help="obs column with patient outcome (for Problem 7 outcome AUC).",
    )
    scfm_options.add_argument(
        "--scfm_scaling_fractions",
        nargs="+",
        type=float,
        default=[0.01, 0.10, 0.50, 1.00],
        help="Subsample fractions for the scaling-law saturation metric. "
        "Default: 0.01 0.10 0.50 1.00.",
    )
    scfm_options.add_argument(
        "--scfm_n_seeds",
        type=int,
        default=5,
        help="Number of subsample seeds for the stability metric.",
    )
    scfm_options.add_argument(
        "--scfm_noise_sigma",
        type=float,
        default=0.10,
        help="Noise sigma for the gaussian-robustness metric (reserved).",
    )
    scfm_options.add_argument(
        "--scfm_min_domain_size",
        type=int,
        default=50,
        help="Minimum number of cells per domain for cross-domain evaluation.",
    )
    scfm_options.add_argument(
        "--scfm_weights",
        type=str,
        default=None,
        help="Optional JSON file overriding composite-score weights.",
    )
    scfm_options.add_argument(
        "--scfm_thresholds",
        type=str,
        default=None,
        help="Optional YAML file overriding the diagnostic thresholds.",
    )
    scfm_options.add_argument(
        "--scfm_fast",
        action="store_true",
        help=(
            "When no --scfm_* flag is set, run the diagnostic "
            "engine on only the first (ref, pred, batch, scfm, "
            "baseline) combination the column detector finds, "
            "instead of the default 'all combinations' sweep. "
            "Faster, but may miss real signal if the auto-detected "
            "'best' combo is wrong. Default: off (run all "
            "combinations)."
        ),
    )
    scfm_options.add_argument(
        "--scfm_max_combos",
        type=int,
        default=27,
        help=(
            "Maximum number of (ref, pred, batch, scfm, baseline) "
            "combinations to evaluate when running in 'all "
            "combinations' mode. Default: 27 (3 ref * 1 pred * "
            "1 batch * 3 scfm * 3 baseline). Capped to keep "
            "compute time bounded on atlases with many columns."
        ),
    )
    return parser


def get_version():
    """
    Get version of checkatlas from checkatlas/VERSION file
    :return: checkatlas version
    """
    version_file = files(__package__).joinpath("VERSION")
    with open(version_file) as file:
        version = file.readline()
        return version
