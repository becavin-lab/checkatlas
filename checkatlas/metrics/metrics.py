import logging
import os
import time

import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    import rpy2.robjects as ro
    import rpy2.robjects as robjects
except (ImportError, OSError):
    ro = None
    robjects = None
from sklearn.preprocessing import LabelEncoder

from . import annot, batch_correction, cluster, dimred
from ._cache import (
    compute_fingerprint,
    load_dimred_cache,
    load_knn,
    load_triangular,
    save_dimred_cache,
    save_knn,
    save_triangular,
)
from ._jax_utils import _GPU_AVAILABLE, _JAX_AVAILABLE, _get_ndarray, pdist_squareform
from ._neighbors import NeighborResults, _clear_neighbors_cache, compute_neighbors
from ._triangular import TriangularMatrix, store_upper_triangle

METRICS_CLUST = cluster.__all__
METRICS_ANNOT = annot.__all__
METRICS_BATCH = batch_correction.__all__
METRICS_DIMRED = dimred.__all__


def _pivot_annot_to_wide(df, atlas_name):
    """Pivot long-format annotation results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Atlas Name, Metric Name,
    Reference/Input 1, Prediction/Input 2, Value, Time (s)] into a wide
    table headed by [Annot_Sample, Reference, obs, <metric>_value,
    <metric>_running_time, ...] where each row is one unique
    (Reference/Input 1 × Prediction/Input 2) combination.
    """
    if df.empty:
        return df
    metric_names = sorted(df["Metric Name"].unique())
    wide_rows = []
    for (ref_val, pred_val), group in df.groupby(
        ["Reference/Input 1", "Prediction/Input 2"]
    ):
        row = {
            "Annot_Sample": f"{atlas_name}_{ref_val}_{pred_val}",
            "Reference": ref_val,
            "obs": pred_val,
        }
        for _, r in group.iterrows():
            metric = r["Metric Name"]
            row[metric] = r["Value"]
            row[f"{metric}_running_time"] = r["Time (s)"]
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ["Annot_Sample", "Reference", "obs"] + [
        c for c in value_cols if c in result_df.columns
    ]
    return result_df[col_order]


def _pivot_cluster_to_wide(df, atlas_name):
    """Pivot long-format cluster results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Atlas Name, Metric Name,
    Embedding, Label Key, Value, Time (s)] into a wide table headed by
    [Clust_Sample, obs, <metric>, <metric>_running_time, ...] where each
    row is one unique (Embedding × Label Key) combination.
    """
    if df.empty:
        return df
    metric_names = sorted(df["Metric Name"].unique())
    wide_rows = []
    for (emb_val, label_val), group in df.groupby(["Embedding", "Label Key"]):
        row = {
            "Clust_Sample": f"{atlas_name}_{emb_val}_{label_val}",
            "obs": label_val,
        }
        for _, r in group.iterrows():
            metric = r["Metric Name"]
            row[metric] = r["Value"]
            row[f"{metric}_running_time"] = r["Time (s)"]
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ["Clust_Sample", "obs"] + [
        c for c in value_cols if c in result_df.columns
    ]
    return result_df[col_order]


def _pivot_dimred_to_wide(df, atlas_name):
    """Pivot long-format dimred results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Metric Name, Low Dim Key,
    High Dim Key, Value, Time (s)] into a wide table headed by
    [Dimred_Sample, obsm, <metric>, <metric>_running_time, …] where
    each row is one unique Low Dim Key (embedding).
    """
    if df.empty:
        return df
    metric_names = sorted(df["Metric Name"].unique())
    wide_rows = []
    for low_dim_key, group in df.groupby("Low Dim Key"):
        row = {
            "Dimred_Sample": f"{atlas_name}_{low_dim_key}",
            "obsm": low_dim_key,
        }
        for _, r in group.iterrows():
            metric = r["Metric Name"]
            row[metric] = r["Value"]
            row[f"{metric}_running_time"] = r["Time (s)"]
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ["Dimred_Sample", "obsm"] + [
        c for c in value_cols if c in result_df.columns
    ]
    return result_df[col_order]


def _pivot_batch_correction_to_wide(df, atlas_name):
    """Pivot long-format batch-correction results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Atlas Name, Metric Name,
    Embedding, Batch Key, Value, Time (s)] into a wide table headed
    by [Batch_Sample, Embedding, Batch Key, <metric>, <metric>_running_time, ...]
    where each row is one unique (Embedding, Batch Key) combination.

    The ``Metric Name`` column is the long-format name (e.g. ``iLISI``,
    ``cLISI``, ``kbet``, ``pcr``, ``graph_connectivity``); each value
    gets its own column in the wide output. ``iLISI`` and ``cLISI`` are
    distinguished by which ``Batch Key`` row they live in
    (batch column vs cell-type column).
    """
    if df.empty:
        return df
    metric_names = sorted(df["Metric Name"].unique())
    wide_rows = []
    for (emb, batch), group in df.groupby(["Embedding", "Batch Key"]):
        row = {
            "Batch_Sample": f"{atlas_name}_{emb}_{batch}",
            "Embedding": emb,
            "Batch Key": batch,
        }
        for _, r in group.iterrows():
            metric = r["Metric Name"]
            row[metric] = r["Value"]
            row[f"{metric}_running_time"] = r["Time (s)"]
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ["Batch_Sample", "Embedding", "Batch Key"] + [
        c for c in value_cols if c in result_df.columns
    ]
    return result_df[col_order]


logger = logging.getLogger("checkatlas")

if robjects is not None:
    R_ANNOT = robjects.r(
        "type <- function(seurat, obs_key){ " "return(seurat[[obs_key]][[obs_key]])}"
    )
    R_REDUCTION = robjects.r(
        "reduc <- function(seurat, obsm_key){"
        " return(Embeddings(object = seurat, "
        "reduction = obsm_key))}"
    )
else:
    R_ANNOT = None
    R_REDUCTION = None


def cal_annot(
    adata,
    atlas_name=None,
    metric_list=None,
    all=False,
    file_dir=None,
    n_jobs=-1,
    verbose=True,
    preprocess_context=None,
):
    """
    Comprehensive annotation pipeline for all annotation metrics.

    Args:
        adata (AnnData): Annotated data matrix.
        metric_list (list): List of metric names to calculate.
                           If provided, overrides `all` parameter.
        all (bool): If True, calculate all available annotation metrics.
                    If False, calculate a default subset.
                    Ignored if metric_list is provided.
        file_dir (str): Directory path where the results CSV will be saved.
                       If None, saves to current working directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress information.
        preprocess_context (PreprocessContext, optional): If provided,
            column detection and kNN precomputation are skipped and
            the context's precomputed data (kNN graphs, neighbour
            graphs, batch keys) is reused.

    Returns:
        pd.DataFrame: Results dataframe with columns:
                       [Atlas Name, Metric Name, Reference/Input 1, Prediction/Input 2, Value, Time (s)]
    """
    import inspect

    from ..utils.col_detector import CheckAtlasColumnDetector

    # Set file directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)

    _subsample_indices = None

    if preprocess_context is not None:
        ref_keys = preprocess_context.ref_keys
        pred_keys = preprocess_context.pred_keys
        embedding_keys = (
            getattr(preprocess_context, "annotation_embedding_keys", None)
            or preprocess_context.embedding_keys
        )

        if not ref_keys and not pred_keys:
            logger.warning(
                "No reference or predicted annotation keys in preprocess context. Skipping."
            )
            return pd.DataFrame()

        # ── Warn when one annotation key type is missing ──────────
        _ref_pred_names = [
            "adj_mutual_info",
            "adj_rand_index",
            "fowlkes_mallows",
            "isolated_f1_score",
            "mutual_info",
            "normalized_mutual_info",
            "rand_index",
            "vmeasure",
        ]
        if not pred_keys and ref_keys:
            logger.warning(
                "No predicted annotation keys detected. "
                "Skipping ref-vs-pred metrics (%s). "
                "Running embedding-label metrics on reference keys only.",
                ", ".join(_ref_pred_names),
            )
        elif not ref_keys and pred_keys:
            logger.warning(
                "No reference annotation keys detected. "
                "Skipping ref-vs-pred metrics (%s). "
                "Running embedding-label metrics on predicted keys only.",
                ", ".join(_ref_pred_names),
            )

        if verbose:
            print("Using precomputed context — skipping column detection")
            print(f"  Reference keys: {ref_keys}")
            print(f"  Predicted keys: {pred_keys}")
            print(f"  Embedding keys: {embedding_keys}")

        # ── Load precomputed distance matrices from cluster cache ──
        # Large atlases (>300k cells) use core-set subsampling; the
        # .tri files were built on the subsample, not the full dataset.
        _subsample_indices = getattr(preprocess_context, "subsample_indices", None)
        _ss_n = len(_subsample_indices) if _subsample_indices is not None else 0
        precomputed_dists = {}
        _safe = lambda s: s.replace("/", "_").replace(" ", "_")
        for emb in embedding_keys:
            tri_path = os.path.join(
                preprocess_context.cluster_dir,
                f"dist_{_safe(emb)}.tri",
            )
            npy_path = tri_path.replace(".tri", ".npy")
            if os.path.exists(tri_path):
                if _ss_n:
                    n_cells = _ss_n
                elif emb in adata.obsm:
                    n_cells = adata.obsm[emb].shape[0]
                else:
                    n_cells = adata.n_obs
                file_n = TriangularMatrix.n_from_file(tri_path)
                if file_n == n_cells:
                    precomputed_dists[emb] = TriangularMatrix(
                        n=n_cells, filepath=tri_path, mode="r"
                    )
                else:
                    logger.warning(
                        "Skipping stale distance matrix %s "
                        "(file N=%d does not match expected N=%d)",
                        tri_path,
                        file_n,
                        n_cells,
                    )
            elif os.path.exists(npy_path):
                precomputed_dists[emb] = np.load(npy_path)
    else:
        # Detect columns
        detector = CheckAtlasColumnDetector(adata)
        params = detector.detect_all_parameters()

        ref_keys = [x[0] for x in params["annotation"]["reference"]]
        pred_keys = [x[0] for x in params["annotation"]["predicted"]]
        embedding_keys = [x[0] for x in params["clustering"]["embeddings"]]
        precomputed_dists = {}

    # ── Warn when one annotation key type is missing ──────────────
    # Ref-vs-pred metrics (8 total) need both reference and predicted
    # columns.  Embedding-label metrics (average_silhouette_width,
    # dunn_index) only need one label set per embedding.
    _ref_pred_names = [
        "adj_mutual_info",
        "adj_rand_index",
        "fowlkes_mallows",
        "isolated_f1_score",
        "mutual_info",
        "normalized_mutual_info",
        "rand_index",
        "vmeasure",
    ]
    if not pred_keys and ref_keys:
        logger.warning(
            "No predicted annotation keys detected. "
            "Skipping ref-vs-pred metrics (%s). "
            "Running embedding-label metrics on reference keys only.",
            ", ".join(_ref_pred_names),
        )
    elif not ref_keys and pred_keys:
        logger.warning(
            "No reference annotation keys detected. "
            "Skipping ref-vs-pred metrics (%s). "
            "Running embedding-label metrics on predicted keys only.",
            ", ".join(_ref_pred_names),
        )

    # Define metrics to run
    if metric_list is not None:
        metrics_list = [m for m in metric_list if m in METRICS_ANNOT]
    elif all:
        metrics_list = METRICS_ANNOT
    else:
        metrics_list = [
            "adj_rand_index",
            "normalized_mutual_info",
            "adj_mutual_info",
        ]
        metrics_list = [m for m in metrics_list if m in METRICS_ANNOT]

    results = []
    atlas_name = atlas_name

    # Categorize metrics based on their input requirements
    # Ref vs Pred
    ref_pred_metrics = [
        "adj_mutual_info",
        "adj_rand_index",
        "fowlkes_mallows",
        "isolated_f1_score",
        "mutual_info",
        "normalized_mutual_info",
        "rand_index",
        "vmeasure",
    ]

    # Embedding + Labels
    emb_label_metrics = ["average_silhouette_width", "dunn_index"]

    # Bio Conservation (adata_before, adata_after) - Skipping for single adata pipeline
    # unless we define strategy.
    bio_metrics = ["cell_cycle_conservation", "highly_variable_genes"]

    # Create progress bar with custom format
    pbar = tqdm(
        metrics_list,
        desc="Calculating Annotation Metrics",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for metric in pbar:
        # Start timing for this metric
        metric_start_time = time.time()

        # Update progress bar with current metric name
        pbar.set_description(f"Processing: {metric}")

        metric_module = getattr(annot, metric)

        try:
            # 1. Ref vs Pred Metrics
            if metric in ref_pred_metrics:
                if not ref_keys or not pred_keys:
                    continue
                for ref in ref_keys:
                    for pred in pred_keys:
                        # Skip if ref == pred
                        if ref == pred:
                            continue
                        try:
                            labels_true = adata.obs[ref]
                            labels_pred = adata.obs[pred]
                            l_pred, l_true = annotation_to_num(labels_pred, labels_true)

                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if "n_jobs" in sig.parameters:
                                kw["n_jobs"] = n_jobs
                            if "verbose" in sig.parameters:
                                kw["verbose"] = False
                            val = metric_module.run(l_pred, l_true, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append(
                                {
                                    "Atlas Name": atlas_name,
                                    "Metric Name": metric,
                                    "Reference/Input 1": ref,
                                    "Prediction/Input 2": pred,
                                    "Value": val,
                                    "Time (s)": round(pair_elapsed, 3),
                                }
                            )
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate {metric} for {ref} vs {pred}: {e}"
                            )

            # 2. Embedding + Labels (ASW, Dunn)
            elif metric in emb_label_metrics:
                if not embedding_keys:
                    continue
                targets = list(set(ref_keys + pred_keys))
                if not targets:
                    continue
                for emb in embedding_keys:
                    if emb not in adata.obsm:
                        continue
                    X_emb = adata.obsm[emb]
                    if _subsample_indices is not None:
                        X_emb = X_emb[_subsample_indices]
                    for label in targets:
                        try:
                            labels = adata.obs[label]
                            if _subsample_indices is not None:
                                labels = labels.iloc[_subsample_indices].values
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if "n_jobs" in sig.parameters:
                                kw["n_jobs"] = n_jobs
                            if "verbose" in sig.parameters:
                                kw["verbose"] = False
                            if (
                                "precomputed_dists" in sig.parameters
                                and emb in precomputed_dists
                            ):
                                kw["precomputed_dists"] = precomputed_dists[emb]
                            val = metric_module.run(X_emb, labels, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append(
                                {
                                    "Atlas Name": atlas_name,
                                    "Metric Name": metric,
                                    "Reference/Input 1": emb,
                                    "Prediction/Input 2": label,
                                    "Value": val,
                                    "Time (s)": round(pair_elapsed, 3),
                                }
                            )
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate {metric} for {emb} vs {label}: {e}"
                            )

            # 3. Bio Conservation
            elif metric in bio_metrics:
                logger.info(f"Skipping {metric} as it requires 'adata_before'.")

            else:
                logger.warning(f"Metric {metric} not categorized in pipeline.")

        except Exception as e:
            logger.error(f"Error running metric {metric}: {e}")

        metric_elapsed = time.time() - metric_start_time
        pbar.set_postfix_str(f"Time: {metric_elapsed:.2f}s", refresh=True)

    df = pd.DataFrame(results)

    # Save MultiQC-compatible wide format to file_dir if provided
    if not df.empty and file_dir is not None and atlas_name is not None:
        os.makedirs(file_dir, exist_ok=True)
        wide_df = _pivot_annot_to_wide(df, atlas_name)
        wide_path = os.path.join(file_dir, f"{atlas_name}.tsv")
        wide_df.to_csv(wide_path, sep="\t", index=False)
        logger.info("MultiQC-compatible annotation table saved to %s", wide_path)

    return df


def cal_batch_correction(
    adata,
    atlas_name=None,
    metric_list=None,
    all=False,
    file_dir=None,
    n_jobs=-1,
    verbose=True,
    preprocess_context=None,
):
    """
    Comprehensive batch-correction / integration-quality pipeline.

    Computes per-(embedding, batch-key) scores for the four
    integration metrics:

      - ``kbet``  : k-nearest-neighbour batch-effect test (rejection
        rate, lower is better)
      - ``lisi``  : Local Inverse Simpson's Index, evaluated in two
        flavours: ``iLISI`` (batch labels, higher is better) and
        ``cLISI`` (cell-type labels, lower is better)
      - ``pcr``   : principal-component regression of batch labels
        (lower is better)
      - ``graph_connectivity`` : fraction of cells in the largest
        connected component of the per-label kNN graph (higher is
        better)

    Inputs (auto-detected via :class:`CheckAtlasColumnDetector` or
    supplied through a :class:`PreprocessContext`):

      - ``batch_keys``   : obs columns with batch / donor / sample /
        covariate labels (used for kBET, iLISI, PCR)
      - ``ref_keys + pred_keys`` : ref / pred cell-type columns (used
        for cLISI and graph_connectivity)
      - ``embedding_keys``: adata.obsm keys with > 2 components

    Parameters
    ----------
    adata : AnnData
    atlas_name : str, optional
        Atlas name for the long-format ``Atlas Name`` column.
    metric_list : list of str, optional
        Subset of metrics to compute. Defaults to
        :data:`METRICS_BATCH` when ``all`` is True.
    all : bool
        If True, run every metric in :data:`METRICS_BATCH`.
    file_dir : str, optional
        Directory for the per-atlas wide-format TSV
        (``<atlas_name>.tsv``). Falls back to ``os.getcwd()``.
    n_jobs : int
    verbose : bool
    preprocess_context : PreprocessContext, optional
        When provided, column detection is skipped and the cached
        kNN graphs / neighbour graphs are reused.

    Returns
    -------
    pd.DataFrame
        Long-format table with columns
        ``[Atlas Name, Metric Name, Embedding, Batch Key,
        Value, Time (s)]``.

        ``Embedding`` is the obsm key the metric was evaluated on
        (or ``"X"`` for the raw expression matrix). ``Batch Key``
        is the obs column whose labels drove the metric (batch
        column for iLISI / kBET / PCR; cell-type column for cLISI
        and graph_connectivity).
    """
    import inspect

    from ..utils.col_detector import CheckAtlasColumnDetector
    from ._jax_utils import _GPU_AVAILABLE as _cal_gpu
    from ._jax_utils import _JAX_AVAILABLE as _cal_jax
    from ._neighbors import NeighborResults
    from ._neighbors import compute_neighbors as _cal_knn

    _USE_JAX = _cal_jax and _cal_gpu

    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)

    emb_nn: dict = {}

    if preprocess_context is not None:
        ref_keys = preprocess_context.ref_keys
        pred_keys = preprocess_context.pred_keys
        embedding_keys = (
            getattr(preprocess_context, "batch_correction_embedding_keys", None)
            or getattr(preprocess_context, "annotation_embedding_keys", None)
            or preprocess_context.embedding_keys
        )
        batch_keys = preprocess_context.batch_keys
        if not batch_keys:
            batch_keys = [
                col for col in adata.obs.columns if "batch" in col.lower()
            ]

        if verbose:
            print("Using precomputed context — skipping column detection")
            print(f"  Reference keys: {ref_keys}")
            print(f"  Predicted keys: {pred_keys}")
            print(f"  Embedding keys: {embedding_keys}")
            print(f"  Batch keys: {batch_keys}")

        for emb in embedding_keys:
            _safe = emb.replace("/", "_").replace(" ", "_")
            if (
                emb in preprocess_context.knn_paths
                or _safe in preprocess_context.knn_paths
            ):
                from ._cache import load_knn

                loaded = load_knn(
                    preprocess_context.annotation_dir, f"knn_{_safe}"
                )
                if loaded is not None:
                    emb_nn[emb] = NeighborResults(
                        indices=loaded[0], distances=loaded[1]
                    )

        for emb, payload in preprocess_context.neighbor_graphs.items():
            key_added = payload.get("key_added", f"neighbors_{emb}")
            if key_added not in adata.uns:
                adata.uns[key_added] = payload.get("uns_entry", {})
            conn_key = (
                payload.get("uns_entry", {})
                .get("connectivities_key", "connectivities")
            )
            dist_key = (
                payload.get("uns_entry", {})
                .get("distances_key", "distances")
            )
            if conn_key in payload and conn_key not in adata.obsp:
                adata.obsp[conn_key] = payload["connectivities"]
            if dist_key in payload and dist_key not in adata.obsp:
                adata.obsp[dist_key] = payload["distances"]
    else:
        detector = CheckAtlasColumnDetector(adata)
        params = detector.detect_all_parameters()
        ref_keys = [x[0] for x in params["annotation"]["reference"]]
        pred_keys = [x[0] for x in params["annotation"]["predicted"]]
        embedding_keys = [x[0] for x in params["clustering"]["embeddings"]]
        batch_keys = [x[0] for x in params.get("batch", [])]
        if not batch_keys:
            batch_keys = [
                col for col in adata.obs.columns if "batch" in col.lower()
            ]

    if not embedding_keys and not batch_keys and not ref_keys and not pred_keys:
        logger.warning(
            "No batch, reference, predicted, or embedding keys available. "
            "Skipping batch-correction metrics."
        )
        return pd.DataFrame(
            columns=[
                "Atlas Name",
                "Metric Name",
                "Embedding",
                "Batch Key",
                "Value",
                "Time (s)",
            ]
        )

    if metric_list is not None:
        metrics_list = [m for m in metric_list if m in METRICS_BATCH]
    elif all:
        metrics_list = METRICS_BATCH
    else:
        metrics_list = ["kbet", "iLISI", "graph_connectivity"]
        metrics_list = [m for m in metrics_list if m in METRICS_BATCH]

    if not metrics_list:
        return pd.DataFrame(
            columns=[
                "Atlas Name",
                "Metric Name",
                "Embedding",
                "Batch Key",
                "Value",
                "Time (s)",
            ]
        )

    results: list[dict] = []

    # LISI gets special handling: the ``lisi`` module dispatches to
    # iLISI (batch) and cLISI (cell-type) based on which labels are
    # passed in. We treat the two flavours as separate ``Metric Name``
    # values to preserve the long-format compatibility the scfm
    # diagnostic engine and MultiQC expect.
    batch_metrics = ["kbet", "pcr", "graph_connectivity"]

    pbar = tqdm(
        metrics_list,
        desc="Calculating Batch-Correction Metrics",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for metric in pbar:
        metric_start_time = time.time()
        pbar.set_description(f"Processing: {metric}")

        metric_module = getattr(batch_correction, metric)

        try:
            # ── LISI: iLISI on batch keys, cLISI on cell-type keys ─
            if metric == "lisi":
                if not emb_nn and _USE_JAX:
                    for emb in embedding_keys:
                        X_emb = np.asarray(adata.obsm[emb], dtype=np.float64)
                        emb_nn[emb] = _cal_knn(
                            X_emb, n_neighbors=90, backend="auto"
                        )

                # iLISI: per batch column
                if batch_keys:
                    for batch in batch_keys:
                        try:
                            if embedding_keys:
                                for emb in embedding_keys:
                                    X_emb = adata.obsm[emb]
                                    pair_start = time.time()
                                    if emb in emb_nn:
                                        val = metric_module.run_with_neighbors(
                                            emb_nn[emb],
                                            adata.obs[batch],
                                            verbose=False,
                                        )
                                    else:
                                        sig = inspect.signature(metric_module.run)
                                        kw = {}
                                        if "n_jobs" in sig.parameters:
                                            kw["n_jobs"] = n_jobs
                                        if "verbose" in sig.parameters:
                                            kw["verbose"] = False
                                        val = metric_module.run(
                                            X_emb, adata.obs[batch], **kw
                                        )
                                    pair_elapsed = time.time() - pair_start
                                    results.append(
                                        {
                                            "Atlas Name": atlas_name,
                                            "Metric Name": "iLISI",
                                            "Embedding": emb,
                                            "Batch Key": batch,
                                            "Value": val,
                                            "Time (s)": round(pair_elapsed, 3),
                                        }
                                    )
                            else:
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {}
                                if "n_jobs" in sig.parameters:
                                    kw["n_jobs"] = n_jobs
                                if "verbose" in sig.parameters:
                                    kw["verbose"] = False
                                val = metric_module.run(
                                    adata.X, adata.obs[batch], **kw
                                )
                                pair_elapsed = time.time() - pair_start
                                results.append(
                                    {
                                        "Atlas Name": atlas_name,
                                        "Metric Name": "iLISI",
                                        "Embedding": "X",
                                        "Batch Key": batch,
                                        "Value": val,
                                        "Time (s)": round(pair_elapsed, 3),
                                    }
                                )
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate iLISI for {batch}: {e}"
                            )

                # cLISI: per cell-type column (ref or pred)
                targets = list(set(ref_keys + pred_keys))
                for label in targets:
                    try:
                        if embedding_keys:
                            for emb in embedding_keys:
                                X_emb = adata.obsm[emb]
                                pair_start = time.time()
                                if emb in emb_nn:
                                    val = metric_module.run_with_neighbors(
                                        emb_nn[emb],
                                        adata.obs[label],
                                        verbose=False,
                                    )
                                else:
                                    sig = inspect.signature(metric_module.run)
                                    kw = {}
                                    if "n_jobs" in sig.parameters:
                                        kw["n_jobs"] = n_jobs
                                    if "verbose" in sig.parameters:
                                        kw["verbose"] = False
                                    val = metric_module.run(
                                        X_emb, adata.obs[label], **kw
                                    )
                                pair_elapsed = time.time() - pair_start
                                results.append(
                                    {
                                        "Atlas Name": atlas_name,
                                        "Metric Name": "cLISI",
                                        "Embedding": emb,
                                        "Batch Key": label,
                                        "Value": val,
                                        "Time (s)": round(pair_elapsed, 3),
                                    }
                                )
                        else:
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if "n_jobs" in sig.parameters:
                                kw["n_jobs"] = n_jobs
                            if "verbose" in sig.parameters:
                                kw["verbose"] = False
                            val = metric_module.run(adata.X, adata.obs[label], **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append(
                                {
                                    "Atlas Name": atlas_name,
                                    "Metric Name": "cLISI",
                                    "Embedding": "X",
                                    "Batch Key": label,
                                    "Value": val,
                                    "Time (s)": round(pair_elapsed, 3),
                                }
                            )
                    except Exception as e:
                        logger.warning(
                            f"Failed to calculate cLISI for {label}: {e}"
                        )

            # ── kBET / PCR: per (embedding, batch key) ─────────────
            elif metric in ("kbet", "pcr"):
                if not batch_keys:
                    continue
                for batch in batch_keys:
                    if embedding_keys:
                        for emb in embedding_keys:
                            try:
                                X_emb = adata.obsm[emb]
                                pair_start = time.time()
                                if (
                                    metric == "kbet"
                                    and hasattr(metric_module, "run_with_neighbors")
                                ):
                                    if emb in emb_nn:
                                        nn = emb_nn[emb]
                                        if nn.n_neighbors > 25:
                                            nn = nn.subset_neighbors(25)
                                        val = metric_module.run_with_neighbors(
                                            nn, adata.obs[batch], verbose=False
                                        )
                                    elif _USE_JAX:
                                        X_arr = np.asarray(X_emb, dtype=np.float64)
                                        nn = _cal_knn(
                                            X_arr, n_neighbors=25, backend="auto"
                                        )
                                        val = metric_module.run_with_neighbors(
                                            nn, adata.obs[batch], verbose=False
                                        )
                                    else:
                                        sig = inspect.signature(metric_module.run)
                                        kw = {}
                                        if "n_jobs" in sig.parameters:
                                            kw["n_jobs"] = n_jobs
                                        if "verbose" in sig.parameters:
                                            kw["verbose"] = False
                                        val = metric_module.run(
                                            X_emb, adata.obs[batch], **kw
                                        )
                                else:
                                    sig = inspect.signature(metric_module.run)
                                    kw = {}
                                    if "n_jobs" in sig.parameters:
                                        kw["n_jobs"] = n_jobs
                                    if "verbose" in sig.parameters:
                                        kw["verbose"] = False
                                    val = metric_module.run(
                                        X_emb, adata.obs[batch], **kw
                                    )
                                pair_elapsed = time.time() - pair_start
                                results.append(
                                    {
                                        "Atlas Name": atlas_name,
                                        "Metric Name": metric,
                                        "Embedding": emb,
                                        "Batch Key": batch,
                                        "Value": val,
                                        "Time (s)": round(pair_elapsed, 3),
                                    }
                                )
                            except Exception as e:
                                logger.warning(
                                    f"Failed to calculate {metric} for "
                                    f"{emb} vs {batch}: {e}"
                                )
                    else:
                        try:
                            pair_start = time.time()
                            if (
                                metric == "kbet"
                                and hasattr(metric_module, "run_with_neighbors")
                            ):
                                if "X" in emb_nn:
                                    nn = emb_nn["X"]
                                    if nn.n_neighbors > 25:
                                        nn = nn.subset_neighbors(25)
                                    val = metric_module.run_with_neighbors(
                                        nn, adata.obs[batch], verbose=False
                                    )
                                elif _USE_JAX:
                                    X_arr = np.asarray(adata.X, dtype=np.float64)
                                    if hasattr(adata.X, "toarray"):
                                        X_arr = adata.X.toarray().astype(np.float64)
                                    nn = _cal_knn(
                                        X_arr, n_neighbors=25, backend="auto"
                                    )
                                    val = metric_module.run_with_neighbors(
                                        nn, adata.obs[batch], verbose=False
                                    )
                                else:
                                    sig = inspect.signature(metric_module.run)
                                    kw = {}
                                    if "n_jobs" in sig.parameters:
                                        kw["n_jobs"] = n_jobs
                                    if "verbose" in sig.parameters:
                                        kw["verbose"] = False
                                    val = metric_module.run(
                                        adata.X, adata.obs[batch], **kw
                                    )
                            else:
                                sig = inspect.signature(metric_module.run)
                                kw = {}
                                if "n_jobs" in sig.parameters:
                                    kw["n_jobs"] = n_jobs
                                if "verbose" in sig.parameters:
                                    kw["verbose"] = False
                                val = metric_module.run(
                                    adata.X, adata.obs[batch], **kw
                                )
                            pair_elapsed = time.time() - pair_start
                            results.append(
                                {
                                    "Atlas Name": atlas_name,
                                    "Metric Name": metric,
                                    "Embedding": "X",
                                    "Batch Key": batch,
                                    "Value": val,
                                    "Time (s)": round(pair_elapsed, 3),
                                }
                            )
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate {metric} for {batch}: {e}"
                            )

            # ── graph_connectivity: per (embedding, cell-type) ─────
            elif metric == "graph_connectivity":
                targets = list(set(ref_keys + pred_keys))
                if not targets:
                    logger.info(
                        "graph_connectivity: no ref / pred cell-type column "
                        "detected; skipping."
                    )
                    continue
                if embedding_keys:
                    for emb in embedding_keys:
                        key_added = f"neighbors_{emb}"
                        try:
                            import scanpy as sc

                            sc.pp.neighbors(adata, use_rep=emb, key_added=key_added)
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate neighbors for {emb}: {e}"
                            )
                            continue

                        for label in targets:
                            try:
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {
                                    "neighbors_key": key_added,
                                    "label_key": label,
                                }
                                if "n_jobs" in sig.parameters:
                                    kw["n_jobs"] = n_jobs
                                if "verbose" in sig.parameters:
                                    kw["verbose"] = False
                                val = metric_module.run(adata, **kw)
                                pair_elapsed = time.time() - pair_start
                                results.append(
                                    {
                                        "Atlas Name": atlas_name,
                                        "Metric Name": "graph_connectivity",
                                        "Embedding": emb,
                                        "Batch Key": label,
                                        "Value": val,
                                        "Time (s)": round(pair_elapsed, 3),
                                    }
                                )
                            except Exception as e:
                                logger.warning(
                                    f"Failed to calculate graph_connectivity "
                                    f"for {emb} vs {label}: {e}"
                                )
                else:
                    for label in targets:
                        try:
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {"label_key": label}
                            if "n_jobs" in sig.parameters:
                                kw["n_jobs"] = n_jobs
                            if "verbose" in sig.parameters:
                                kw["verbose"] = False
                            val = metric_module.run(adata, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append(
                                {
                                    "Atlas Name": atlas_name,
                                    "Metric Name": "graph_connectivity",
                                    "Embedding": "Default",
                                    "Batch Key": label,
                                    "Value": val,
                                    "Time (s)": round(pair_elapsed, 3),
                                }
                            )
                        except Exception as e:
                            logger.warning(
                                f"Failed to calculate graph_connectivity "
                                f"for Default vs {label}: {e}"
                            )

            else:
                logger.warning(
                    f"Metric {metric} not categorized in batch_correction pipeline."
                )

        except Exception as e:
            logger.error(f"Error running batch-correction metric {metric}: {e}")

        metric_elapsed = time.time() - metric_start_time
        pbar.set_postfix_str(f"Time: {metric_elapsed:.2f}s", refresh=True)

    df = pd.DataFrame(results)

    if not df.empty and file_dir is not None and atlas_name is not None:
        os.makedirs(file_dir, exist_ok=True)
        wide_df = _pivot_batch_correction_to_wide(df, atlas_name)
        wide_path = os.path.join(file_dir, f"{atlas_name}.tsv")
        wide_df.to_csv(wide_path, sep="\t", index=False)
        logger.info(
            "MultiQC-compatible batch-correction table saved to %s", wide_path
        )

    try:
        from .batch_correction import kbet as _kbet

        _kbet._clear_knn_cache()
    except Exception:
        pass
    try:
        from .batch_correction import lisi as _lisi

        _lisi._clear_knn_cache()
    except Exception:
        pass

    return df


def cal_cluster(
    adata,
    atlas_name=None,
    metric_list=None,
    all_metrics=True,
    file_dir=None,
    n_jobs=-1,
    verbose=True,
    seed=42,
    preprocess_context=None,
):
    """
    Comprehensive clustering assessment pipeline.

    Calculates all clustering metrics to evaluate the quality of cluster assignments
    against embedding representations. Uses CheckAtlasColumnDetector to auto-detect
    embedding keys and cluster label columns.

    Args:
        adata (AnnData): Annotated data matrix.
        atlas_name (str): Name of the atlas for labeling results.
        metric_list (list, optional): List of metric names to calculate.
                           If provided, overrides ``all_metrics``.
        all_metrics (bool): If True, calculate all available cluster metrics.
                            If False, calculate a default subset.
                            Ignored if ``metric_list`` is provided.
        file_dir (str): Directory path where the results CSV will be saved.
                       If None, saves to current working directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress.
        seed (int): Random seed for reproducibility.
        preprocess_context (PreprocessContext, optional): If provided,
            column detection and distance-matrix precomputation are
            skipped and the context's precomputed data is reused.

    Returns:
        pd.DataFrame: Results dataframe with columns:
                      [Atlas Name, Metric Name, Embedding, Label Key, Value, Time (s)]
    """
    import gc
    import inspect

    from scipy.sparse import issparse

    from ..utils.col_detector import CheckAtlasColumnDetector

    # Set file directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)

    precomputed_dists = {}
    _ss_indices = None

    if preprocess_context is not None:
        embedding_keys = preprocess_context.cluster_embedding_keys or preprocess_context.embedding_keys
        label_keys = preprocess_context.cluster_label_keys

        # Exclude <4-dim embeddings (UMAP, t‑SNE) — useless for
        # kNN-based cluster metrics and double the workload.
        embedding_keys = [
            k for k in embedding_keys
            if k == "X" or (k in adata.obsm and adata.obsm[k].shape[1] > 3)
        ]

        if not label_keys:
            logger.warning(
                "No cluster labels in preprocess context. Skipping."
            )
            return pd.DataFrame()

        if verbose:
            print("Using precomputed context — skipping column detection")
            print(f"  Embeddings: {embedding_keys}")
            print(f"  Cluster labels: {label_keys}")

        # Load precomputed distance matrices for silhouette.
        # Large atlases use core-set subsampling; account for .tri
        # files built on the subset.
        _ss_indices = getattr(preprocess_context, "subsample_indices", None)
        _ss_n = len(_ss_indices) if _ss_indices is not None else 0
        _safe = lambda s: s.replace("/", "_").replace(" ", "_")
        for emb in embedding_keys:
            tri_path = os.path.join(
                preprocess_context.cluster_dir, f"dist_{_safe(emb)}.tri"
            )
            npy_path = tri_path.replace(".tri", ".npy")
            if os.path.exists(tri_path):
                from ._triangular import TriangularMatrix

                if _ss_n:
                    n_cells = _ss_n
                elif emb == "X":
                    n_cells = adata.X.shape[0]
                elif emb in adata.obsm:
                    n_cells = adata.obsm[emb].shape[0]
                else:
                    n_cells = adata.n_obs
                file_n = TriangularMatrix.n_from_file(tri_path)
                if file_n == n_cells:
                    precomputed_dists[emb] = TriangularMatrix(
                        n=n_cells, filepath=tri_path, mode="r"
                    )
                else:
                    logger.warning(
                        "Skipping stale distance matrix %s "
                        "(file N=%d does not match expected N=%d)",
                        tri_path,
                        file_n,
                        n_cells,
                    )
            elif os.path.exists(npy_path):
                precomputed_dists[emb] = np.load(npy_path)
    else:
        # Detect columns using CheckAtlasColumnDetector
        if verbose:
            print("Detecting embeddings and cluster labels...")
        detector = CheckAtlasColumnDetector(adata)
        params = detector.detect_all_parameters()

        embedding_keys = [x[0] for x in params["clustering"]["embeddings"]]
        label_keys = [x[0] for x in params["clustering"]["cluster_labels"]]
        ref_keys = [x[0] for x in params["annotation"]["reference"]]
        label_keys = label_keys + [r for r in ref_keys if r not in label_keys]

    if verbose:
        print(f"  Detected embeddings: {embedding_keys}")
        print(f"  Detected cluster labels: {label_keys}")

    if not embedding_keys:
        logger.warning(
            "No embeddings detected in adata.obsm. Trying adata.X as fallback."
        )
        embedding_keys = ["X"]

    if not label_keys:
        logger.warning(
            "No cluster labels detected in adata.obs. Cannot run cluster metrics."
        )
        return pd.DataFrame()

    # Define metrics to run
    if metric_list is not None:
        metrics_list = [m for m in metric_list if m in METRICS_CLUST]
    elif all_metrics:
        metrics_list = METRICS_CLUST
    else:
        # Default subset: fastest and most commonly used
        metrics_list = [
            "silhouette",
            "davies_bouldin",
            "calinski_harabasz",
        ]
        metrics_list = [m for m in metrics_list if m in METRICS_CLUST]

    results = []

    # Total combinations
    total_combos = len(embedding_keys) * len(label_keys) * len(metrics_list)
    if verbose:
        print(
            f"\nRunning {len(metrics_list)} metrics × "
            f"{len(embedding_keys)} embeddings × "
            f"{len(label_keys)} label keys = {total_combos} calculations\n"
        )

    # Create progress bar
    pbar = tqdm(
        total=total_combos,
        desc="Calculating Cluster Metrics",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
        disable=not verbose,
    )

    for emb_key in embedding_keys:
        # Extract embedding data
        if emb_key == "X":
            X_emb = adata.X
            if issparse(X_emb):
                X_emb = X_emb.toarray()
        else:
            if emb_key not in adata.obsm:
                logger.warning(
                    f"Embedding '{emb_key}' not found in adata.obsm. Skipping."
                )
                pbar.update(len(label_keys) * len(metrics_list))
                continue
            X_emb = np.asarray(adata.obsm[emb_key])

        if _ss_indices is not None:
            X_emb = X_emb[_ss_indices]

        for label_key in label_keys:
            if label_key not in adata.obs.columns:
                logger.warning(
                    f"Label key '{label_key}' not found in adata.obs. Skipping."
                )
                pbar.update(len(metrics_list))
                continue

            labels = np.asarray(adata.obs[label_key])
            if _ss_indices is not None:
                labels = labels[_ss_indices]

            # Check we have at least 2 clusters
            n_unique = len(np.unique(labels))
            if n_unique < 2:
                logger.warning(f"Label '{label_key}' has < 2 clusters. Skipping.")
                pbar.update(len(metrics_list))
                continue

            for metric_name in metrics_list:
                metric_start_time = time.time()
                pbar.set_description(f"{emb_key}|{label_key}: {metric_name}")

                try:
                    metric_module = getattr(cluster, metric_name)

                    # Build kwargs dynamically based on metric's signature
                    sig = inspect.signature(metric_module.run)
                    metric_params = sig.parameters

                    kwargs = {}

                    if "n_jobs" in metric_params:
                        kwargs["n_jobs"] = n_jobs
                    if "verbose" in metric_params:
                        kwargs["verbose"] = False
                    if "random_state" in metric_params:
                        kwargs["random_state"] = seed
                    if "seed" in metric_params:
                        kwargs["seed"] = seed
                    if "max_samples" in metric_params:
                        kwargs["max_samples"] = None  # disable subsampling

                    # Pass precomputed distances for silhouette when available
                    if "precomputed_dists" in metric_params and emb_key in precomputed_dists:
                        kwargs["precomputed_dists"] = precomputed_dists[emb_key]

                    # Call the metric
                    value = metric_module.run(X_emb, labels, **kwargs)

                    elapsed = time.time() - metric_start_time

                    results.append(
                        {
                            "Atlas Name": atlas_name,
                            "Metric Name": metric_name,
                            "Embedding": emb_key,
                            "Label Key": label_key,
                            "Value": value,
                            "Time (s)": round(elapsed, 3),
                        }
                    )

                    pbar.set_postfix_str(f"={value:.4f} ({elapsed:.1f}s)", refresh=True)

                except Exception as e:
                    elapsed = time.time() - metric_start_time
                    logger.warning(f"Failed {metric_name} for {emb_key}/{label_key}: {e}")
                    results.append(
                        {
                            "Atlas Name": atlas_name,
                            "Metric Name": metric_name,
                            "Embedding": emb_key,
                            "Label Key": label_key,
                            "Value": np.nan,
                            "Time (s)": round(elapsed, 3),
                        }
                    )

                pbar.update(1)

    pbar.close()

    # Create DataFrame
    df = pd.DataFrame(results)

    # Summary
    if verbose and not df.empty:
        total_time = df["Time (s)"].sum()
        print(f"\nTotal computation time: {total_time:.2f}s")
        print(
            f"Results: {len(df)} measurements across {len(df['Metric Name'].unique())} metrics"
        )

    gc.collect()
    return df


def calc_metric_cluster_scanpy(metric, adata, obs_key, obsm_key_representation):
    if metric in METRICS_CLUST:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        # Get the right module for the metric
        metric_module = getattr(cluster, metric)
        annotations = adata.obs[obs_key]
        if obsm_key_representation != "":
            count_repr = adata.obsm[obsm_key_representation]
            # execute the run function from metric_module
            metric_value = metric_module.run(count_repr, annotations)
            running_time = time.time() - start_time
            logger.debug(f"{metric} calc finished, duration {running_time}")
            return metric_value, running_time

        else:
            original_count = adata.X.toarray()
            # execute the run function from metric_module
            metric_value = metric_module.run(original_count, annotations)
            running_time = time.time() - start_time
            logger.debug(f"{metric} calc finished, duration {running_time}")
            return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_cluster_seurat(metric, seurat, obs_key, obsm_key_representation):
    if metric in METRICS_CLUST:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(cluster, metric)
        annotations = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        count_repr = ro.conversion.rpy2py(R_REDUCTION(seurat, obsm_key_representation))
        metric_value = metric_module.run(count_repr, annotations)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_annot_scanpy(metric, adata, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        # Get the right module for the metric
        metric_module = getattr(annot, metric)
        annotation = adata.obs[obs_key]
        ref_annotation = adata.obs[ref_obs]
        annotation, ref_annotation = annotation_to_num(annotation, ref_annotation)
        metric_value = metric_module.run(annotation, ref_annotation)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def calc_metric_annot_seurat(metric, seurat, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(annot, metric)
        annotation = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        ref_annotation = ro.conversion.rpy2py(R_ANNOT(seurat, ref_obs))
        # annotation, ref_annotation = annotation_to_num(
        #    annotation, ref_annotation
        # )
        metric_value = metric_module.run(annotation, ref_annotation)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def _safe_close_memmap(memmap_obj):
    """Force-close a numpy memmap so the backing .dat file can be deleted.

    ``del memmap_obj`` does not reliably release the file handle in time
    for a subsequent ``os.remove()``.  This helper walks the buffer chain
    to close the underlying mmap explicitly.
    """
    import mmap as _mmap_mod

    if memmap_obj is None:
        return
    # Walk the memmap object's inheritance chain to find the raw mmap
    for _attr in ("_mmap", "base"):
        inner = getattr(memmap_obj, _attr, None)
        if inner is not None and isinstance(inner, _mmap_mod.mmap):
            try:
                inner.close()
            except Exception:
                pass
    # Also try closing the memmap directly if it IS an mmap
    if isinstance(memmap_obj, _mmap_mod.mmap):
        try:
            memmap_obj.close()
        except Exception:
            pass


def _safe_close_triangular(tri):
    """Close the underlying mmap of a TriangularMatrix."""
    if tri is None:
        return
    tri._close_mmap()


def cal_dimred(
    adata,
    atlas_name=None,
    low_dim_keys=None,
    high_dim_key="X",
    metric_list=None,
    k_neighbors=30,
    n_samples=None,
    seed=42,
    n_jobs=-1,
    file_dir=None,
    verbose=True,
    use_cache=True,
):
    """Calculate dimensionality reduction metrics for multiple embeddings.

    For each ``low_dim_key`` the embedding is compared against the
    reference ``high_dim_key`` (defaults to ``adata.X`` — the raw gene
    expression matrix).

    Precomputation (distance matrices + kNN graphs) is performed once for
    the high‑dimensional reference and once per low‑dimensional embedding
    key.  Memory‑mapped files are used for N×N distance matrices on
    datasets with > 10 000 cells to keep RAM usage bounded.

    Parameters
    ----------
    adata : AnnData
    atlas_name : str, optional
        Used only for logging.
    low_dim_keys : list of str, optional
        ``.obsm`` keys to evaluate.  Defaults to all keys returned by
        :meth:`CheckAtlasColumnDetector.get_dimred_embeddings`.
    high_dim_key : str
        Reference key.  ``'X'`` means ``adata.X``.  (Default ``'X'``)
    metric_list : list of str, optional
        Metrics to compute.  Defaults to :data:`METRICS_DIMRED`.
    k_neighbors : int
    n_samples : int or None
        Number of cells to subsample.  ``None`` = use all cells.
    seed : int
    n_jobs : int
    file_dir : str, optional
        Directory for temporary memmap files.  Falls back to system
        ``/tmp``.
    verbose : bool

    Returns
    -------
    pd.DataFrame
        Long‑format table with columns
        ``[Metric Name, Low Dim Key, High Dim Key, Value, Time (s)]``.
        Caller is responsible for pivoting and saving.
    """
    import gc
    import inspect
    import tempfile
    import uuid

    from sklearn.metrics import pairwise_distances
    from sklearn.neighbors import NearestNeighbors

    from ..utils.col_detector import CheckAtlasColumnDetector

    if metric_list is None:
        metric_list = METRICS_DIMRED
    metric_list = [m for m in metric_list if m in METRICS_DIMRED]
    if not metric_list:
        return pd.DataFrame()

    if low_dim_keys is None:
        detector = CheckAtlasColumnDetector(adata)
        low_dim_keys = detector.get_dimred_embeddings()

    low_dim_keys = [k for k in low_dim_keys if k != high_dim_key]
    if not low_dim_keys:
        logger.warning(
            "No low-dim embeddings to evaluate (only high-dim key '%s' present). "
            "Skipping dimred metrics.",
            high_dim_key,
        )
        return pd.DataFrame()

    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        np.random.seed(seed)
        sample_indices = np.random.choice(n_obs, n_samples, replace=False)
        n_cells = n_samples
    else:
        sample_indices = np.arange(n_obs)
        n_cells = n_obs

    if high_dim_key == "X":
        high_dim_data = adata.X[sample_indices]
        if hasattr(high_dim_data, "toarray"):
            high_dim_data = high_dim_data.toarray()
    else:
        if high_dim_key not in adata.obsm_keys():
            logger.warning(
                "High-dim key '%s' not found in adata.obsm. "
                "Skipping dimred metrics.",
                high_dim_key,
            )
            return pd.DataFrame()
        high_dim_data = adata.obsm[high_dim_key][sample_indices]

    high_n_features = high_dim_data.shape[1]

    use_memmap = n_cells > 10000
    all_memmap_files = []

    if file_dir:
        temp_dir = file_dir
    else:
        # Try default, fall back to /data if system temp is on a small partition
        default_tmp = os.path.join(os.path.expanduser("~"), ".checkatlas", "tmp")
        temp_dir = default_tmp
        try:
            import shutil

            _usage = shutil.disk_usage(default_tmp)
            if _usage.free < 50 * 1024**3:  # < 50 GB free
                _fallback = "/data/tmp" if os.path.exists("/data/tmp") else "/tmp"
                # Also try nextflow TMPDIR if set
                _nf_tmp = os.environ.get("TMPDIR", "")
                if _nf_tmp and os.path.exists(_nf_tmp):
                    _fallback = _nf_tmp
                temp_dir = os.path.join(_fallback, ".checkatlas_tmp")
                if verbose:
                    print(f"  Low disk space on default temp; using {temp_dir}")
        except Exception:
            pass
    os.makedirs(temp_dir, exist_ok=True)

    # ── Persistent cache check ───────────────────────────────────
    _from_cache = False
    _cache_low_dim = {}
    if use_cache and atlas_name:
        _fp = compute_fingerprint(
            n_cells=n_cells,
            n_features=high_n_features,
            embedding_keys=low_dim_keys,
            embedding_shapes={
                k: tuple(adata.obsm[k][sample_indices].shape) for k in low_dim_keys
            },
            k_neighbors=k_neighbors,
            source_path=getattr(adata, "filename", None),
        )
        _cached = load_dimred_cache(temp_dir, _fp, n_cells, low_dim_keys)
        if _cached is not None:
            high_dim_dists = _cached["high_dim_dists"]
            high_knn_dists = _cached["high_knn_dists"]
            high_knn_indices = _cached["high_knn_indices"]
            _cache_low_dim = _cached["low_dim"]
            _from_cache = True
            # Large atlases: the cached data was built on a core-set
            # subsample.  Re-extract high_dim_data and adjust n_cells
            # so all downstream computation uses the same subset.
            _cached_subsample = _cached.get("subsample_indices")
            if _cached_subsample is not None:
                sample_indices = _cached_subsample
                n_cells = len(sample_indices)
                if high_dim_key == "X":
                    high_dim_data = adata.X[sample_indices]
                    if hasattr(high_dim_data, "toarray"):
                        high_dim_data = high_dim_data.toarray()
                else:
                    high_dim_data = adata.obsm[high_dim_key][sample_indices]
            if verbose:
                print("  [CACHE HIT] Reusing precomputed distances & kNN")

    run_id = str(uuid.uuid4())[:8]
    chunk_size = min(1000, n_cells)

    if verbose:
        print(
            f"  Reference: {high_dim_key}  "
            f"({n_cells:,} cells  ×  {high_n_features} features)"
        )
        print(f"  Embeddings: {low_dim_keys}")
        print(f"  Metrics: {metric_list}")
        backend = "GPU (JAX)" if (_JAX_AVAILABLE and _GPU_AVAILABLE) else "CPU"
        print(f"  Backend: {backend}")

    # ── GPU path: single-kernel distance matrix + kNN (Phase 3) ────
    # Memory: N² float32 ≈ 4·N² bytes → with intermediate buffers ≈ 4·N²·3.5 bytes
    # A100 40 GB: safe N ≤ 50 000 (≈ 32.6 GB total)
    # For 50k < N ≤ 150k: chunked GPU + memmap on /data (avoid disk IO bottleneck)
    _GPU_SINGLE_SHOT = _JAX_AVAILABLE and _GPU_AVAILABLE and (n_cells <= 50000)
    _GPU_CHUNKED = _JAX_AVAILABLE and _GPU_AVAILABLE and (50000 < n_cells <= 150000)

    _DIST_METRICS = frozenset(
        (
            "kruskal_stress",
            "spearman_rho",
            "dCor",
            "trustworthiness",  # OOM-safe via the per-row memmap path
            "continuity",       # OOM-safe via the per-row memmap path
        )
    )

    _low_dim_precomputed = {}  # collect for saving to cache later

    if _from_cache:
        # ── Recompute high-dim distances on cache hit if missing ─
        # Old cache files may have been written by the CPU path,
        # which never persisted high-dists (only TriangularMatrix
        # memmaps are saved by save_dimred_cache).  When the metric
        # list includes distance-based metrics (kruskal_stress,
        # spearman_rho, dCor, trust, continuity) we need the full
        # high-dim distance matrix.  Recompute it on the GPU and
        # persist to the cache so future runs skip this entirely.
        if (
            high_dim_dists is None
            and bool(set(metric_list) & _DIST_METRICS)
        ):
            _high_tri_path = (
                os.path.join(temp_dir, "high_dists.tri")
                if use_cache and atlas_name
                else None
            )

            if _GPU_SINGLE_SHOT:
                if verbose:
                    print(
                        "  [CACHE] Reconstructing high-dim distances on GPU "
                        "(single-shot, missing from prior cache)..."
                    )
                _dists_np = pdist_squareform(high_dim_data)
                if _high_tri_path is not None:
                    _high_tri = TriangularMatrix(
                        n=n_cells, filepath=_high_tri_path, mode="w+"
                    )
                    _chunk = min(5000, n_cells)
                    for _i in range(0, n_cells, _chunk):
                        _end = min(_i + _chunk, n_cells)
                        store_upper_triangle(
                            _high_tri._data, _dists_np[_i:_end, :],
                            _i, 0, n_cells,
                        )
                    _high_tri.flush()
                    high_dim_dists = _high_tri
                    del _dists_np
                else:
                    high_dim_dists = _dists_np
                gc.collect()
            elif _GPU_CHUNKED:
                if verbose:
                    print(
                        "  [CACHE] Reconstructing high-dim distances on GPU "
                        "(chunked, missing from prior cache)..."
                    )
                if _high_tri_path is None:
                    _high_tri_path = os.path.join(
                        temp_dir, f"high_dists_{run_id}.tri"
                    )
                high_dim_dists = TriangularMatrix(
                    n=n_cells, filepath=_high_tri_path, mode="w+"
                )
                import numpy as _np
                _knn_res = compute_neighbors(
                    _np.asarray(high_dim_data, dtype=_np.float64),
                    n_neighbors=k_neighbors + 1,
                    backend="auto",
                    tri_memmap=high_dim_dists,
                )
                high_knn_dists = _knn_res.distances
                high_knn_indices = _knn_res.indices
                high_dim_dists.flush()
                gc.collect()
            else:
                # No GPU: CPU fallback.  Slow but functional.
                if verbose:
                    print(
                        "  [CACHE] Reconstructing high-dim distances on CPU "
                        "(no GPU available)..."
                    )
                from sklearn.metrics import pairwise_distances
                high_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
                for _i in range(0, n_cells, chunk_size):
                    _end = min(_i + chunk_size, n_cells)
                    high_dim_dists[_i:_end, :] = pairwise_distances(
                        high_dim_data[_i:_end], high_dim_data, n_jobs=n_jobs
                    )
                np.fill_diagonal(high_dim_dists, 0)

    elif _GPU_SINGLE_SHOT:
        import jax.numpy as jnp

        high_dim_dists = pdist_squareform(high_dim_data)  # GPU matmul → (n,n) float32
        high_dists_jax = jnp.asarray(high_dim_dists, dtype=jnp.float32)

        # kNN from precomputed distance matrix on GPU
        import jax

        high_knn_dists_jax, high_knn_indices_jax = jax.lax.approx_min_k(
            high_dists_jax, k=k_neighbors + 1
        )
        high_knn_dists = _get_ndarray(high_knn_dists_jax)
        high_knn_indices = _get_ndarray(high_knn_indices_jax)

        # ── Persist GPU distance matrix for cache reuse ───────
        # save_dimred_cache only stores TriangularMatrix memmaps;
        # convert the GPU numpy array so subsequent runs don't
        # have to recompute this expensive (N² × D) matrix.
        if use_cache and atlas_name:
            _high_tri_path = os.path.join(temp_dir, "high_dists.tri")
            _high_tri = TriangularMatrix(
                n=n_cells, filepath=_high_tri_path, mode="w+"
            )
            _chunk = min(5000, n_cells)
            for _i in range(0, n_cells, _chunk):
                _end = min(_i + _chunk, n_cells)
                store_upper_triangle(
                    _high_tri._data, high_dim_dists[_i:_end, :],
                    _i, 0, n_cells,
                )
            _high_tri.flush()
            high_dim_dists = _high_tri

    elif _GPU_CHUNKED:
        # ── Chunked GPU path: fused kNN + distance matrix ──
        # kNN: streaming GPU, auto-detect for large atlases
        # Distances: upper‑triangle float16 memmap written in same pass
        import jax.numpy as jnp

        _qchunk = 15000  # query rows per chunk
        _rchunk = 10000  # reference columns per chunk

        if verbose:
            print(
                f"  Chunked GPU/streaming: {n_cells:,} cells, "
                f"query={_qchunk}, ref={_rchunk}"
            )
            print(f"  Storing in: {temp_dir}")

        # Create TriangularMatrix BEFORE kNN pass to fuse the operations
        high_dists_path = os.path.join(temp_dir, f"high_dists_{run_id}.tri")
        all_memmap_files.append(high_dists_path)
        high_dim_dists = TriangularMatrix(n=n_cells, filepath=high_dists_path, mode="w+")

        if verbose:
            print("  Computing kNN + distances (fused GPU)...")
        high_knn_results = compute_neighbors(
            np.asarray(high_dim_data, dtype=np.float64),
            n_neighbors=k_neighbors + 1,
            backend="auto",
            tri_memmap=high_dim_dists,
        )
        high_knn_dists = high_knn_results.distances
        high_knn_indices = high_knn_results.indices
        high_dim_dists.flush()
        gc.collect()
    else:
        # ─── CPU path (unchanged) ───────────────────────────────────
        # ── High‑dim kNN ──
        nbrs_high = NearestNeighbors(n_neighbors=k_neighbors + 1, n_jobs=n_jobs).fit(
            high_dim_data
        )
        high_knn_dists, high_knn_indices = nbrs_high.kneighbors(high_dim_data)

        # ── High‑dim distance matrix (only when distance‑based metrics
        #     are in the list) ──
        need_high_dists = bool(set(metric_list) & _DIST_METRICS)
        high_dim_dists = None

        if need_high_dists:
            if verbose:
                print(
                    "  Precomputing high-dim distance matrix "
                    f"({n_cells:,}×{n_cells:,})…"
                )
            if use_memmap:
                high_dists_path = os.path.join(temp_dir, f"high_dists_{run_id}.tri")
                all_memmap_files.append(high_dists_path)
                high_dim_dists = TriangularMatrix(
                    n=n_cells, filepath=high_dists_path, mode="w+"
                )
                for i in tqdm(
                    range(0, n_cells, chunk_size),
                    desc="High-Dim Distances",
                    disable=not verbose,
                ):
                    end = min(i + chunk_size, n_cells)
                    block = pairwise_distances(
                        high_dim_data[i:end], high_dim_data, n_jobs=n_jobs
                    )
                    store_upper_triangle(high_dim_dists._data, block, i, 0, n_cells)
                    high_dim_dists.flush()
            else:
                high_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
                for i in tqdm(
                    range(0, n_cells, chunk_size),
                    desc="High-Dim Distances",
                    disable=not verbose,
                ):
                    end = min(i + chunk_size, n_cells)
                    high_dim_dists[i:end, :] = pairwise_distances(
                        high_dim_data[i:end], high_dim_data, n_jobs=n_jobs
                    )
            gc.collect()

    # ── end of precomputation block ──

    total_combos = len(low_dim_keys) * len(metric_list)
    if verbose:
        print(
            f"\n  {len(metric_list)} metrics × {len(low_dim_keys)} embeddings"
            f" = {total_combos} calculations\n"
        )

    results = []
    pbar = tqdm(
        total=total_combos,
        desc="Calculating Dimred Metrics",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} " "[{elapsed}<{remaining}]",
        disable=not verbose,
    )

    for low_dim_key in low_dim_keys:
        low_dim_data = adata.obsm[low_dim_key][sample_indices]

        low_dim_dists = None
        low_dists_path = None

        if _from_cache and low_dim_key in _cache_low_dim:
            # ── Load from cache ───────────────────────────────
            _c = _cache_low_dim[low_dim_key]
            low_dim_dists = _c.get("dists")
            low_knn_dists = _c["knn_dists"]
            low_knn_indices = _c["knn_indices"]

            # ── Recompute missing distance matrix ─────────────
            # Cache may have kNN but not distances (GPU numpy
            # arrays are not persisted — save_dimred_cache only
            # stores TriangularMatrix memmaps).  Compute the
            # distance matrix once so every metric can share it.
            if low_dim_dists is None and low_knn_dists is not None:
                if _GPU_SINGLE_SHOT:
                    import jax.numpy as _jnp

                    _dists_np = pdist_squareform(low_dim_data)
                    if use_cache and atlas_name:
                        _safe = low_dim_key.replace("/", "_").replace(" ", "_")
                        _tri_path = os.path.join(
                            temp_dir, f"low_dists_{_safe}.tri"
                        )
                        _tri = TriangularMatrix(
                            n=n_cells, filepath=_tri_path, mode="w+"
                        )
                        _chunk = min(5000, n_cells)
                        for _i in range(0, n_cells, _chunk):
                            _end = min(_i + _chunk, n_cells)
                            store_upper_triangle(
                                _tri._data, _dists_np[_i:_end, :],
                                _i, 0, n_cells,
                            )
                        _tri.flush()
                        low_dim_dists = _tri
                        del _dists_np
                    else:
                        low_dim_dists = _dists_np
                elif n_cells > 10000:
                    _safe = low_dim_key.replace("/", "_").replace(" ", "_")
                    _tri_path = os.path.join(
                        temp_dir, f"low_dists_{_safe}_{run_id}.tri"
                    )
                    _tri = TriangularMatrix(
                        n=n_cells, filepath=_tri_path, mode="w+"
                    )
                    for _i in range(0, n_cells, chunk_size):
                        _end = min(_i + chunk_size, n_cells)
                        _block = pairwise_distances(
                            low_dim_data[_i:_end], low_dim_data,
                            n_jobs=n_jobs,
                        )
                        store_upper_triangle(
                            _tri._data, _block, _i, 0, n_cells
                        )
                    _tri.flush()
                    low_dim_dists = _tri
                else:
                    need_d = bool(set(metric_list) & _DIST_METRICS)
                    if need_d:
                        low_dim_dists = np.zeros(
                            (n_cells, n_cells), dtype=np.float32
                        )
                        for _i in range(0, n_cells, chunk_size):
                            _end = min(_i + chunk_size, n_cells)
                            low_dim_dists[_i:_end, :] = pairwise_distances(
                                low_dim_data[_i:_end], low_dim_data,
                                n_jobs=n_jobs,
                            )
                gc.collect()

        elif _GPU_SINGLE_SHOT:
            # ── GPU path: single-kernel distance + kNN ──────────
            import jax
            import jax.numpy as jnp

            low_dim_dists = pdist_squareform(low_dim_data)
            low_dists_jax = jnp.asarray(low_dim_dists, dtype=jnp.float32)
            low_knn_dists_jax, low_knn_indices_jax = jax.lax.approx_min_k(
                low_dists_jax, k=k_neighbors + 1
            )
            low_knn_dists = _get_ndarray(low_knn_dists_jax)
            low_knn_indices = _get_ndarray(low_knn_indices_jax)

            # ── Persist GPU distance matrix for cache reuse ──
            if use_cache and atlas_name:
                _safe = low_dim_key.replace("/", "_").replace(" ", "_")
                _tri_path = os.path.join(
                    temp_dir, f"low_dists_{_safe}.tri"
                )
                _tri = TriangularMatrix(
                    n=n_cells, filepath=_tri_path, mode="w+"
                )
                _chunk = min(5000, n_cells)
                for _i in range(0, n_cells, _chunk):
                    _end = min(_i + _chunk, n_cells)
                    store_upper_triangle(
                        _tri._data, low_dim_dists[_i:_end, :],
                        _i, 0, n_cells,
                    )
                _tri.flush()
                low_dim_dists = _tri

        elif _GPU_CHUNKED:
            # ── GPU for low-dim when features ≤ 200 dims ─────────
            # Low-dim embeddings fit on GPU easily.
            # High-dim raw data needs CPU memmap for N×N storage.
            low_ndim = low_dim_data.shape[1]
            if low_ndim <= 200 and n_cells <= 50000:
                # One-shot GPU: small feature dim → N² fits
                import jax
                import jax.numpy as jnp

                low_dim_dists = pdist_squareform(low_dim_data)
                low_dists_jax = jnp.asarray(low_dim_dists, dtype=jnp.float32)
                low_knn_dists_jax, low_knn_indices_jax = jax.lax.approx_min_k(
                    low_dists_jax, k=k_neighbors + 1
                )
                low_knn_dists = _get_ndarray(low_knn_dists_jax)
                low_knn_indices = _get_ndarray(low_knn_indices_jax)

                # ── Persist GPU distance matrix for cache reuse ──
                if use_cache and atlas_name:
                    _safe = low_dim_key.replace("/", "_").replace(" ", "_")
                    _tri_path = os.path.join(
                        temp_dir, f"low_dists_{_safe}.tri"
                    )
                    _tri = TriangularMatrix(
                        n=n_cells, filepath=_tri_path, mode="w+"
                    )
                    _chunk = min(5000, n_cells)
                    for _i in range(0, n_cells, _chunk):
                        _end = min(_i + _chunk, n_cells)
                        store_upper_triangle(
                            _tri._data, low_dim_dists[_i:_end, :],
                            _i, 0, n_cells,
                        )
                    _tri.flush()
                    low_dim_dists = _tri

            elif low_ndim <= 200:
                # Fused GPU streaming kNN + optional distance matrix
                need_low_dists = bool(set(metric_list) & _DIST_METRICS)
                if need_low_dists:
                    low_dists_path = os.path.join(
                        temp_dir,
                        f"low_dists_{low_dim_key.replace('/', '_')}_{run_id}.tri",
                    )
                    all_memmap_files.append(low_dists_path)
                    low_dim_dists = TriangularMatrix(
                        n=n_cells, filepath=low_dists_path, mode="w+"
                    )
                else:
                    low_dim_dists = None

                kNN = compute_neighbors(
                    np.asarray(low_dim_data, dtype=np.float64),
                    n_neighbors=k_neighbors + 1,
                    backend="auto",
                    tri_memmap=low_dim_dists,
                )
                low_knn_dists = kNN.distances
                low_knn_indices = kNN.indices
                if need_low_dists and low_dim_dists is not None:
                    low_dim_dists.flush()
                gc.collect()
            else:
                # High-dim X — CPU memmap (sklearn) → float16 upper‑triangle
                nbrs_low = NearestNeighbors(
                    n_neighbors=k_neighbors + 1, n_jobs=n_jobs
                ).fit(low_dim_data)
                low_knn_dists, low_knn_indices = nbrs_low.kneighbors(low_dim_data)
                low_dists_path = os.path.join(
                    temp_dir,
                    f"low_dists_{low_dim_key.replace('/', '_')}_{run_id}.tri",
                )
                all_memmap_files.append(low_dists_path)
                low_dim_dists = TriangularMatrix(
                    n=n_cells, filepath=low_dists_path, mode="w+"
                )
                for i in tqdm(
                    range(0, n_cells, chunk_size),
                    desc=f"Low-Dim Distances ({low_dim_key})",
                    disable=not verbose,
                ):
                    end = min(i + chunk_size, n_cells)
                    block = pairwise_distances(
                        low_dim_data[i:end], low_dim_data, n_jobs=n_jobs
                    )
                    store_upper_triangle(low_dim_dists._data, block, i, 0, n_cells)
                    low_dim_dists.flush()
                gc.collect()
        else:
            # ── CPU path ────────────────────────────────────────
            nbrs_low = NearestNeighbors(n_neighbors=k_neighbors + 1, n_jobs=n_jobs).fit(
                low_dim_data
            )
            low_knn_dists, low_knn_indices = nbrs_low.kneighbors(low_dim_data)

            need_low_dists = bool(set(metric_list) & _DIST_METRICS)
            if need_low_dists:
                if use_memmap:
                    low_dists_path = os.path.join(
                        temp_dir,
                        f"low_dists_{low_dim_key.replace('/', '_')}_{run_id}.tri",
                    )
                    all_memmap_files.append(low_dists_path)
                    low_dim_dists = TriangularMatrix(
                        n=n_cells, filepath=low_dists_path, mode="w+"
                    )
                    for i in tqdm(
                        range(0, n_cells, chunk_size),
                        desc=f"Low-Dim Distances ({low_dim_key})",
                        disable=not verbose,
                    ):
                        end = min(i + chunk_size, n_cells)
                        block = pairwise_distances(
                            low_dim_data[i:end], low_dim_data, n_jobs=n_jobs
                        )
                        store_upper_triangle(low_dim_dists._data, block, i, 0, n_cells)
                        low_dim_dists.flush()
                else:
                    low_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
                    for i in tqdm(
                        range(0, n_cells, chunk_size),
                        desc=f"Low-Dim Distances ({low_dim_key})",
                        disable=not verbose,
                    ):
                        end = min(i + chunk_size, n_cells)
                        low_dim_dists[i:end, :] = pairwise_distances(
                            low_dim_data[i:end], low_dim_data, n_jobs=n_jobs
                        )
                gc.collect()

        # ── Capture precomputed low-dim data for cache ────────
        if use_cache and atlas_name:
            _low_dim_precomputed[low_dim_key] = {
                "dists": low_dim_dists,
                "knn_indices": low_knn_indices,
                "knn_dists": low_knn_dists,
            }

        for metric_name in metric_list:
            pbar.set_description(f"{low_dim_key}: {metric_name}")
            t_start = time.time()

            # ── Inspect metric signature to decide on to_dense ──────
            metric_module = getattr(dimred, metric_name)
            sig = inspect.signature(metric_module.run)
            params = sig.parameters

            # ── TriangularMatrix is passed through as-is ────────────
            # The _rank_penalty helper handles TriangularMatrix inputs
            # natively via the per-row CPU read / GPU chunked paths
            # (see checkatlas/metrics/dimred/_rank_penalty.py).  For
            # large memmap-backed matrices, that per-row path is the
            # fastest documented backend — faster than materialising a
            # full dense N×N float32 buffer via to_dense() (which would
            # be ≈ 10 GB for N=50_000 and a hard OOM for N=100_000+).
            # We therefore do NOT call to_dense() here, regardless of N.

            try:
                kwargs = {}
                if "low_dim_key" in params:
                    kwargs["low_dim_key"] = low_dim_key
                if "high_dim_key" in params:
                    kwargs["high_dim_key"] = high_dim_key
                if "k_neighbors" in params:
                    kwargs["k_neighbors"] = k_neighbors
                if "n_samples" in params:
                    kwargs["n_samples"] = n_samples
                if "n_jobs" in params:
                    kwargs["n_jobs"] = n_jobs
                if "seed" in params:
                    kwargs["seed"] = seed
                if "verbose" in params:
                    kwargs["verbose"] = False
                if "rank_backend" in params:
                    # Trust/continuity support "auto" | "jax_single_shot"
                    # | "jax_chunked" | "cpu".  "auto" picks the
                    # fastest backend based on N + GPU availability.
                    kwargs["rank_backend"] = "auto"

                if "precomputed_high_knn" in params:
                    kwargs["precomputed_high_knn"] = high_knn_indices
                if "precomputed_low_knn" in params:
                    kwargs["precomputed_low_knn"] = low_knn_indices
                if "precomputed_high_knn_dists" in params:
                    kwargs["precomputed_high_knn_dists"] = high_knn_dists
                if "precomputed_low_knn_dists" in params:
                    kwargs["precomputed_low_knn_dists"] = low_knn_dists
                if "precomputed_high_dists" in params:
                    kwargs["precomputed_high_dists"] = high_dim_dists
                if "precomputed_low_dists" in params:
                    kwargs["precomputed_low_dists"] = low_dim_dists

                value = metric_module.run(adata, **kwargs)
                elapsed = round(time.time() - t_start, 3)

                results.append(
                    {
                        "Metric Name": metric_name,
                        "Low Dim Key": low_dim_key,
                        "High Dim Key": high_dim_key,
                        "Value": value,
                        "Time (s)": elapsed,
                    }
                )
            except Exception as e:
                elapsed = round(time.time() - t_start, 3)
                logger.warning("Failed %s on %s: %s", metric_name, low_dim_key, e)
                results.append(
                    {
                        "Metric Name": metric_name,
                        "Low Dim Key": low_dim_key,
                        "High Dim Key": high_dim_key,
                        "Value": np.nan,
                        "Time (s)": elapsed,
                    }
                )

            pbar.update(1)

        if low_dim_dists is not None:
            # When caching, defer the close until AFTER save_dimred_cache
            # has moved/renamed the file.  Closing here invalidates the
            # mmap that save_triangular still needs to operate on.
            if not use_cache:
                try:
                    if isinstance(low_dim_dists, TriangularMatrix):
                        _safe_close_triangular(low_dim_dists)
                    else:
                        _safe_close_memmap(low_dim_dists)
                except Exception:
                    pass
            if not use_cache:
                if low_dists_path and os.path.exists(low_dists_path):
                    try:
                        os.remove(low_dists_path)
                    except Exception:
                        pass
                    if low_dists_path in all_memmap_files:
                        all_memmap_files.remove(low_dists_path)
        gc.collect()

    pbar.close()

    # ── Save to persistent cache ──
    # On a fresh run we save the full precomputation.  On a cache hit
    # we still save IF we had to recompute anything (e.g. the high-dim
    # distances were missing from the prior cache); in that case the
    # fingerprint is unchanged but the augmented cache will be valid
    # for future runs.
    _save_cache = use_cache and atlas_name and high_dim_dists is not None
    if _save_cache:
        _low_info = {}
        for emb in low_dim_keys:
            _low_info[emb] = {
                "dists": _low_dim_precomputed.get(emb, {}).get("dists"),
                "knn_indices": _low_dim_precomputed.get(emb, {}).get("knn_indices"),
                "knn_dists": _low_dim_precomputed.get(emb, {}).get("knn_dists"),
            }
        try:
            save_dimred_cache(
                temp_dir,
                _fp,
                high_dim_dists=(
                    high_dim_dists
                    if isinstance(high_dim_dists, TriangularMatrix)
                    else None
                ),
                high_knn_dists=high_knn_dists,
                high_knn_indices=high_knn_indices,
                low_dim_data=_low_info,
                low_dim_keys=low_dim_keys,
            )
        except Exception as exc:
            logger.warning("Failed to save dimred cache: %s", exc)

    # Now safe to close any memmaps — save_dimred_cache has either
    # moved them to their canonical location or copied the data.
    if high_dim_dists is not None:
        try:
            if isinstance(high_dim_dists, TriangularMatrix):
                _safe_close_triangular(high_dim_dists)
            else:
                _safe_close_memmap(high_dim_dists)
        except Exception:
            pass
        gc.collect()

    # Close any low-dim memmaps (deferred from the per-embedding loop
    # so save_dimred_cache could read them above).
    for emb in low_dim_keys:
        _ld = _low_dim_precomputed.get(emb, {}).get("dists")
        if _ld is not None and isinstance(_ld, TriangularMatrix):
            try:
                _safe_close_triangular(_ld)
            except Exception:
                pass

    # Delete leftover temp files (only when NOT caching)
    if not use_cache:
        for fpath in all_memmap_files:
            try:
                if os.path.exists(fpath):
                    os.remove(fpath)
            except OSError as exc:
                logger.warning("Failed to clean up temp file %s: %s", fpath, exc)

    if "_dimred_cache" in adata.uns:
        del adata.uns["_dimred_cache"]

    gc.collect()

    df = pd.DataFrame(results)

    if verbose and not df.empty:
        total_time = df["Time (s)"].sum()
        n_results = len(df)
        print(f"\nTotal computation time: {total_time:.2f}s")
        print(f"Results: {n_results} measurements across " f"{len(metric_list)} metrics")

    return df


def calc_metric_dimred(metric, adata, obsm_key):
    """
    Calculate a single dimensionality reduction metric (legacy function).

    For comprehensive assessment, use cal_dimred() instead.
    """
    if metric in METRICS_DIMRED:
        start_time = time.time()
        # Get the right module for the metric
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(dimred, metric)
        print(metric_module)
        high_dim_counts = adata.X
        low_dim_counts = adata.obsm[obsm_key]
        # execute the run function from metric_module
        metric_value = metric_module.run(high_dim_counts, low_dim_counts)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"dimensionality reduction metric."
            f"\nList of dim. red. metrics: {METRICS_DIMRED}"
        )
        return -1


def run_all_metrics(
    adata=None,
    atlas_path=None,
    atlas_name=None,
    file_dir=None,
    n_jobs=-1,
    verbose=True,
    seed=42,
    # Annotation params
    run_annotation=True,
    all_annot_metrics=True,
    # Clustering params
    run_clustering=True,
    all_cluster_metrics=True,
    # Dimred params
    run_dimred=True,
    all_dimred_metrics=True,
    low_dim_key="X_umap",
    high_dim_key="X",
    k_neighbors=30,
    n_samples=None,
    use_memmap=True,
    temp_dir=None,
    chunk_size=1000,
):
    """
    Unified metrics pipeline for CheckAtlas.

    The standalone entry point that runs ALL metric tasks (annotation, clustering,
    dimensionality reduction) on a single-cell atlas and produces a consolidated
    results CSV. Designed for biologists and bioinformaticians — just provide
    an atlas path and get comprehensive quality metrics.

    Usage:
        >>> from checkatlas.metrics.metrics import run_all_metrics
        >>> results = run_all_metrics(atlas_path="my_atlas.h5ad")

        >>> # Or with a preloaded AnnData:
        >>> results = run_all_metrics(adata=my_adata, atlas_name="my_atlas")

    Args:
        adata (AnnData, optional): Preloaded AnnData object. If None, loads from atlas_path.
        atlas_path (str, optional): Path to .h5ad atlas file. Used if adata is None.
        atlas_name (str, optional): Name for the atlas (used in output). If None,
                                   derived from atlas_path filename.
        file_dir (str, optional): Directory for output CSV. Defaults to current directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress and summary.
        seed (int): Random seed for reproducibility.

        run_annotation (bool): Whether to run annotation metrics.
        all_annot_metrics (bool): If True, run all annotation metrics.

        run_clustering (bool): Whether to run clustering metrics.
        all_cluster_metrics (bool): If True, run all clustering metrics.

        run_dimred (bool): Whether to run dimensionality reduction metrics.
        all_dimred_metrics (bool): If True, run all dimred metrics.
        low_dim_key (str): Key in adata.obsm for low-dimensional embedding.
        high_dim_key (str): Key for high-dimensional data ('X' = adata.X).
        k_neighbors (int): Number of neighbors for kNN-based metrics.
        n_samples (int, optional): Number of samples for dimred metrics.
        use_memmap (bool): Whether to use memory-mapped files for large distance matrices.
        temp_dir (str, optional): Directory for temporary memmap files.
        chunk_size (int): Chunk size for batched distance computation.

    Returns:
        pd.DataFrame: Consolidated results with columns:
                      [Atlas Name, Task, Metric Name, Input 1, Input 2, Value, Time (s)]
    """
    import gc

    import scanpy as sc

    overall_start = time.time()

    # ────────────────────────────────────────────────────────────────────
    # 1. Load atlas
    # ────────────────────────────────────────────────────────────────────
    if adata is None and atlas_path is not None:
        if verbose:
            print(f"Loading atlas from: {atlas_path}")
        adata = sc.read_h5ad(atlas_path)
        if atlas_name is None:
            atlas_name = os.path.splitext(os.path.basename(atlas_path))[0]
    elif adata is not None:
        if atlas_name is None:
            atlas_name = "atlas"
    else:
        raise ValueError("Provide either `adata` (AnnData) or `atlas_path` (str).")

    if verbose:
        print(f"\n{'='*70}")
        print(f"  CheckAtlas — Unified Metrics Pipeline")
        print(f"  Atlas: {atlas_name}")
        print(f"  Cells: {adata.n_obs:,}  |  Genes: {adata.n_vars:,}")
        print(f"  Tasks: ", end="")
        tasks = []
        if run_annotation:
            tasks.append("Annotation")
        if run_clustering:
            tasks.append("Clustering")
        if run_dimred:
            tasks.append("Dimred")
        print(" → ".join(tasks))
        print(f"{'='*70}\n")

    # Set output directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)

    # Consolidated results
    all_results = []
    task_summary = {}

    # ────────────────────────────────────────────────────────────────────
    # 2. Run Annotation Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_annotation:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 1/3: Annotation Metrics")
            print(f"{'─'*50}")

        try:
            df_annot = cal_annot(
                adata,
                atlas_name=atlas_name,
                all=all_annot_metrics,
                file_dir=file_dir,
            )

            if not df_annot.empty:
                # Normalize columns to unified schema
                unified = pd.DataFrame(
                    {
                        "Atlas Name": df_annot["Atlas Name"],
                        "Task": "annotation",
                        "Metric Name": df_annot["Metric Name"],
                        "Input 1": df_annot["Reference/Input 1"],
                        "Input 2": df_annot["Prediction/Input 2"],
                        "Value": df_annot["Value"],
                        "Time (s)": df_annot.get("Time (s)", np.nan),
                    }
                )
                all_results.append(unified)

            task_elapsed = time.time() - task_start
            n_metrics = len(df_annot) if not df_annot.empty else 0
            task_summary["Annotation"] = {
                "count": n_metrics,
                "time": round(task_elapsed, 2),
            }
            if verbose:
                print(f"  ✓ Annotation: {n_metrics} results in {task_elapsed:.1f}s")

        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Annotation pipeline failed: {e}")
            task_summary["Annotation"] = {
                "count": 0,
                "time": round(task_elapsed, 2),
                "error": str(e),
            }
            if verbose:
                print(f"  ✗ Annotation failed: {e}")

    # ────────────────────────────────────────────────────────────────────
    # 3. Run Clustering Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_clustering:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 2/3: Clustering Metrics")
            print(f"{'─'*50}")

        try:
            df_cluster = cal_cluster(
                adata,
                atlas_name=atlas_name,
                all_metrics=all_cluster_metrics,
                file_dir=file_dir,
                n_jobs=n_jobs,
                verbose=verbose,
                seed=seed,
            )

            if not df_cluster.empty:
                unified = pd.DataFrame(
                    {
                        "Atlas Name": df_cluster["Atlas Name"],
                        "Task": "clustering",
                        "Metric Name": df_cluster["Metric Name"],
                        "Input 1": df_cluster["Embedding"],
                        "Input 2": df_cluster["Label Key"],
                        "Value": df_cluster["Value"],
                        "Time (s)": df_cluster["Time (s)"],
                    }
                )
                all_results.append(unified)

            task_elapsed = time.time() - task_start
            n_metrics = len(df_cluster) if not df_cluster.empty else 0
            task_summary["Clustering"] = {
                "count": n_metrics,
                "time": round(task_elapsed, 2),
            }
            if verbose:
                print(f"  ✓ Clustering: {n_metrics} results in {task_elapsed:.1f}s")

        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Clustering pipeline failed: {e}")
            task_summary["Clustering"] = {
                "count": 0,
                "time": round(task_elapsed, 2),
                "error": str(e),
            }
            if verbose:
                print(f"  ✗ Clustering failed: {e}")

    # ────────────────────────────────────────────────────────────────────
    # 4. Run Dimensionality Reduction Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_dimred:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 3/3: Dimensionality Reduction Metrics")
            print(f"{'─'*50}")

        try:
            metric_list = METRICS_DIMRED if all_dimred_metrics else None
            df_dimred = cal_dimred(
                adata,
                atlas_name=atlas_name,
                high_dim_key=high_dim_key,
                metric_list=metric_list,
                k_neighbors=k_neighbors,
                n_samples=n_samples,
                seed=seed,
                n_jobs=n_jobs,
                file_dir=file_dir,
                verbose=verbose,
            )

            if not df_dimred.empty:
                unified = pd.DataFrame(
                    {
                        "Atlas Name": atlas_name,
                        "Task": "dimred",
                        "Metric Name": df_dimred["Metric Name"],
                        "Input 1": df_dimred["Low Dim Key"],
                        "Input 2": df_dimred["High Dim Key"],
                        "Value": df_dimred["Value"],
                        "Time (s)": df_dimred["Time (s)"],
                    }
                )
                all_results.append(unified)

            task_elapsed = time.time() - task_start
            n_metrics = len(df_dimred) if not df_dimred.empty else 0
            task_summary["Dimred"] = {
                "count": n_metrics,
                "time": round(task_elapsed, 2),
            }
            if verbose:
                print(f"  ✓ Dimred: {n_metrics} results in {task_elapsed:.1f}s")

        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Dimred pipeline failed: {e}")
            task_summary["Dimred"] = {
                "count": 0,
                "time": round(task_elapsed, 2),
                "error": str(e),
            }
            if verbose:
                print(f"  ✗ Dimred failed: {e}")

    # ────────────────────────────────────────────────────────────────────
    # 5. Consolidate & Save
    # ────────────────────────────────────────────────────────────────────
    if all_results:
        df_all = pd.concat(all_results, ignore_index=True)
    else:
        df_all = pd.DataFrame(
            columns=[
                "Atlas Name",
                "Task",
                "Metric Name",
                "Input 1",
                "Input 2",
                "Value",
                "Time (s)",
            ]
        )

    # Save unified CSV
    if not df_all.empty:
        filename = os.path.join(file_dir, f"checkatlas_all_metrics_{atlas_name}.csv")
        df_all.to_csv(filename, index=False)
        if verbose:
            print(f"\nSaved unified results to: {filename}")
        logger.info(f"Saved unified metrics to {filename}")

    # ────────────────────────────────────────────────────────────────────
    # 6. Print Summary
    # ────────────────────────────────────────────────────────────────────
    overall_elapsed = time.time() - overall_start

    if verbose:
        print(f"\n{'='*70}")
        print(f"  CheckAtlas — Results Summary")
        print(f"{'='*70}")
        print(f"  {'Task':<20} {'Metrics':>10} {'Time (s)':>12} {'Status':>10}")
        print(f"  {'─'*52}")

        for task_name, info in task_summary.items():
            status = "✗ Error" if "error" in info else "✓ Done"
            print(
                f"  {task_name:<20} {info['count']:>10} {info['time']:>12.2f} {status:>10}"
            )

        print(f"  {'─'*52}")
        total_count = sum(info["count"] for info in task_summary.values())
        print(f"  {'TOTAL':<20} {total_count:>10} {overall_elapsed:>12.2f}")
        print(f"{'='*70}")

        if not df_all.empty:
            print(f"\n  Unique metrics computed: {df_all['Metric Name'].nunique()}")
            print(f"  Total measurements: {len(df_all)}")

            # Top 5 slowest metrics
            if "Time (s)" in df_all.columns:
                slowest = df_all.nlargest(5, "Time (s)")[
                    ["Task", "Metric Name", "Time (s)"]
                ]
                if not slowest.empty:
                    print(f"\n  Slowest metrics:")
                    for _, row in slowest.iterrows():
                        print(
                            f"    {row['Task']}/{row['Metric Name']}: {row['Time (s)']:.3f}s"
                        )

        print()

    gc.collect()
    return df_all


def annotation_to_num(annotation, ref_annotation):
    """
    Transforms the annotations from categorical to numerical

    Parameters
    ----------
    adata
    partition_key
    reference

    Returns
    -------

    """
    annotation = annotation.to_numpy()
    ref_annotation = ref_annotation.to_numpy()
    le = LabelEncoder()
    le.fit(annotation)
    annotation = le.transform(annotation)
    le2 = LabelEncoder()
    le2.fit(ref_annotation)
    ref_annotation = le2.transform(ref_annotation)
    return annotation, ref_annotation
