import argparse
import logging
import os
import re
import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata import _io as _io
from sklearn.utils.fixes import _object_dtype_isnan
from tqdm import tqdm

try:
    from . import cellranger, check
    from .metrics import metrics
    from .metrics._cache import load_dimred_cache, save_dimred_cache, save_knn
    from .metrics._jax_utils import _GPU_AVAILABLE, _JAX_AVAILABLE, cdist_chunked
    from .metrics._neighbors import NeighborResults, compute_neighbors
    from .metrics._preprocess_context import (
        PreprocessContext,
        load_context,
        make_preprocess_fingerprint,
        save_context,
    )
    from .utils import files, folders
    from .utils import log_format
    from .utils.checkatlas_arguments import MAX_N_JOBS, cap_n_jobs
    from .utils.col_detector import CheckAtlasColumnDetector
except ImportError:
    from checkatlas.utils import files, folders
    from checkatlas.utils import log_format
    from checkatlas.utils.checkatlas_arguments import MAX_N_JOBS, cap_n_jobs
    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

"""
Atlas module for AnnData (Scanpy)
All the function to screen the atlases
"""

ANNDATA_TYPE = "AnnData"
ANNDATA_EXTENSION = ".h5ad"

CELLINDEX_HEADER = "cell_index"

# Default obs columns to search for clustering/annotation keys
OBS_CLUSTERS = [
    "leiden",
    "louvain",
    "celltype",
    "cell_type",
    "cell_annotation",
    "cluster",
    "seurat_clusters",
    "annotation",
]

# Default obs columns for QC display
OBS_QC = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "n_genes",
    "n_counts",
]

logger = logging.getLogger("checkatlas")

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
sc.settings.verbosity = 0


_PRECOMPUTE_PROCESSES = frozenset(
    (
        "preprocess",
        "metric_cluster",
        "metric_annot",
        "metric_batch_correction",
        "metric_dimred",
        "metric",
        "analyse",
    )
)


def _resolve_n_jobs(args) -> int:
    """Read ``args.n_jobs`` (or fall back to :data:`MAX_N_JOBS`) and cap it.

    This is the single entry point every other function in this module
    should use to obtain the n_jobs value, so the cap is enforced in
    exactly one place and cannot be bypassed.
    """
    raw = getattr(args, "n_jobs", None) if args is not None else None
    return cap_n_jobs(raw)


def _metrics_enabled(metric_list):
    """Return True if *metric_list* is non-empty and not ['none']."""
    return bool(metric_list) and metric_list != ["none"]


def _should_precompute(args) -> bool:
    """Return True if at least one metric category is non-['none']."""
    if args is None:
        return False
    return (
        _metrics_enabled(args.metric_cluster)
        or _metrics_enabled(args.metric_annot)
        or _metrics_enabled(getattr(args, "metric_batch_correction", None))
        or _metrics_enabled(args.metric_dimred)
    )


def _wants_task(args, task: str) -> bool:
    """Return True if *args* implies the given *task* precomputation is needed."""
    if args is None:
        return False
    process = getattr(args, "process", "preprocess")
    if process not in _PRECOMPUTE_PROCESSES:
        return False
    if task == "cluster":
        return _metrics_enabled(args.metric_cluster)
    elif task == "annot":
        return _metrics_enabled(args.metric_annot)
    elif task == "batch_correction":
        return _metrics_enabled(getattr(args, "metric_batch_correction", None))
    elif task == "dimred":
        return _metrics_enabled(args.metric_dimred)
    return False


def preprocess_atlas(atlas_info: dict, args=None) -> AnnData:
    """
    Read adata, clean it, and run task-specific precomputations
    (column detection, kNN graphs, distance matrices, etc.) based
    on which metric categories are requested in *args*.

    Precomputed artefacts are persisted under
    ``checkatlas_files/temp/<atlas>/<task>/`` so that downstream
    child processes (including Nextflow metric steps) can skip
    redundant computation.

    When *args* is ``None`` the function behaves exactly like the
    legacy version: read + clean only.

    Returns the cleaned AnnData object.
    """
    adata = read_atlas(atlas_info)
    adata = clean_scanpy_atlas(adata, atlas_info)

    if not _should_precompute(args):
        return adata

    run_cluster = _wants_task(args, "cluster")
    run_annot = _wants_task(args, "annot")
    run_batch = _wants_task(args, "batch_correction")
    run_dimred = _wants_task(args, "dimred")

    if not (run_cluster or run_annot or run_batch or run_dimred):
        return adata

    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    source_path = atlas_info.get(check.ATLAS_PATH_KEY, None)

    # Per-atlas timings (consumed by the summary table in __main__).
    # Stashed on the argparse Namespace so we don't have to change
    # every internal signature; this is a private convention used
    # only by the four create_metric_* functions and the four
    # _precompute_* helpers in this module.
    timings = getattr(args, "_checkatlas_timings", None)
    if timings is None:
        timings = {}
        try:
            args._checkatlas_timings = timings
        except (AttributeError, TypeError):
            # args is a plain object (tests) or unmodifiable; fall
            # back to a local dict so the per-precompute footer
            # still works.
            timings = {}

    # ── 1. Column detection (once) ──────────────────────────────────
    detector = CheckAtlasColumnDetector(adata)
    params = detector.detect_all_parameters()

    ref_keys = [c for c, _ in params["annotation"]["reference"]]
    pred_keys = [c for c, _ in params["annotation"]["predicted"]]
    cluster_label_keys = [c for c, _ in params["clustering"]["cluster_labels"]]
    batch_keys = [c for c, _ in params.get("batch", [])]
    if not batch_keys:
        batch_keys = [col for col in adata.obs.columns if "batch" in col.lower()]

    embedding_keys = []
    cluster_embedding_keys = []
    annotation_embedding_keys = []
    batch_correction_embedding_keys = []
    for emb, meta in params["clustering"]["embeddings"]:
        embedding_keys.append(emb)
        cluster_embedding_keys.append(emb)
        n_comp = meta.get("n_components", 0)
        if n_comp > 2:
            annotation_embedding_keys.append(emb)
            batch_correction_embedding_keys.append(emb)
    # Also add ALL .obsm keys for dimred/cluster use;
    # for annotation / batch_correction only include embeddings with > 2 components
    all_obsm_keys = adata.obsm_keys()
    for key in all_obsm_keys:
        if key not in embedding_keys:
            embedding_keys.append(key)
        if key not in cluster_embedding_keys:
            cluster_embedding_keys.append(key)
        if key not in annotation_embedding_keys:
            try:
                if adata.obsm[key].shape[1] > 2:
                    annotation_embedding_keys.append(key)
                    batch_correction_embedding_keys.append(key)
            except Exception:
                pass

    k_max = 90  # covers LISI (90) and kBET (25 via subset)
    fingerprint = make_preprocess_fingerprint(
        adata,
        embedding_keys=embedding_keys,
        cluster_label_keys=cluster_label_keys,
        batch_keys=batch_keys,
        k_neighbors=k_max,
        source_path=source_path,
        annotation_embedding_keys=annotation_embedding_keys,
        ref_keys=ref_keys,
        pred_keys=pred_keys,
        batch_correction_embedding_keys=batch_correction_embedding_keys,
    )

    # ── 2. Early exit: cached context still valid ──────────────────
    temp_parent = folders.get_folder(args.path, folders.TEMP)
    existing = load_context(atlas_name, temp_parent, fingerprint)
    if existing is not None:
        log_format.preprocess_header(
            logger,
            atlas_name,
            adata.n_obs,
            adata.n_vars,
            run_dimred=run_dimred,
            run_annot=run_annot,
            run_batch=run_batch,
            run_cluster=run_cluster,
            cache_hit=True,
        )
        timings["preprocess:cache_hit"] = 0.0
        log_format.preprocess_footer(
            logger,
            atlas_name=atlas_name,
            timings=timings,
            cache_hit=True,
        )
        return adata

    # ── 3. Build fresh context ─────────────────────────────────────
    ctx = PreprocessContext(
        atlas_name=atlas_name,
        fingerprint=fingerprint,
        ref_keys=ref_keys,
        pred_keys=pred_keys,
        embedding_keys=embedding_keys,
        annotation_embedding_keys=annotation_embedding_keys,
        batch_correction_embedding_keys=batch_correction_embedding_keys,
        cluster_embedding_keys=cluster_embedding_keys,
        cluster_label_keys=cluster_label_keys,
        batch_keys=batch_keys,
        temp_parent_dir=temp_parent,
        dimred_dir=os.path.join(temp_parent, atlas_name, folders.DIMRED),
        annotation_dir=os.path.join(temp_parent, atlas_name, folders.ANNOTATION),
        batch_correction_dir=os.path.join(
            temp_parent, atlas_name, folders.BATCH_CORRECTION
        ),
        cluster_dir=os.path.join(temp_parent, atlas_name, folders.CLUSTER),
    )
    for d in (
        ctx.dimred_dir,
        ctx.annotation_dir,
        ctx.batch_correction_dir,
        ctx.cluster_dir,
    ):
        os.makedirs(d, exist_ok=True)

    # ── 4. Task-specific precomputation ────────────────────────────
    n_jobs = _resolve_n_jobs(args)

    # Banner up front so the user can see the precompute plan
    # before the (potentially slow) artefact build starts.
    log_format.preprocess_header(
        logger,
        atlas_name,
        adata.n_obs,
        adata.n_vars,
        run_dimred=run_dimred,
        run_annot=run_annot,
        run_batch=run_batch,
        run_cluster=run_cluster,
        cache_hit=False,
    )

    if run_dimred:
        import time as _time
        _t0 = _time.time()
        _precompute_dimred(adata, ctx, k_neighbors=30, n_jobs=n_jobs)
        timings["preprocess:dimred"] = _time.time() - _t0

    if run_annot:
        import time as _time
        _t0 = _time.time()
        _precompute_annot(adata, ctx, k_neighbors=90, n_jobs=n_jobs)
        timings["preprocess:annot"] = _time.time() - _t0

    # batch_correction shares kNN / neighbour-graph precompute with
    # the annotation task. The shared kNN_paths + neighbour_graphs
    # are reused by the per-task engines; the only thing this helper
    # adds is the batch_correction_dir entry in the context (above)
    # so downstream code can find a per-task cache root.
    if run_batch:
        import time as _time
        _t0 = _time.time()
        _precompute_batch_correction(adata, ctx, k_neighbors=90, n_jobs=n_jobs)
        timings["preprocess:batch_correction"] = _time.time() - _t0

    if run_cluster:
        import time as _time
        _t0 = _time.time()
        _precompute_cluster(adata, ctx, k_neighbors=30, n_jobs=n_jobs)
        timings["preprocess:cluster"] = _time.time() - _t0

    # ── 5. Persist context ─────────────────────────────────────────
    save_context(ctx)

    log_format.preprocess_footer(
        logger,
        atlas_name=atlas_name,
        timings=timings,
        cache_hit=False,
    )

    return adata


# ═══════════════════════════════════════════════════════════════════════
# Precompute helpers — one per task category
# ═══════════════════════════════════════════════════════════════════════


def _precompute_dimred(
    adata: AnnData,
    ctx: PreprocessContext,
    k_neighbors: int = 90,
    n_jobs: int = -1,
    chunk_size: int = 1000,
) -> None:
    """Precompute high-dim and per-embedding dimred distance matrices + kNN.

    Writes to ``ctx.dimred_dir`` in the format expected by
    :func:`save_dimred_cache` so that :func:`cal_dimred` can reuse
    the results via :func:`load_dimred_cache`.
    """
    logger.info("Precomputing dimred: %d embedding(s)", len([k for k in ctx.embedding_keys if k != "X"]))

    n_obs = adata.n_obs
    sample_indices = np.arange(n_obs)
    n_cells = n_obs

    high_dim_key = "X"
    if high_dim_key == "X":
        high_dim_data = adata.X[sample_indices]
        if hasattr(high_dim_data, "toarray"):
            high_dim_data = high_dim_data.toarray()
    else:
        high_dim_data = adata.obsm[high_dim_key][sample_indices]

    high_n_features = high_dim_data.shape[1]
    use_memmap = n_cells > 10000

    # ── Fingerprint for cache validation ────────────────────────
    emb_to_eval = [k for k in ctx.embedding_keys if k != high_dim_key]
    if not emb_to_eval:
        logger.warning("No embeddings to precompute for dimred")
        return

    fp = make_preprocess_fingerprint(
        adata,
        embedding_keys=emb_to_eval,
        cluster_label_keys=ctx.cluster_label_keys,
        batch_keys=ctx.batch_keys,
        k_neighbors=k_neighbors,
        source_path=ctx.fingerprint.get("source_path"),
    )

    cached = load_dimred_cache(ctx.dimred_dir, fp, n_cells, emb_to_eval)
    if cached is not None:
        logger.info("Dimred cache exists — skipping recomputation")
        return

    # ── High-dim kNN ─────────────────────────────────────────────
    tqdm.write(f"  Computing high-dim kNN (k={k_neighbors + 1})...")
    try:
        high_knn = compute_neighbors(
            np.asarray(high_dim_data, dtype=np.float64),
            n_neighbors=k_neighbors + 1,
            backend="auto",
        )
        high_knn_dists = high_knn.distances
        high_knn_indices = high_knn.indices
    except Exception as exc:
        logger.warning("High-dim kNN failed: %s", exc)
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1, n_jobs=n_jobs).fit(
            high_dim_data
        )
        high_knn_dists, high_knn_indices = nbrs.kneighbors(high_dim_data)

    # ── High-dim distance matrix ──────────────────────────────────
    _DIST_METRICS = frozenset(
        ("kruskal_stress", "spearman_rho", "dCor", "trustworthiness", "continuity")
    )
    need_high_dists = bool(set(ctx.fingerprint.get("metric_dimred", _DIST_METRICS)) & _DIST_METRICS)

    high_dim_dists = None
    if need_high_dists:
        _use_gpu = _JAX_AVAILABLE and _GPU_AVAILABLE
        backend_label = "GPU (JAX)" if _use_gpu else f"CPU (n_jobs={n_jobs})"
        tqdm.write(
            "  Computing high-dim distance matrix"
            f" ({n_cells:,} x {n_cells:,})"
            " [shared with dimred + cluster metrics]..."
        )
        tqdm.write(f"  Backend: {backend_label}")
        if use_memmap:
            import uuid

            run_id = str(uuid.uuid4())[:8]
            from .metrics._triangular import TriangularMatrix, store_upper_triangle

            high_dists_path = os.path.join(ctx.dimred_dir, f"high_dists_{run_id}.tri")
            high_dim_dists = TriangularMatrix(
                n=n_cells, filepath=high_dists_path, mode="w+"
            )
            x_cluster_dists = None
            if ctx.cluster_dir:
                x_cluster_path = os.path.join(ctx.cluster_dir, "dist_X.tri")
                x_cluster_dists = TriangularMatrix(
                    n=n_cells, filepath=x_cluster_path, mode="w+"
                )

            if _use_gpu:
                _last_flush = [0]

                def _store_block(block, qs, rs):
                    store_upper_triangle(high_dim_dists._data, block, qs, rs, n_cells)
                    if x_cluster_dists is not None:
                        store_upper_triangle(x_cluster_dists._data, block, qs, rs, n_cells)
                    _last_flush[0] += 1
                    if _last_flush[0] % 20 == 0:
                        high_dim_dists.flush()
                        if x_cluster_dists is not None:
                            x_cluster_dists.flush()

                cdist_chunked(
                    np.asarray(high_dim_data, dtype=np.float32),
                    np.asarray(high_dim_data, dtype=np.float32),
                    metric="euclidean",
                    block_callback=_store_block,
                )
                high_dim_dists.flush()
                if x_cluster_dists is not None:
                    x_cluster_dists.flush()
            else:
                for i in tqdm(
                    range(0, n_cells, chunk_size),
                    desc="  High-dim distance matrix [CPU]",
                    unit="chunk",
                ):
                    end = min(i + chunk_size, n_cells)
                    from sklearn.metrics import pairwise_distances

                    block = pairwise_distances(
                        high_dim_data[i:end], high_dim_data, n_jobs=n_jobs
                    )
                    store_upper_triangle(high_dim_dists._data, block, i, 0, n_cells)
                    high_dim_dists.flush()
                    if x_cluster_dists is not None:
                        store_upper_triangle(x_cluster_dists._data, block, i, 0, n_cells)
                        x_cluster_dists.flush()
            if x_cluster_dists is not None:
                del x_cluster_dists
        else:
            from sklearn.metrics import pairwise_distances

            high_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
            for i in tqdm(
                range(0, n_cells, chunk_size),
                desc="  High-dim distance matrix [CPU]",
                unit="chunk",
            ):
                end = min(i + chunk_size, n_cells)
                high_dim_dists[i:end, :] = pairwise_distances(
                    high_dim_data[i:end], high_dim_data, n_jobs=n_jobs
                )
            if ctx.cluster_dir:
                np.save(
                    os.path.join(ctx.cluster_dir, "dist_X.npy"),
                    high_dim_dists.astype(np.float32),
                )

    # ── Per-embedding low-dim kNN + distances ────────────────────
    low_dim_data_cache = {}
    for emb_key in tqdm(emb_to_eval, desc="  Low-dim precompute (shared with cluster metrics)"):
        low_dim_data = adata.obsm[emb_key][sample_indices]
        low_n_cells = low_dim_data.shape[0]

        tqdm.write(f"  kNN for {emb_key}...")
        try:
            low_knn = compute_neighbors(
                np.asarray(low_dim_data, dtype=np.float64),
                n_neighbors=k_neighbors + 1,
                backend="auto",
            )
            low_knn_dists = low_knn.distances
            low_knn_indices = low_knn.indices
        except Exception:
            from sklearn.neighbors import NearestNeighbors

            nbrs_low = NearestNeighbors(
                n_neighbors=k_neighbors + 1, n_jobs=n_jobs
            ).fit(low_dim_data)
            low_knn_dists, low_knn_indices = nbrs_low.kneighbors(low_dim_data)

        low_dists = None
        if need_high_dists:
            if use_memmap:
                import uuid

                run_id = str(uuid.uuid4())[:8]
                safe_name = emb_key.replace("/", "_").replace(" ", "_")
                from .metrics._triangular import TriangularMatrix, store_upper_triangle

                low_dists_path = os.path.join(
                    ctx.dimred_dir, f"low_dists_{safe_name}_{run_id}.tri"
                )
                low_dists = TriangularMatrix(
                    n=low_n_cells, filepath=low_dists_path, mode="w+"
                )
                cluster_dists = None
                if ctx.cluster_dir:
                    cluster_path = os.path.join(
                        ctx.cluster_dir, f"dist_{safe_name}.tri"
                    )
                    cluster_dists = TriangularMatrix(
                        n=low_n_cells, filepath=cluster_path, mode="w+"
                    )

                if _use_gpu:
                    _last_flush_low = [0]

                    def _store_low_block(block, qs, rs):
                        store_upper_triangle(low_dists._data, block, qs, rs, low_n_cells)
                        if cluster_dists is not None:
                            store_upper_triangle(cluster_dists._data, block, qs, rs, low_n_cells)
                        _last_flush_low[0] += 1
                        if _last_flush_low[0] % 20 == 0:
                            low_dists.flush()
                            if cluster_dists is not None:
                                cluster_dists.flush()

                    tqdm.write(f"    Distances ({emb_key}) [GPU]...")
                    cdist_chunked(
                        np.asarray(low_dim_data, dtype=np.float32),
                        np.asarray(low_dim_data, dtype=np.float32),
                        metric="euclidean",
                        block_callback=_store_low_block,
                    )
                    low_dists.flush()
                    if cluster_dists is not None:
                        cluster_dists.flush()
                else:
                    for i in tqdm(
                        range(0, low_n_cells, chunk_size),
                        desc=f"    Distances ({emb_key}) [dimred + cluster] [CPU]",
                        unit="chunk",
                        leave=False,
                    ):
                        end = min(i + chunk_size, low_n_cells)
                        from sklearn.metrics import pairwise_distances

                        block = pairwise_distances(
                            low_dim_data[i:end], low_dim_data, n_jobs=n_jobs
                        )
                        store_upper_triangle(low_dists._data, block, i, 0, low_n_cells)
                        low_dists.flush()
                        if cluster_dists is not None:
                            store_upper_triangle(cluster_dists._data, block, i, 0, low_n_cells)
                            cluster_dists.flush()
                if cluster_dists is not None:
                    del cluster_dists
            else:
                from sklearn.metrics import pairwise_distances

                low_dists = np.zeros((low_n_cells, low_n_cells), dtype=np.float32)
                for i in tqdm(
                    range(0, low_n_cells, chunk_size),
                    desc=f"    Distances ({emb_key}) [dimred + cluster] [CPU]",
                    unit="chunk",
                    leave=False,
                ):
                    end = min(i + chunk_size, low_n_cells)
                    low_dists[i:end, :] = pairwise_distances(
                        low_dim_data[i:end], low_dim_data, n_jobs=n_jobs
                    )
                if ctx.cluster_dir:
                    safe_name = emb_key.replace("/", "_").replace(" ", "_")
                    np.save(
                        os.path.join(ctx.cluster_dir, f"dist_{safe_name}.npy"),
                        low_dists.astype(np.float32),
                    )

        low_dim_data_cache[emb_key] = {
            "dists": low_dists,
            "knn_indices": low_knn_indices,
            "knn_dists": low_knn_dists,
        }

    # ── Persist to cache ─────────────────────────────────────────
    try:
        save_dimred_cache(
            ctx.dimred_dir,
            fp,
            high_dim_dists=high_dim_dists,
            high_knn_dists=high_knn_dists,
            high_knn_indices=high_knn_indices,
            low_dim_data=low_dim_data_cache,
            low_dim_keys=emb_to_eval,
        )
        logger.info("Dimred precomputation saved to %s", ctx.dimred_dir)
    except Exception as exc:
        logger.warning("Failed to save dimred cache: %s", exc)

    # Clean up non-memmap arrays
    if high_dim_dists is not None and not hasattr(high_dim_dists, "_filepath"):
        del high_dim_dists
    import gc

    gc.collect()


def _precompute_annot(
    adata: AnnData,
    ctx: PreprocessContext,
    k_neighbors: int = 90,
    n_jobs: int = -1,
) -> None:
    """Precompute kNN graphs and neighbour graphs for annotation metrics.

    - Per-embedding kNN at k=90 (covers LISI's 90 and kBET's 25 via subset)
    - Per-embedding ``sc.pp.neighbors`` for ``graph_connectivity``

    Saves ``.npz`` files to ``ctx.annotation_dir`` and stores paths in
    ``ctx.knn_paths`` and neighbour-graph CSRs in ``ctx.neighbor_graphs``.
    """
    logger.info("Precomputing annotation: %d embedding(s)", len(ctx.annotation_embedding_keys))
    import gc

    safe_name = lambda s: s.replace("/", "_").replace(" ", "_")

    for emb in tqdm(ctx.annotation_embedding_keys, desc="  Annotation kNN + neighbors"):
        if emb not in adata.obsm:
            continue
        X_emb = np.asarray(adata.obsm[emb], dtype=np.float64)
        n_cells = X_emb.shape[0]
        k_eff = min(k_neighbors + 1, n_cells - 1)

        # ── kNN graph ──────────────────────────────────────────
        try:
            knn = compute_neighbors(X_emb, n_neighbors=k_eff, backend="auto")
        except Exception:
            from sklearn.neighbors import NearestNeighbors

            nbrs = NearestNeighbors(n_neighbors=k_eff, n_jobs=n_jobs).fit(X_emb)
            dists, idx = nbrs.kneighbors(X_emb)
            knn = NeighborResults(indices=idx, distances=dists)

        sn = safe_name(emb)
        npz_path = os.path.join(ctx.annotation_dir, f"knn_{sn}.npz")
        save_knn(
            ctx.annotation_dir, f"knn_{sn}", knn.indices, knn.distances
        )
        ctx.knn_paths[emb] = npz_path

        # ── Neighbour graph for graph_connectivity ──────────────
        key_added = f"neighbors_{emb}"
        try:
            sc.pp.neighbors(adata, use_rep=emb, key_added=key_added)
            neighbor_payload = {}
            if key_added in adata.uns:
                neighbor_payload["uns_entry"] = {
                    k: v
                    for k, v in adata.uns[key_added].items()
                    if k in ("params", "connectivities_key", "distances_key")
                }
                conn_key = adata.uns[key_added].get("connectivities_key", "connectivities")
                dist_key = adata.uns[key_added].get("distances_key", "distances")
                if conn_key in adata.obsp:
                    neighbor_payload["connectivities"] = adata.obsp[conn_key]
                if dist_key in adata.obsp:
                    neighbor_payload["distances"] = adata.obsp[dist_key]
            neighbor_payload["key_added"] = key_added
            ctx.neighbor_graphs[emb] = neighbor_payload
        except Exception as exc:
            logger.warning("Failed sc.pp.neighbors for %s: %s", emb, exc)

        gc.collect()

    logger.info("Annotation precomputation complete (%d kNN graphs)", len(ctx.knn_paths))


def _precompute_batch_correction(
    adata: AnnData,
    ctx: PreprocessContext,
    k_neighbors: int = 90,
    n_jobs: int = -1,
) -> None:
    """Precompute kNN / neighbour-graph artefacts for the batch-correction task.

    The batch-correction task (``cal_batch_correction``) consumes the
    same kNN graphs (k=90) and ``sc.pp.neighbors`` neighbour graphs
    that the annotation task already produces.  To keep the cache
    footprint low and to honour the *single source of truth*
    contract, this helper reuses the artefacts that
    ``_precompute_annot`` writes under ``ctx.annotation_dir``.

    Concretely:

    * The kNN ``.npz`` files are already populated in ``ctx.knn_paths``
      and live under ``ctx.annotation_dir``.  ``cal_batch_correction``
      loads them from there via the ``preprocess_context`` plumbing.
    * The neighbour graphs are already populated in
      ``ctx.neighbor_graphs``.  ``cal_batch_correction`` re-injects
      them into ``adata.uns`` / ``adata.obsp`` before invoking
      ``graph_connectivity.run``.

    This helper exists as a deliberate hook point: if a future
    batch-correction-specific precomputation is needed (e.g. an
    per-batch distance matrix), it goes here, and the per-atlas
    cache directory is ``ctx.batch_correction_dir`` (already created
    by :func:`preprocess_atlas`).
    """
    if not ctx.batch_correction_embedding_keys:
        ctx.batch_correction_embedding_keys = (
            ctx.annotation_embedding_keys or ctx.embedding_keys
        )

    logger.info(
        "Precomputing batch-correction: %d embedding(s) (shared with annotation)",
        len(ctx.batch_correction_embedding_keys),
    )

    if ctx.knn_paths:
        logger.info(
            "Reusing %d precomputed kNN graphs from annotation cache",
            len(ctx.knn_paths),
        )
    else:
        logger.info(
            "No annotation kNN cache found; batch-correction metrics "
            "will build their own kNN on the fly."
        )

    if ctx.neighbor_graphs:
        logger.info(
            "Reusing %d precomputed neighbour graphs from annotation cache",
            len(ctx.neighbor_graphs),
        )


def _precompute_cluster(
    adata: AnnData,
    ctx: PreprocessContext,
    k_neighbors: int = 90,
    n_jobs: int = -1,
    chunk_size: int = 1000,
) -> None:
    """Precompute per-embedding distance matrices for cluster metrics
    (primarily ``silhouette`` via ``precomputed_dists``).  Saves as
    upper-triangle float16 ``.tri`` files to ``ctx.cluster_dir``.
    """
    if not ctx.cluster_embedding_keys:
        ctx.cluster_embedding_keys = ctx.embedding_keys

    import gc

    from sklearn.metrics import pairwise_distances

    from .metrics._triangular import TriangularMatrix, store_upper_triangle

    safe_name = lambda s: s.replace("/", "_").replace(" ", "_")

    skippable = 0
    for emb in ctx.cluster_embedding_keys:
        tri_path = os.path.join(ctx.cluster_dir, f"dist_{safe_name(emb)}.tri")
        npy_path = tri_path.replace(".tri", ".npy")
        if os.path.exists(tri_path) or os.path.exists(npy_path):
            skippable += 1

    if skippable > 0:
        logger.info(
            "Reusing %d distance matrix(s) already computed by dimred precompute",
            skippable,
        )

    pending = len(ctx.cluster_embedding_keys) - skippable
    if pending == 0:
        logger.info(
            "All %d cluster distance matrix(s) already available"
            " from dimred precompute; skipping cluster precompute",
            len(ctx.cluster_embedding_keys),
        )
        return

    _use_gpu = _JAX_AVAILABLE and _GPU_AVAILABLE

    logger.info(
        "Precomputing cluster distances: %d embedding(s) (%d already cached from dimred)",
        len(ctx.cluster_embedding_keys),
        skippable,
    )

    for emb in tqdm(ctx.cluster_embedding_keys, desc="  Cluster distances"):
        if emb == "X":
            X_emb = adata.X
            if hasattr(X_emb, "toarray"):
                X_emb = X_emb.toarray()
        else:
            if emb not in adata.obsm:
                continue
            X_emb = np.asarray(adata.obsm[emb])

        n_cells = X_emb.shape[0]
        tri_path = os.path.join(ctx.cluster_dir, f"dist_{safe_name(emb)}.tri")
        npy_path = tri_path.replace(".tri", ".npy")

        # Skip if already computed (e.g. by dimred precompute)
        if os.path.exists(tri_path) or os.path.exists(npy_path):
            continue

        if n_cells > 10000:
            tri = TriangularMatrix(n=n_cells, filepath=tri_path, mode="w+")

            if _use_gpu:
                _last_flush_cluster = [0]

                def _store_cluster_block(block, qs, rs):
                    store_upper_triangle(tri._data, block, qs, rs, n_cells)
                    _last_flush_cluster[0] += 1
                    if _last_flush_cluster[0] % 20 == 0:
                        tri.flush()

                cdist_chunked(
                    np.asarray(X_emb, dtype=np.float32),
                    np.asarray(X_emb, dtype=np.float32),
                    metric="euclidean",
                    block_callback=_store_cluster_block,
                )
                tri.flush()
            else:
                for i in tqdm(
                    range(0, n_cells, chunk_size),
                    desc=f"    Distances ({emb}) [CPU]",
                    unit="chunk",
                    leave=False,
                ):
                    end = min(i + chunk_size, n_cells)
                    block = pairwise_distances(X_emb[i:end], X_emb, n_jobs=n_jobs)
                    store_upper_triangle(tri._data, block, i, 0, n_cells)
                tri.flush()
        else:
            dists = pairwise_distances(X_emb, n_jobs=n_jobs)
            np.save(npy_path, dists.astype(np.float32))

        gc.collect()

    logger.info("Cluster precomputation complete")


def detect_scanpy(atlas_path: str) -> dict:
    if atlas_path.endswith(ANNDATA_EXTENSION):
        atlas_info = dict()
        atlas_info[check.ATLAS_NAME_KEY] = os.path.splitext(os.path.basename(atlas_path))[
            0
        ]
        atlas_info[check.ATLAS_TYPE_KEY] = ANNDATA_TYPE
        atlas_info[check.ATLAS_EXTENSION_KEY] = ANNDATA_EXTENSION
        atlas_info[check.ATLAS_PATH_KEY] = atlas_path
        return atlas_info
    else:
        return dict()


def read_atlas(atlas_info: dict) -> AnnData:
    """
    Read Scanpy or Cellranger data : .h5ad or .h5

    Args:
        atlas_info (dict): info dict about the atlas

    Returns:
        AnnData: scanpy object from .h5ad
    """
    logger.info(
        f"Load {atlas_info[check.ATLAS_NAME_KEY]} "
        f"in {atlas_info[check.ATLAS_PATH_KEY]}"
    )
    try:
        if atlas_info[check.ATLAS_TYPE_KEY] == cellranger.CELLRANGER_TYPE_CURRENT:
            logger.debug(
                "Read Cellranger >= v3 results " f"{atlas_info[check.ATLAS_PATH_KEY]}"
            )
            adata = cellranger.read_cellranger_current(atlas_info)
        elif atlas_info[check.ATLAS_TYPE_KEY] == cellranger.CELLRANGER_TYPE_OBSOLETE:
            logger.debug(
                "Read Cellranger < v3 results " f"{atlas_info[check.ATLAS_PATH_KEY]}"
            )
            adata = cellranger.read_cellranger_obsolete(atlas_info)
        else:
            logger.debug(f"Read Scanpy file {atlas_info[check.ATLAS_PATH_KEY]}")
            adata = sc.read_h5ad(atlas_info[check.ATLAS_PATH_KEY])
        return adata
    except _io.utils.AnnDataReadError:
        logger.warning(
            "AnnDataReadError, cannot read: " f"{atlas_info[check.ATLAS_PATH_KEY]}"
        )
        return dict()


def clean_scanpy_atlas(adata: AnnData, atlas_info: dict) -> AnnData:
    """
    Clean the Scanpy object to be sure to get all information out of it

    - Make var names unique
    - Make var unique for Raw matrix
    - If OBS_CLUSTERS are present and in int32 -> be sure to
    transform them in categorical

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas

    Returns:
        AnnData: cleaned atlas
    """
    logger.debug(f"Clean scanpy: {atlas_info[check.ATLAS_NAME_KEY]}")
    # Make var names unique
    list_var = adata.var_names
    if len(set(list_var)) == len(list_var):
        logger.debug("Var names unique")
    else:
        logger.debug("Var names not unique, ran : adata.var_names_make_unique()")
        adata.var_names_make_unique()
        # Test a second time if it is unique (sometimes it helps)
        list_var = adata.var_names
        if len(set(list_var)) == len(list_var):
            logger.debug("Var names unique")
        else:
            logger.debug("Var names not unique, ran : adata.var_names_make_unique()")
            adata.var_names_make_unique()
            # If it is still not unique, create unique var_names "by hand"
            list_var = adata.var_names
            if len(set(list_var)) == len(list_var):
                logger.debug("Var names unique")
            else:
                logger.debug("Var names not unique, ran : adata.var_names_make_unique()")
                adata.var.index = [
                    x + "_" + str(i)
                    for i, x in zip(range(len(adata.var)), adata.var_names)
                ]
                list_var = adata.var_names
                if len(set(list_var)) == len(list_var):
                    logger.debug("Var names unique")
    # Make var unique for Raw matrix
    if adata.raw is not None:
        list_var = adata.raw.var_names
        if len(set(list_var)) == len(list_var):
            logger.debug("Var names for Raw unique, transform ")
        else:
            logger.debug("Var names for Raw not unique")
            adata.raw.var.index = [
                x + "_" + str(i)
                for i, x in zip(range(len(adata.raw.var)), adata.raw.var_names)
            ]
            list_var = adata.raw.var_names
            if len(set(list_var)) == len(list_var):
                logger.debug("Var names for Raw unique")

    # If OBS_CLUSTERS are present and in int32 -> be sure to
    # transform them in categorical
    for obs_key in adata.obs_keys():
        for obs_key_celltype in OBS_CLUSTERS:
            if obs_key_celltype in obs_key:
                if (
                    adata.obs[obs_key].dtype == np.int32
                    or adata.obs[obs_key].dtype == np.int64
                ):
                    adata.obs[obs_key] = pd.Categorical(adata.obs[obs_key])
    return adata


def get_viable_obs_qc(adata: AnnData, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_QC values
    Extract sorted obs_keys in same order then OBS_QC

    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obs_keys
    """
    obs_keys = list()
    for obs_key in adata.obs_keys():
        if obs_key in args.qc_display:
            obs_keys.append(obs_key)
    return obs_keys


def get_viable_obs_annot(adata: AnnData, args: argparse.Namespace) -> list:
    """
    Search in obs_keys a match to OBS_CLUSTERS values
    ! Remove obs_key with only one category !
    Extract sorted obs_keys in same order then OBS_CLUSTERS

    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obs_keys
    """
    obs_keys = list()
    # Get keys from OBS_CLUSTERS
    for obs_key in adata.obs_keys():
        for obs_key_celltype in args.obs_cluster:
            if obs_key_celltype in obs_key:
                if isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype):
                    obs_keys.append(obs_key)
    # Remove keys with only one category and no NaN in the array
    obs_keys_final = list()
    for obs_key in obs_keys:
        annotations = adata.obs[obs_key]
        if not _object_dtype_isnan(annotations).any():
            categories_temp = annotations.cat.categories
            # remove nan if found
            categories = categories_temp.dropna()
            if True in categories.isin(["nan"]):
                index = categories.get_loc("nan")
                categories = categories.delete(index)
            # Add obs_key with more than one category (with Nan removed)
            if len(categories) != 1:
                logger.debug(f"Add obs_key {obs_key} with cat {categories_temp}")
                obs_keys_final.append(obs_key)
    return sorted(obs_keys_final)


def get_viable_obsm(adata: AnnData, args: argparse.Namespace) -> list:
    """
    TO DO
    Search viable obsm for dimensionality reduction metric
    calc.
    ! No filter on osbm is appled for now !
    Args:
        adata (AnnData): atlas to analyse
        args (argparse.Namespace): list of arguments from checkatlas workflow

    Returns:
        list: obsm_keys
    """
    obsm_keys = list()
    # for obsm_key in adata.obsm_keys():
    #   if obsm_key in args.obsm_dimred:
    obsm_keys = adata.obsm_keys()
    logger.debug(f"Add obsm {obsm_keys}")
    return obsm_keys


def create_summary_table(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Create a table with all summarizing variables

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (str): info dict of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    atlas_type = atlas_info[check.ATLAS_TYPE_KEY]
    atlas_path = atlas_info[check.ATLAS_PATH_KEY]
    logger.debug(f"Create Summary table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.SUMMARY, check.TSV_EXTENSION, args.path
    )
    # Create summary table
    header = [
        "AtlasFileType",
        "NbCells",
        "NbGenes",
        "AnnData.raw",
        "AnnData.X",
        "File_extension",
        "File_path",
    ]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    df_summary["AtlasFileType"][atlas_name] = atlas_type
    df_summary["NbCells"][atlas_name] = adata.n_obs
    df_summary["NbGenes"][atlas_name] = adata.n_vars
    df_summary["AnnData.raw"][atlas_name] = adata.raw is not None
    df_summary["AnnData.X"][atlas_name] = adata.X is not None
    df_summary["File_extension"][atlas_name] = atlas_name
    df_summary["File_path"][atlas_name] = atlas_path
    df_summary.to_csv(csv_path, index=False, sep="\t")


def create_anndata_table(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Create an html table with all AnnData arguments
    The html code will make all elements of the table visible in MultiQC
    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info dict on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]

    logger.debug(f"Create Adata table for {atlas_name}")
    csv_path = files.get_file_path(
        atlas_name, folders.ANNDATA, check.TSV_EXTENSION, args.path
    )
    # Create AnnData table
    header = ["atlas_obs", "obsm", "var", "varm", "uns"]
    df_summary = pd.DataFrame(index=[atlas_name], columns=header)
    # html_element = "<span class=\"label label-primary\">"
    # new_line = ''
    # for value in list(adata.obs.columns):
    #     new_line += html_element + value + "</span><br>"
    #     print(new_line)
    df_summary["atlas_obs"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.obs.columns)) + "</code>"
    )
    df_summary["obsm"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.obsm_keys())) + "</code>"
    )
    df_summary["var"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.var_keys())) + "</code>"
    )
    df_summary["varm"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.varm_keys())) + "</code>"
    )
    df_summary["uns"][atlas_name] = (
        "<code>" + "</code><br><code>".join(list(adata.uns_keys())) + "</code>"
    )
    df_summary.to_csv(csv_path, index=False, quoting=False, sep="\t")


def create_qc_tables(adata: AnnData, atlas_info: dict, args: argparse.Namespace) -> None:
    """
    Display the atlas QC table
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    qc_path = files.get_file_path(atlas_name, folders.QC, check.TSV_EXTENSION, args.path)
    logger.debug(f"Create QC tables for {atlas_name}")
    qc_genes = []
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("mt")
        logger.debug(f"Mitochondrial genes in {atlas_name} for QC")
    else:
        logger.debug(f"No mitochondrial genes in {atlas_name} for QC")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    if len(adata.var[adata.var["mt"]]) != 0:
        qc_genes.append("ribo")
        logger.debug(f"Ribosomal genes in {atlas_name} for QC")
    else:
        logger.debug(f"No ribosomal genes in {atlas_name} for QC")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_genes,
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    df_annot = adata.obs[get_viable_obs_qc(adata, args)]
    # Rank cell by qc metric
    for header in df_annot.columns:
        if header != CELLINDEX_HEADER:
            new_header = f"cellrank_{header}"
            df_annot = df_annot.sort_values(header, ascending=False)
            df_annot.loc[:, [new_header]] = range(1, adata.n_obs + 1)

    # Sample QC table when more cells than args.plot_celllimit are present
    df_annot = atlas_sampling(df_annot, "QC", args)
    df_annot.loc[:, [CELLINDEX_HEADER]] = range(1, len(df_annot) + 1)
    df_annot.to_csv(qc_path, index=False, quoting=False, sep="\t")


def create_qc_plots(adata: AnnData, atlas_info: dict, args: argparse.Namespace) -> None:
    """
    Display the atlas QC plot
    Search for the OBS variable which correspond to the toal_RNA, total_UMI,
     MT_ratio, RT_ratio

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.settings.figdir = folders.get_workingdir(args.path)
    sc.set_figure_params(dpi_save=80)
    qc_path = os.sep + atlas_name + check.QC_FIG_EXTENSION
    logger.debug(f"Create QC violin plot for {atlas_name}")
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    sc.pl.violin(
        adata,
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
        ],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=qc_path,
    )


def create_umap_fig(adata: AnnData, atlas_info: dict, args: argparse.Namespace) -> None:
    """
    Display the UMAP of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.set_figure_params(dpi_save=150)
    # Search if umap reduction exists
    obsm_keys = get_viable_obsm(adata, args)
    r = re.compile(".*umap*.")
    obsm_umap_keys = list(filter(r.match, obsm_keys))
    if len(obsm_umap_keys) > 0:
        obsm_umap = obsm_umap_keys[0]
        logger.debug(f"Create UMAP figure for {atlas_name} with obsm={obsm_umap}")
        # Set the umap to display
        if isinstance(adata.obsm[obsm_umap], pd.DataFrame):
            # Transform to numpy if it is a pandas dataframe
            adata.obsm["X_umap"] = adata.obsm[obsm_umap].to_numpy()
        else:
            adata.obsm["X_umap"] = adata.obsm[obsm_umap]
        # Setting up figures directory
        sc.settings.figdir = folders.get_workingdir(args.path)
        umap_path = os.sep + atlas_name + check.UMAP_EXTENSION
        # Exporting umap
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.umap(adata, color=obs_keys[0], show=False, save=umap_path)
        else:
            sc.pl.umap(adata, show=False, save=umap_path)


def create_tsne_fig(adata: AnnData, atlas_info: dict, args: argparse.Namespace) -> None:
    """
    Display the TSNE of celltypes
    Search for the OBS variable which correspond to the celltype annotation

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info on the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    sc.set_figure_params(dpi_save=150)
    # Search if tsne reduction exists
    obsm_keys = get_viable_obsm(adata, args)
    r = re.compile(".*tsne*.")
    obsm_tsne_keys = list(filter(r.match, obsm_keys))
    if len(obsm_tsne_keys) > 0:
        obsm_tsne = obsm_tsne_keys[0]
        logger.debug(f"Create t-SNE figure for {atlas_name} with obsm={obsm_tsne}")
        # Set the t-sne to display
        if isinstance(adata.obsm[obsm_tsne], pd.DataFrame):
            # Transform to numpy if it is a pandas dataframe
            adata.obsm["X_tsne"] = adata.obsm[obsm_tsne].to_numpy()
        else:
            adata.obsm["X_tsne"] = adata.obsm[obsm_tsne]
        # Setting up figures directory
        sc.settings.figdir = sc.settings.figdir = folders.get_workingdir(args.path)
        tsne_path = os.sep + atlas_name + check.TSNE_EXTENSION
        # Exporting tsne
        obs_keys = get_viable_obs_annot(adata, args)
        if len(obs_keys) != 0:
            sc.pl.tsne(adata, color=obs_keys[0], show=False, save=tsne_path)
        else:
            sc.pl.tsne(adata, show=False, save=tsne_path)


def _try_load_context(atlas_info: dict, args) -> PreprocessContext | None:
    """Build a minimal fingerprint from the currently-loaded atlas and
    attempt to load a previously-saved :class:`PreprocessContext`.

    Returns ``None`` if the context file is missing, the fingerprint
    mismatches, or the atlas source file has changed.
    """
    try:
        atlas_name = atlas_info[check.ATLAS_NAME_KEY]
        temp_parent = folders.get_folder(args.path, folders.TEMP)
        source_path = atlas_info.get(check.ATLAS_PATH_KEY, None)

        # Minimal fingerprint — only source file identity.
        # fingerprint_match allows missing fields, so we only fill
        # the fields we can cheaply determine without a full reload.
        fp = {
            "n_cells": 0,
            "n_features": 0,
            "embedding_keys": [],
            "embedding_shapes": {},
            "k_neighbors": 90,
            "cluster_label_keys": [],
            "batch_keys": [],
        }
        if source_path and os.path.exists(source_path):
            fp["source_mtime"] = os.path.getmtime(source_path)
            fp["source_path"] = source_path
        return load_context(atlas_name, temp_parent, fp)
    except Exception:
        return None


def create_metric_cluster(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calculate all clustering metrics via the comprehensive ``cal_cluster``
    pipeline.  The pipeline auto-detects embeddings and cluster-label columns,
    runs every metric listed in ``METRICS_CLUST`` across all combinations,
    and writes results as a tab-separated file in the cluster folder.

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    if not args.metric_cluster:
        logger.info("Skipping clustering metrics (no metrics requested)")
        return

    if args.metric_cluster == ["none"]:
        logger.info("Skipping clustering metrics (--metric_cluster none)")
        return

    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    cluster_dir = folders.get_folder(args.path, folders.CLUSTER)

    log_format.task_header(
        logger,
        atlas_name=atlas_name,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        task="cluster",
        keys={
            "cluster labels": list(
                adata.obs.select_dtypes(include="category").columns
            ),
            "embeddings": list(adata.obsm_keys()),
        },
    )

    import time as _time
    _t0 = _time.time()
    preprocess_ctx = _try_load_context(atlas_info, args)

    df = metrics.cal_cluster(
        adata,
        atlas_name=atlas_name,
        metric_list=args.metric_cluster,
        file_dir=cluster_dir,
        n_jobs=_resolve_n_jobs(args),
        verbose=True,
        seed=42,
        preprocess_context=preprocess_ctx,
    )
    _elapsed = _time.time() - _t0
    _csv_path = None

    if not df.empty:
        csv_path = files.get_file_path(
            atlas_name,
            folders.CLUSTER,
            check.TSV_EXTENSION,
            args.path,
        )
        wide_df = metrics._pivot_cluster_to_wide(df, atlas_name)
        wide_df.to_csv(csv_path, index=False, sep="\t")
        _csv_path = csv_path
    else:
        logger.warning("No clustering metrics calculated for %s", atlas_name)

    _timings = getattr(args, "_checkatlas_timings", None)
    if _timings is not None:
        _timings["metric:cluster"] = (_elapsed, len(df))
    log_format.task_footer(
        logger,
        task="cluster",
        atlas_name=atlas_name,
        elapsed=_elapsed,
        n_rows=len(df),
        output_path=_csv_path,
    )


def create_metric_annot(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calculate all annotation metrics via the comprehensive ``cal_annot``
    pipeline (ref-vs-pred, embedding-based, batch/integration, graph).

    The pipeline auto-detects reference/predicted columns, embedding keys,
    and batch labels; runs every metric listed in ``METRICS_ANNOT``; and
    writes results as tab-separated files in the annotation folder.

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    if args.metric_annot == ["none"]:
        logger.info("Skipping annotation metrics (--metric_annot none)")
        return

    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    annotation_dir = folders.get_folder(args.path, folders.ANNOTATION)

    # Detect the keys the metric engine will use so the per-task
    # header reflects the same columns the cal_annot pipeline
    # picked.  This is a one-shot column-detector call: cheap
    # (sub-second on a 30k-cell atlas) and saves the user from
    # having to read the cal_annot source to know what ran.
    _detector = CheckAtlasColumnDetector(adata)
    _params = _detector.detect_all_parameters()
    _ref_keys = [c for c, _ in _params["annotation"]["reference"]]
    _pred_keys = [c for c, _ in _params["annotation"]["predicted"]]
    _emb_keys = [
        k for k, meta in _params["clustering"]["embeddings"]
        if meta.get("n_components", 0) > 2
    ]

    log_format.task_header(
        logger,
        atlas_name=atlas_name,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        task="annot",
        keys={
            "ref keys": _ref_keys,
            "pred keys": _pred_keys,
            "embeddings (>2 comp)": _emb_keys,
        },
    )

    import time as _time
    _t0 = _time.time()
    preprocess_ctx = _try_load_context(atlas_info, args)

    df = metrics.cal_annot(
        adata,
        atlas_name=atlas_name,
        metric_list=args.metric_annot,
        all=True,
        file_dir=annotation_dir,
        n_jobs=_resolve_n_jobs(args),
        verbose=True,
        preprocess_context=preprocess_ctx,
    )
    _elapsed = _time.time() - _t0
    _csv_path = None

    if not df.empty:
        csv_path = files.get_file_path(
            atlas_name,
            folders.ANNOTATION,
            check.TSV_EXTENSION,
            args.path,
        )
        wide_df = metrics._pivot_annot_to_wide(df, atlas_name)
        wide_df.to_csv(csv_path, index=False, sep="\t")
        _csv_path = csv_path
    else:
        logger.warning("No annotation metrics calculated for %s", atlas_name)

    _timings = getattr(args, "_checkatlas_timings", None)
    if _timings is not None:
        _timings["metric:annot"] = (_elapsed, len(df))
    log_format.task_footer(
        logger,
        task="annot",
        atlas_name=atlas_name,
        elapsed=_elapsed,
        n_rows=len(df),
        output_path=_csv_path,
    )


def create_metric_dimred(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calculate all dimred metrics via the comprehensive ``cal_dimred``
    pipeline.  The pipeline auto-detects all ``.obsm`` embedding keys,
    compares each against ``adata.X`` as the high‑dimensional reference,
    runs every metric listed by ``--metric_dimred``, and writes results as
    a tab‑separated file in the dimred folder compatible with MultiQC.

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    if args.metric_dimred == ["none"]:
        logger.info("Skipping dimred metrics (--metric_dimred none)")
        return

    atlas_name = atlas_info[check.ATLAS_NAME_KEY]

    log_format.task_header(
        logger,
        atlas_name=atlas_name,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        task="dimred",
        keys={
            "embeddings (all .obsm)": list(adata.obsm_keys()),
        },
    )

    import time as _time
    _t0 = _time.time()

    dimred_dir = folders.get_folder(args.path, folders.DIMRED)
    # Per-atlas persistent cache under temp/
    cache_dir = os.path.join(
        folders.get_folder(args.path, folders.TEMP), atlas_name, "dimred"
    )

    df = metrics.cal_dimred(
        adata,
        atlas_name=atlas_name,
        metric_list=args.metric_dimred,
        file_dir=cache_dir,
        use_cache=True,
        n_jobs=_resolve_n_jobs(args),
        verbose=True,
        seed=42,
    )
    _elapsed = _time.time() - _t0
    _csv_path = None

    if not df.empty:
        csv_path = files.get_file_path(
            atlas_name,
            folders.DIMRED,
            check.TSV_EXTENSION,
            args.path,
        )
        wide_df = metrics._pivot_dimred_to_wide(df, atlas_name)
        wide_df.to_csv(csv_path, index=False, sep="\t")
        _csv_path = csv_path
    else:
        logger.warning("No dimred metrics calculated for %s", atlas_name)

    _timings = getattr(args, "_checkatlas_timings", None)
    if _timings is not None:
        _timings["metric:dimred"] = (_elapsed, len(df))
    log_format.task_footer(
        logger,
        task="dimred",
        atlas_name=atlas_name,
        elapsed=_elapsed,
        n_rows=len(df),
        output_path=_csv_path,
    )


def create_metric_batch_correction(
    adata: AnnData, atlas_info: dict, args: argparse.Namespace
) -> None:
    """
    Calculate all batch-correction metrics via the comprehensive
    ``cal_batch_correction`` pipeline.  The pipeline auto-detects
    batch / donor / sample columns, embedding keys, and (for cLISI)
    cell-type labels; runs every metric listed by
    ``--metric_batch_correction``; and writes results as a
    tab-separated file in the batch_correction folder compatible
    with MultiQC.

    Args:
        adata (AnnData): atlas to analyse
        atlas_info (dict): info of the atlas
        args (argparse.Namespace): list of arguments from checkatlas workflow
    """
    metric_list = getattr(args, "metric_batch_correction", None)
    if metric_list == ["none"]:
        logger.info(
            "Skipping batch-correction metrics (--metric_batch_correction none)"
        )
        return

    atlas_name = atlas_info[check.ATLAS_NAME_KEY]
    batch_correction_dir = folders.get_folder(
        args.path, folders.BATCH_CORRECTION
    )

    # Detect the keys the batch-correction pipeline will consume
    # so the per-task header reflects the same columns the
    # cal_batch_correction engine picked.
    _detector = CheckAtlasColumnDetector(adata)
    _params = _detector.detect_all_parameters()
    _batch_keys = [c for c, _ in _params.get("batch", [])]
    if not _batch_keys:
        _batch_keys = [c for c in adata.obs.columns if "batch" in c.lower()]
    _ref_keys = [c for c, _ in _params["annotation"]["reference"]]
    _pred_keys = [c for c, _ in _params["annotation"]["predicted"]]
    _emb_keys = [
        k for k, meta in _params["clustering"]["embeddings"]
        if meta.get("n_components", 0) > 2
    ]

    log_format.task_header(
        logger,
        atlas_name=atlas_name,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        task="batch_correction",
        keys={
            "batch keys": _batch_keys,
            "ref keys (cLISI / graph_connectivity)": _ref_keys,
            "pred keys (cLISI / graph_connectivity)": _pred_keys,
            "embeddings (>2 comp)": _emb_keys,
        },
    )

    import time as _time
    _t0 = _time.time()
    preprocess_ctx = _try_load_context(atlas_info, args)

    df = metrics.cal_batch_correction(
        adata,
        atlas_name=atlas_name,
        metric_list=metric_list,
        all=True,
        file_dir=batch_correction_dir,
        n_jobs=_resolve_n_jobs(args),
        verbose=True,
        preprocess_context=preprocess_ctx,
    )
    _elapsed = _time.time() - _t0

    csv_path = files.get_file_path(
        atlas_name,
        folders.BATCH_CORRECTION,
        check.TSV_EXTENSION,
        args.path,
    )

    if not df.empty:
        wide_df = metrics._pivot_batch_correction_to_wide(df, atlas_name)
        wide_df.to_csv(csv_path, index=False, sep="\t")
    else:
        # Always write a header-only TSV so the per-atlas file
        # exists for downstream consumers (MultiQC, scfm).
        metric_list = metric_list or metrics.METRICS_BATCH
        header = ["Batch_Sample", "Embedding", "Batch Key"] + list(
            metric_list
        )
        pd.DataFrame(columns=header).to_csv(csv_path, index=False, sep="\t")
        logger.warning(
            "No batch-correction metrics calculated for %s; wrote "
            "empty TSV with header at %s",
            atlas_name,
            csv_path,
        )

    _timings = getattr(args, "_checkatlas_timings", None)
    if _timings is not None:
        _timings["metric:batch_correction"] = (_elapsed, len(df))
    log_format.task_footer(
        logger,
        task="batch_correction",
        atlas_name=atlas_name,
        elapsed=_elapsed,
        n_rows=len(df),
        output_path=csv_path,
    )


def atlas_sampling(
    df_annot: pd.DataFrame, type_df: str, args: argparse.Namespace
) -> pd.DataFrame:
    """
    If args.plot_celllimit != 0 and args.plot_celllimit < len(df_annot)
    The atlas qC table will be sampled for MultiQC

    Args:
        df_annot (pd.DataFrame): Table to sample
        type_df (str): type of table
        args (argparse.Namespace): arguments of checkatlas workflow

    Returns:
        pd.DataFrame: Sampled QC table
    """
    if args.plot_celllimit != 0 and args.plot_celllimit < len(df_annot):
        logger.debug(f"Sample {type_df} table with {len(df_annot)} cells")
        df_annot = df_annot.sample(args.plot_celllimit)
        logger.debug(f"{type_df} table sampled to {len(df_annot)} cells")
    return df_annot


# Public API functions for column detection


def col_annotation_ref(
    adata: AnnData, min_score: float = 0.5, return_all: bool = False
) -> Optional[str]:
    """
    Detect reference (ground truth) annotation column in AnnData object.

    This function uses intelligent semantic and statistical analysis to identify
    the most likely reference/ground truth cell type annotation column.

    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return list of all candidates with scores. Default: False

    Returns:
        str or List[Tuple[str, float]] or None:
            - If return_all=False: Best reference column name, or None if none found
            - If return_all=True: List of (column_name, score) tuples sorted by score

    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> ref_col = atlas.col_annotation_ref(adata)
        >>> print(f"Reference column: {ref_col}")
        >>>
        >>> # Get all candidates with scores
        >>> all_refs = atlas.col_annotation_ref(adata, return_all=True)
        >>> for col, score in all_refs:
        ...     print(f"{col}: {score:.3f}")
    """

    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(
        min_reference_score=min_score, min_predicted_score=0.3
    )

    ref_candidates = results["annotation"]["reference"]

    if return_all:
        return ref_candidates
    else:
        return ref_candidates[0][0] if ref_candidates else None


def col_annotation_pred(
    adata: AnnData,
    min_score: float = 0.5,
    return_all: bool = False,
    max_results: int = 5,
) -> Optional[List[str]]:
    """
    Detect predicted/cluster annotation columns in AnnData object.

    This function identifies columns containing cluster labels or automated
    cell type predictions (e.g., leiden, louvain, seurat_clusters, celltypist).

    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return with scores. Default: False
        max_results (int): Maximum number of columns to return. Default: 5

    Returns:
        List[str] or List[Tuple[str, float]] or None:
            - If return_all=False: List of column names sorted by confidence
            - If return_all=True: List of (column_name, score) tuples
            - None if no columns found

    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> pred_cols = atlas.col_annotation_pred(adata)
        >>> print(f"Predicted columns: {pred_cols}")
        >>>
        >>> # Get with scores
        >>> pred_with_scores = atlas.col_annotation_pred(adata, return_all=True)
        >>> for col, score in pred_with_scores:
        ...     print(f"{col}: {score:.3f}")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(
        min_reference_score=0.3, min_predicted_score=min_score
    )

    pred_candidates = results["annotation"]["predicted"][:max_results]

    if not pred_candidates:
        return None

    if return_all:
        return pred_candidates
    else:
        return [col for col, score in pred_candidates]


def col_cluster(
    adata: AnnData,
    min_score: float = 0.5,
    return_all: bool = False,
    max_results: int = 5,
) -> Optional[List[str]]:
    """
    Detect cluster label columns in AnnData object.

    Uses the dedicated cluster-label detector (Leiden, Louvain, k‑means,
    Seurat clusters, PhenoGraph, etc.) which applies semantic 70 % +
    statistical 30 % scoring tuned for algorithmic clustering outputs.

    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        min_score (float): Minimum confidence score threshold (0-1). Default: 0.5
        return_all (bool): If True, return with scores. Default: False
        max_results (int): Maximum number of columns to return. Default: 5

    Returns:
        List[str] or List[Tuple[str, float]] or None:
            - If return_all=False: List of column names sorted by confidence
            - If return_all=True: List of (column_name, score) tuples
            - None if no columns found

    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> cluster_cols = atlas.col_cluster(adata)
        >>> print(f"Cluster columns: {cluster_cols}")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters(min_cluster_score=min_score)

    clust_candidates = results["clustering"]["cluster_labels"][:max_results]

    if not clust_candidates:
        return None

    if return_all:
        return clust_candidates
    else:
        return [col for col, _score in clust_candidates]


def col_dimred(
    adata: AnnData, return_all: bool = False, max_results: int = 10
) -> Optional[List[str] | List[dict[str, Any]]]:
    """
    Detect dimensionality reduction representations in AnnData.obsm.

    This function identifies embedding keys like X_pca, X_umap, X_tsne, etc.

    Args:
        adata (AnnData): Scanpy AnnData object to analyze
        return_all (bool): If True, return with metadata. Default: False
        max_results (int): Maximum number of representations to return. Default: 10

    Returns:
        List[str] or List[Dict] or None:
            - If return_all=False: List of obsm keys (e.g., ['X_umap', 'X_pca'])
            - If return_all=True: List of dicts with 'key', 'shape', 'n_components'
            - None if no representations found

    Example:
        >>> import scanpy as sc
        >>> import checkatlas.atlas as atlas
        >>> adata = sc.read_h5ad("atlas.h5ad")
        >>> dimred_keys = atlas.col_dimred(adata)
        >>> print(f"Dimensionality reductions: {dimred_keys}")
        >>>
        >>> # Get with metadata
        >>> dimred_detailed = atlas.col_dimred(adata, return_all=True)
        >>> for emb in dimred_detailed:
        ...     print(f"{emb['key']}: {emb['n_components']} components")
    """
    detector = CheckAtlasColumnDetector(adata)
    results = detector.detect_all_parameters()

    embeddings = results["clustering"]["embeddings"][:max_results]

    if not embeddings:
        return None

    if return_all:
        return [
            {
                "key": key,
                "shape": meta["shape"],
                "n_components": meta["n_components"],
            }
            for key, meta in embeddings
        ]
    else:
        return [key for key, meta in embeddings]
