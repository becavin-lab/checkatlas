"""
Persistent preprocess context for CheckAtlas.

Stores column-detector results, kNN graph paths, distance-matrix paths,
neighbor graphs and other precomputed data so that downstream metric
processes can skip redundant computation.  The context is keyed by a
fingerprint that includes cell counts, embedding shapes and the source
file mtime so a stale cache is automatically invalidated.
"""

from __future__ import annotations

import logging
import os
import pickle
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger("checkatlas")


@dataclass
class PreprocessContext:
    """Container for all precomputed data a metric pipeline needs.

    Attributes
    ----------
    fingerprint : dict
        Immutable signature of the atlas state at precomputation time.
    ref_keys, pred_keys : list[str]
        Reference and predicted annotation columns from the column detector.
    embedding_keys : list[str]
        All ``.obsm`` keys (used by annotation and dimred metrics).
    cluster_embedding_keys : list[str]
        Embedding keys with > 3 dimensions (used by clustering and dimred).
    cluster_label_keys : list[str]
        Detected cluster-label columns.
    batch_keys : list[str]
        Detected batch / covariate columns.
    temp_parent_dir : str
        Root of the ``checkatlas_files/temp/`` directory.
    dimred_dir : str
        Cache directory for dimensionality-reduction precomputation.
    annotation_dir : str
        Cache directory for annotation precomputation.
    cluster_dir : str
        Cache directory for clustering precomputation.
    knn_paths : dict[str, str]
        ``{embedding: .npz_path}`` mapping for precomputed kNN graphs.
    neighbor_graphs : dict[str, dict]
        ``{embedding: {connectivities, distances, uns_entry, ...}}``
        mapping for precomputed ``sc.pp.neighbors`` output.
    """

    atlas_name: str = ""
    fingerprint: dict = field(default_factory=dict)

    # ── Column-detector results ──────────────────────────────────
    ref_keys: list[str] = field(default_factory=list)
    pred_keys: list[str] = field(default_factory=list)
    embedding_keys: list[str] = field(default_factory=list)
    cluster_embedding_keys: list[str] = field(default_factory=list)
    cluster_label_keys: list[str] = field(default_factory=list)
    batch_keys: list[str] = field(default_factory=list)

    # ── Filesystem roots ─────────────────────────────────────────
    temp_parent_dir: str = ""
    dimred_dir: str = ""
    annotation_dir: str = ""
    cluster_dir: str = ""

    # ── kNN graph paths ({embedding_key -> .npz path}) ───────────
    knn_paths: dict[str, str] = field(default_factory=dict)

    # ── Precomputed neighbour graphs for graph_connectivity ──────
    neighbor_graphs: dict[str, dict] = field(default_factory=dict)


def make_preprocess_fingerprint(
    adata,
    embedding_keys: list[str],
    cluster_label_keys: list[str],
    batch_keys: list[str],
    k_neighbors: int,
    source_path: str | None = None,
) -> dict:
    """Build a fingerprint uniquely identifying this precomputation.

    Reuses ``checkatlas.metrics._cache.compute_fingerprint`` and adds
    cluster/batch key sets for finer-grained invalidation.
    """
    from ._cache import compute_fingerprint

    emb_shapes = {}
    for k in embedding_keys:
        try:
            emb_shapes[k] = tuple(adata.obsm[k].shape)
        except Exception:
            emb_shapes[k] = (-1, -1)

    fp = compute_fingerprint(
        n_cells=adata.n_obs,
        n_features=getattr(adata, "n_vars", 0),
        embedding_keys=embedding_keys,
        embedding_shapes=emb_shapes,
        k_neighbors=k_neighbors,
        source_path=source_path,
    )
    fp["cluster_label_keys"] = sorted(cluster_label_keys)
    fp["batch_keys"] = sorted(batch_keys)
    return fp


def _non_empty(obj) -> bool:
    """Return True if *obj* is a non-empty collection or truthy scalar."""
    if obj is None:
        return False
    if isinstance(obj, (list, dict, tuple)):
        return len(obj) > 0
    return bool(obj) if isinstance(obj, (int, float)) else True


def fingerprint_match(cached: dict, current: dict) -> bool:
    """Return True if *cached* fingerprint matches *current*.

    Fields are only compared when **both** fingerprints carry a
    non-empty value for the key.  Empty lists / empty dicts in the
    *current* fingerprint are treated as "don't compare" — the
    caller is signalling that it doesn't know the value and is
    willing to trust the cache.
    """
    for key in ("n_cells", "n_features", "k_neighbors"):
        cv = cached.get(key)
        cr = current.get(key)
        if _non_empty(cv) and _non_empty(cr) and cv != cr:
            return False

    # Embedding keys / shapes — only compare if current provides them
    c_keys = current.get("embedding_keys", [])
    if _non_empty(c_keys):
        cached_keys = cached.get("embedding_keys", [])
        if sorted(cached_keys) != sorted(c_keys):
            return False

    c_shapes = current.get("embedding_shapes", {})
    if _non_empty(c_shapes):
        cached_shapes = cached.get("embedding_shapes", {})
        if set(cached_shapes.keys()) != set(c_shapes.keys()):
            return False
        for k in cached_shapes:
            if k in c_shapes and cached_shapes[k] != c_shapes[k]:
                return False

    for key in ("cluster_label_keys", "batch_keys"):
        cv = cached.get(key)
        cr = current.get(key)
        if _non_empty(cv) and _non_empty(cr) and sorted(cv) != sorted(cr):
            return False

    if "source_mtime" in cached and "source_mtime" in current:
        if cached["source_mtime"] != current["source_mtime"]:
            return False

    return True


def save_context(ctx: PreprocessContext) -> str:
    """Save *ctx* as a pickle alongside the precomputed artefact directories.

    Returns the path to the written pickle file.
    """
    ctx_dir = _context_dir(ctx.atlas_name, ctx.temp_parent_dir)
    os.makedirs(ctx_dir, exist_ok=True)
    path = os.path.join(ctx_dir, "context.pkl")
    with open(path, "wb") as f:
        pickle.dump(ctx, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("Preprocess context saved to %s", path)
    return path


def load_context(
    atlas_name: str,
    temp_parent_dir: str,
    current_fp: dict,
) -> PreprocessContext | None:
    """Load the saved context, or None on fingerprint mismatch / missing file.

    The fingerprint and paths are restored; ``adata`` is not restored.
    """
    path = os.path.join(_context_dir(atlas_name, temp_parent_dir), "context.pkl")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as f:
            ctx = pickle.load(f)
        if not fingerprint_match(ctx.fingerprint, current_fp):
            logger.info("Preprocess context fingerprint mismatch — cache invalidated")
            return None
        return ctx
    except Exception as exc:
        logger.warning("Failed to load preprocess context: %s", exc)
        return None


def _context_dir(atlas_name: str, temp_parent_dir: str) -> str:
    return os.path.join(temp_parent_dir, atlas_name, "preprocess")
