"""
Persistent precomputation cache with cache validation for checkatlas.

Precomputed distance matrices (float16 upper‑triangle ``.tri``) and kNN
graphs (``.npz``) are stored per‑atlas under ``temp/<atlas_name>/<task>/``
and reused on subsequent runs when the atlas hasn't changed.

Cache invalidation checks: n_cells, embedding keys + shapes,
source file mtime, and k_neighbors.
"""

import json
import logging
import os
import time
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from ._triangular import TriangularMatrix

import numpy as np

logger = logging.getLogger("checkatlas")


# ═══════════════════════════════════════════════════════════════════════
# Fingerprint
# ═══════════════════════════════════════════════════════════════════════


def compute_fingerprint(
    n_cells: int,
    n_features: int,
    embedding_keys: list,
    embedding_shapes: dict,
    k_neighbors: int,
    source_path: Optional[str] = None,
) -> dict:
    """Create a fingerprint dict that uniquely identifies this precomputation.

    Parameters
    ----------
    n_cells, n_features : int
        Dimensions of the high‑dim data.
    embedding_keys : list[str]
        Keys in ``.obsm``.
    embedding_shapes : dict[str, tuple]
    k_neighbors : int
    source_path : str | None
        Path to the source .h5ad file (for mtime check).

    Returns
    -------
    dict — serialisable to JSON
    """
    fp = {
        "n_cells": n_cells,
        "n_features": n_features,
        "embedding_keys": sorted(embedding_keys),
        "embedding_shapes": {k: list(v) for k, v in sorted(embedding_shapes.items())},
        "k_neighbors": k_neighbors,
    }
    if source_path and os.path.exists(source_path):
        fp["source_mtime"] = os.path.getmtime(source_path)
        fp["source_path"] = source_path
    fp["created_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    return fp


def fingerprint_matches(cached: dict, current: dict) -> bool:
    """Return True if the cached fingerprint matches the current atlas."""
    critical = ["n_cells", "n_features", "embedding_keys", "k_neighbors"]
    for key in critical:
        if cached.get(key) != current.get(key):
            return False
    # Embedding shapes must match
    cached_shapes = cached.get("embedding_shapes", {})
    curr_shapes = current.get("embedding_shapes", {})
    if set(cached_shapes.keys()) != set(curr_shapes.keys()):
        return False
    for k in cached_shapes:
        if cached_shapes[k] != curr_shapes[k]:
            return False
    # Source file mtime
    if "source_mtime" in cached and "source_mtime" in current:
        if cached["source_mtime"] != current["source_mtime"]:
            return False
    return True


def read_fingerprint(cache_dir: str) -> Optional[dict]:
    """Read fingerprint.json from *cache_dir*, or None."""
    fp_path = os.path.join(cache_dir, "fingerprint.json")
    if not os.path.exists(fp_path):
        return None
    try:
        with open(fp_path) as f:
            return json.load(f)
    except Exception:
        return None


def write_fingerprint(cache_dir: str, fp: dict) -> None:
    """Write fingerprint.json."""
    os.makedirs(cache_dir, exist_ok=True)
    fp_path = os.path.join(cache_dir, "fingerprint.json")
    with open(fp_path, "w") as f:
        json.dump(fp, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════
# kNN save / load (tiny — .npz)
# ═══════════════════════════════════════════════════════════════════════


def save_knn(
    cache_dir: str, name: str, indices: np.ndarray, distances: np.ndarray
) -> None:
    """Save kNN indices + distances as ``{name}.npz``."""
    os.makedirs(cache_dir, exist_ok=True)
    np.savez_compressed(
        os.path.join(cache_dir, f"{name}.npz"),
        indices=indices.astype(np.int32),
        distances=distances.astype(np.float16),
    )


def load_knn(cache_dir: str, name: str) -> Optional[tuple]:
    """Load ``(indices, distances)`` from ``{name}.npz``, or None."""
    path = os.path.join(cache_dir, f"{name}.npz")
    if not os.path.exists(path):
        return None
    try:
        data = np.load(path)
        return data["indices"], data["distances"]
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════
# Triangular distance matrix save / load
# ═══════════════════════════════════════════════════════════════════════


def save_triangular(cache_dir: str, name: str, tri: "TriangularMatrix") -> str:
    """Move/rename a TriangularMatrix's backing file to a canonical name.

    The *tri* object must have been created with a temp filepath.
    This function renames it to ``{cache_dir}/{name}.tri`` and
    returns the new path.
    """
    from ._triangular import TriangularMatrix

    os.makedirs(cache_dir, exist_ok=True)
    dest = os.path.join(cache_dir, f"{name}.tri")

    src = tri._filepath
    if src and os.path.exists(src) and src != dest:
        tri._close_mmap()
        time.sleep(0.1)
        if os.path.exists(dest):
            os.remove(dest)
        os.rename(src, dest)
    elif src and src == dest:
        # Already at the right place
        tri.flush()
    return dest


def load_triangular(cache_dir: str, name: str, n: int) -> Optional["TriangularMatrix"]:
    """Open a saved TriangularMatrix for reading. Returns None if missing."""
    from ._triangular import TriangularMatrix

    path = os.path.join(cache_dir, f"{name}.tri")
    if not os.path.exists(path):
        return None
    try:
        return TriangularMatrix(n=n, filepath=path, mode="r")
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════
# Full dimred cache save / load
# ═══════════════════════════════════════════════════════════════════════


def save_dimred_cache(
    cache_dir: str,
    fp: dict,
    high_dim_dists,  # TriangularMatrix or numpy or None
    high_knn_dists: np.ndarray,
    high_knn_indices: np.ndarray,
    low_dim_data: dict,
    low_dim_keys: list,
    subsample_indices: Optional[np.ndarray] = None,
) -> None:
    """Save all precomputed dimred data to *cache_dir*.

    Distance matrices that are numpy arrays (in‑memory, from GPU one‑shot)
    are NOT saved — they're cheap to recompute.  Only TriangularMatrix
    memmaps are persisted.
    """
    os.makedirs(cache_dir, exist_ok=True)

    # Fingerprint
    write_fingerprint(cache_dir, fp)

    # High-dim — only save if it's a TriangularMatrix (persisted memmap)
    if high_dim_dists is not None and hasattr(high_dim_dists, "_filepath"):
        save_triangular(cache_dir, "high_dists", high_dim_dists)
    save_knn(cache_dir, "knn_high", high_knn_indices, high_knn_dists)

    # Subsample indices for large-atlas core-set protocol
    if subsample_indices is not None:
        np.save(
            os.path.join(cache_dir, "subsample_indices.npy"),
            subsample_indices.astype(np.int64),
        )
        fp["subsample_n"] = int(len(subsample_indices))

    # Low-dim per embedding
    for emb in low_dim_keys:
        info = low_dim_data.get(emb, {})
        safe_name = emb.replace("/", "_").replace(" ", "_")
        tri = info.get("dists")
        if tri is not None and hasattr(tri, "_filepath"):
            save_triangular(cache_dir, f"low_dists_{safe_name}", tri)
        if info.get("knn_indices") is not None:
            save_knn(
                cache_dir,
                f"knn_low_{safe_name}",
                info["knn_indices"],
                info["knn_dists"],
            )
        emb_subsample = info.get("subsample_indices")
        if emb_subsample is not None:
            np.save(
                os.path.join(cache_dir, f"subsample_indices_{safe_name}.npy"),
                emb_subsample.astype(np.int64),
            )

    logger.info("Precomputed dimred data cached for reuse at %s", cache_dir)


def load_dimred_cache(
    cache_dir: str,
    fp: dict,
    n_cells: int,
    low_dim_keys: list,
) -> Optional[dict]:
    """Try to load cached dimred precomputation. Returns None on miss.

    Partial cache is valid — distance matrices may be missing (for
    GPU one‑shot atlases where they live in RAM), but kNN data must
    be present.

    When the cache was produced by the core‑set subsampling protocol
    (large atlases > 300k cells) the distance matrices have a smaller
    *n* matching the subsample size read from ``subsample_indices.npy``.
    """
    cached_fp = read_fingerprint(cache_dir)
    if cached_fp is None or not fingerprint_matches(cached_fp, fp):
        return None

    # ── Subsample indices (optional — large-atlas core-set) ──
    subsample_indices: Optional[np.ndarray] = None
    sub_path = os.path.join(cache_dir, "subsample_indices.npy")
    if os.path.exists(sub_path):
        subsample_indices = np.load(sub_path).astype(np.int64)
        dist_n = int(len(subsample_indices))
    else:
        dist_n = n_cells

    # ── High-dim kNN (required) ──
    knn_high = load_knn(cache_dir, "knn_high")
    if knn_high is None:
        logger.debug("Cache miss: missing high‑dim kNN file")
        return None
    high_knn_indices, high_knn_dists = knn_high

    # ── High-dim distance matrix (optional — may be None for GPU one‑shot) ──
    high_dim_dists = load_triangular(cache_dir, "high_dists", n=dist_n)

    # ── Low-dim per embedding (kNN required, distances optional) ──
    low_dim_data = {}
    for emb in low_dim_keys:
        safe_name = emb.replace("/", "_").replace(" ", "_")
        knn_low = load_knn(cache_dir, f"knn_low_{safe_name}")
        if knn_low is None:
            logger.debug("Cache miss: missing low-dim kNN for %s", emb)
            return None
        low_indices, low_dists_knn = knn_low
        low_tri = load_triangular(cache_dir, f"low_dists_{safe_name}", n=dist_n)
        entry: dict = {
            "dists": low_tri,  # may be None
            "knn_indices": low_indices,
            "knn_dists": low_dists_knn,
        }
        emb_sub_path = os.path.join(cache_dir, f"subsample_indices_{safe_name}.npy")
        if os.path.exists(emb_sub_path):
            entry["subsample_indices"] = np.load(emb_sub_path).astype(np.int64)
        low_dim_data[emb] = entry

    logger.info("Cache HIT — reusing precomputed dimred data from %s", cache_dir)
    return {
        "high_dim_dists": high_dim_dists,
        "high_knn_indices": high_knn_indices,
        "high_knn_dists": high_knn_dists,
        "low_dim": low_dim_data,
        "subsample_indices": subsample_indices,
    }
