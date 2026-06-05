"""
Unified kNN abstraction for CheckAtlas metrics.

Replaces four independent, duplicate kNN implementations
(LISI, kBET, DBCVI, cal_dimred) with a single shared module.

Provides:
- ``NeighborResults`` dataclass — container for kNN indices + distances
- ``compute_neighbors()``  — auto-dispatches GPU (JAX) vs CPU (pynndescent)
"""

import logging
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional, Tuple

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from ._jax_utils import _GPU_AVAILABLE, _JAX_AVAILABLE, _get_ndarray, pdist_squareform

logger = logging.getLogger("checkatlas")


# ═══════════════════════════════════════════════════════════════════════
# NeighborResults dataclass
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class NeighborResults:
    """Container for precomputed k-nearest-neighbour results.

    All metrics that require kNN (LISI, kBET, PCR, DBCVI,
    graph_connectivity, silhouette, etc.) consume this object instead
    of computing their own kNN internally.

    Attributes
    ----------
    indices : np.ndarray
        kNN indices, shape ``(n_samples, n_neighbors)``.
    distances : np.ndarray
        kNN distances, shape ``(n_samples, n_neighbors)``.
    """

    indices: np.ndarray
    distances: np.ndarray

    def __post_init__(self):
        if self.indices.shape != self.distances.shape:
            raise ValueError(
                f"Shape mismatch: indices {self.indices.shape} != "
                f"distances {self.distances.shape}"
            )

    @property
    def n_samples(self) -> int:
        """Number of cells in the graph."""
        return self.indices.shape[0]

    @property
    def n_neighbors(self) -> int:
        """Number of neighbours per cell."""
        return self.indices.shape[1]

    def subset_neighbors(self, n: int) -> "NeighborResults":
        """Return a copy using only the first *n* neighbours."""
        if n > self.n_neighbors:
            raise ValueError(f"n={n} > available neighbours={self.n_neighbors}")
        return NeighborResults(
            indices=self.indices[:, :n],
            distances=self.distances[:, :n],
        )

    def to_knn_graph_distances(self) -> csr_matrix:
        """Sparse CSR adjacency matrix weighted by distance."""
        n_samples, n_neighbors = self.indices.shape
        rowptr = np.arange(0, n_samples * n_neighbors + 1, n_neighbors)
        return csr_matrix(
            (self.distances.ravel(), self.indices.ravel(), rowptr),
            shape=(n_samples, n_samples),
        )


# ═══════════════════════════════════════════════════════════════════════
# CPU backend — pynndescent (approximate kNN, linear-ish time)
# ═══════════════════════════════════════════════════════════════════════


def _pynndescent_knn(
    X: np.ndarray, n_neighbors: int, n_jobs: int = -1, random_state: int = 42
) -> NeighborResults:
    """Compute approximate kNN using pynndescent (CPU)."""
    try:
        from pynndescent import NNDescent

        n_trees = min(64, 5 + int(round((X.shape[0]) ** 0.5 / 20.0)))
        n_iters = max(5, int(round(np.log2(X.shape[0]))))
        max_candidates = 60

        index = NNDescent(
            X,
            n_neighbors=n_neighbors,
            random_state=random_state,
            low_memory=True,
            n_jobs=n_jobs,
            compressed=False,
            n_trees=n_trees,
            n_iters=n_iters,
            max_candidates=max_candidates,
        )
        knn_idx, knn_dist = index.neighbor_graph
    except ImportError:
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto", n_jobs=n_jobs)
        nbrs.fit(X)
        knn_dist, knn_idx = nbrs.kneighbors(X)

    return NeighborResults(indices=knn_idx, distances=knn_dist)


# ═══════════════════════════════════════════════════════════════════════
# GPU backend — JAX exact brute-force kNN (Phase 2)
# ═══════════════════════════════════════════════════════════════════════


def _jax_exact_knn(
    X: np.ndarray,
    n_neighbors: int,
    chunk_size: int = 4096,
    # JAX imports are module-local to keep them optional
) -> NeighborResults:
    """Compute **exact** kNN via JAX GPU brute-force distance computation.

    Unlike pynndescent (which builds an approximate index), this path
    computes *every* pairwise distance on the GPU (via a single matrix
    multiplication) and then selects the top-k via ``approx_min_k``.

    For N ≤ 30 000 cells the full N×N distance matrix fits in GPU
    memory in a single pass.  For larger atlases the computation is
    chunked in the query axis.

    **Falls back to pynndescent CPU when the input data size exceeds
    GPU memory** (e.g. 85k cells × 60k genes ≈ 20 GB on GPU).
    """
    import functools

    import jax
    import jax.numpy as jnp

    n, d = X.shape
    data_size_gb = (n * d * 4) / (1024**3)  # float32 bytes

    # ── GPU OOM guard: if input alone exceeds 10 GB, use CPU ────
    if data_size_gb > 10.0:
        from logging import getLogger

        getLogger("checkatlas").debug(
            "JAX kNN: input too large (%.1f GB). Falling back to pynndescent.",
            data_size_gb,
        )
        return _pynndescent_knn(X, n_neighbors, n_jobs=-1)

    # ── Small atlas: one-shot full matrix ──────────────────────────
    # 15k×15k float32 = 900 MB — easily fits in 40 GB A100
    n, d = X.shape
    if n <= 15000:
        db = jnp.asarray(X, dtype=jnp.float32)
        D = pdist_squareform(X)  # returns numpy via our jax_utils path
        D_jax = jnp.asarray(D, dtype=jnp.float32)
        knn_dists, knn_indices = jax.lax.approx_min_k(D_jax, k=n_neighbors)
        return NeighborResults(
            indices=_get_ndarray(knn_indices),
            distances=_get_ndarray(knn_dists),
        )

    # ── Medium atlas: chunked GPU (10k query rows at a time) ──────
    db = jnp.asarray(X, dtype=jnp.float32)
    all_dists = []
    all_idx = []

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        qy = db[start:end]
        # Compute (chunk × n) distance matrix on GPU
        qy_sq = jnp.sum(qy**2, axis=1)  # (chunk,)
        db_sq = jnp.sum(db**2, axis=1)  # (n,)
        dot = jnp.dot(qy, db.T)  # (chunk, n) ← GPU matmul
        D_chunk_sq = qy_sq[:, None] + db_sq[None, :] - 2 * dot
        D_chunk = jnp.sqrt(jnp.maximum(D_chunk_sq, 0.0))

        dist, idx = jax.lax.approx_min_k(D_chunk, k=n_neighbors)
        all_dists.append(_get_ndarray(dist))
        all_idx.append(_get_ndarray(idx))

    return NeighborResults(
        indices=np.concatenate(all_idx, axis=0),
        distances=np.concatenate(all_dists, axis=0),
    )


def _jax_streaming_knn(
    X: np.ndarray,
    n_neighbors: int,
    qchunk: int = 15000,
    rchunk: int = 10000,
    tri_memmap=None,
) -> NeighborResults:
    """GPU kNN for **any atlas size** — chunks both queries and references.

    Unlike :func:`_jax_exact_knn` which loads the full database on GPU,
    this function streams through the data in (query × reference) blocks.
    Each block computes a sub‑matrix of distances on GPU, then merges
    partial top‑k results into a running buffer.

    This is the scib‑metrics *external‑memory kNN* pattern, but with
    two‑axis chunking to avoid ``db = jnp.asarray(X)`` for high‑dim
    gene expression data (e.g. 85 k × 60 k = 20 GB on GPU).

    Parameters
    ----------
    X : shape (n, d)
    n_neighbors : int
    qchunk : int
        Query rows per GPU block.
    rchunk : int
        Reference rows per GPU block.
    tri_memmap : TriangularMatrix or None
        If provided, the upper‑triangle distance matrix is written
        to this memmap during the pass, fusing kNN + distance storage.

    Returns
    -------
    NeighborResults
    """
    import jax
    import jax.numpy as jnp

    n, d = X.shape
    k = n_neighbors

    _do_store = tri_memmap is not None

    running_dist = np.full((n, k), np.inf, dtype=np.float32)
    running_idx = np.full((n, k), -1, dtype=np.int32)

    n_qchunks = (n + qchunk - 1) // qchunk
    n_rchunks = (n + rchunk - 1) // rchunk
    total_blocks = n_qchunks * n_rchunks

    from tqdm import tqdm as _tqdm

    _desc = "GPU kNN + distances (streaming)" if _do_store else "GPU kNN (streaming)"

    with _tqdm(total=total_blocks, desc=_desc, disable=(total_blocks < 10)) as pbar:
        for qs in range(0, n, qchunk):
            qe = min(qs + qchunk, n)
            q = jnp.asarray(X[qs:qe], dtype=jnp.float32)
            q_sq = jnp.sum(q**2, axis=1)

            for rs in range(0, n, rchunk):
                re = min(rs + rchunk, n)
                r = jnp.asarray(X[rs:re], dtype=jnp.float32)
                r_sq = jnp.sum(r**2, axis=1)
                dot = jnp.dot(q, r.T)
                D_sq = q_sq[:, None] + r_sq[None, :] - 2 * dot
                D = jnp.sqrt(jnp.maximum(D_sq, 0.0))
                D_np = _get_ndarray(D)

                if _do_store:
                    from ._triangular import store_upper_triangle

                    store_upper_triangle(tri_memmap._data, D_np, qs, rs, n)

                ref_idx = np.arange(rs, re, dtype=np.int32)[None, :]
                ref_idx_broad = np.broadcast_to(ref_idx, (qe - qs, re - rs))

                merged_dist = np.concatenate([running_dist[qs:qe], D_np], axis=1)
                merged_idx = np.concatenate([running_idx[qs:qe], ref_idx_broad], axis=1)

                topk = np.argpartition(merged_dist, k, axis=1)[:, :k]
                running_dist[qs:qe] = np.take_along_axis(merged_dist, topk, axis=1)
                running_idx[qs:qe] = np.take_along_axis(merged_idx, topk, axis=1)

                pbar.update(1)

    # ── Final sort ───────────────────────────────────────────────
    sort_idx = np.argsort(running_dist, axis=1)
    running_dist = np.take_along_axis(running_dist, sort_idx, axis=1)
    running_idx = np.take_along_axis(running_idx, sort_idx, axis=1)

    return NeighborResults(indices=running_idx, distances=running_dist)


# ═══════════════════════════════════════════════════════════════════════
# Module-level cache (prevents recomputing kNN for identical input)
# ═══════════════════════════════════════════════════════════════════════

_NEIGHBORS_CACHE: dict = {}


def _make_cache_key(X: np.ndarray, n_neighbors: int) -> tuple:
    """Cache key from array shape + partial content hash."""
    flat = np.asarray(X).ravel()
    n = len(flat)
    chunks = [
        flat[:2000].tobytes(),
        flat[max(0, n // 2 - 1000) : n // 2 + 1000].tobytes(),
        flat[-2000:].tobytes() if n > 2000 else b"",
    ]
    return (X.shape, hash(b"".join(chunks)), n_neighbors)


def _clear_neighbors_cache():
    """Clear the kNN cache — call at end of pipeline to free memory."""
    _NEIGHBORS_CACHE.clear()


# ═══════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════


def compute_neighbors(
    X: np.ndarray,
    n_neighbors: int,
    backend: str = "auto",
    n_jobs: int = -1,
    use_cache: bool = True,
    tri_memmap=None,
) -> NeighborResults:
    """Compute k-nearest neighbours with automatic backend selection.

    Parameters
    ----------
    X : np.ndarray
        Data matrix, shape ``(n_samples, n_features)``.
    n_neighbors : int
        Number of neighbours to find (including self).
    backend : str
        ``"auto"`` — use GPU if JAX + CUDA available, else pynndescent.
        ``"jax"``  — force JAX GPU path (raises if unavailable).
        ``"pynndescent"`` — force CPU path.
    n_jobs : int
        Number of parallel jobs for CPU backend.  Ignored on GPU.
    use_cache : bool
        Cache kNN results per unique input to avoid recomputation
        across multiple metrics that share the same embedding.
    tri_memmap : TriangularMatrix or None
        If provided, the upper‑triangle distance matrix is written
        to this memmap during the streaming GPU pass (fusing
        kNN computation with distance storage).

    Returns
    -------
    NeighborResults
    """
    X_arr = np.asarray(X, dtype=np.float64)

    # ── Cache check ────────────────────────────────────────────────
    if use_cache:
        key = _make_cache_key(X_arr, n_neighbors)
        if key in _NEIGHBORS_CACHE:
            logger.debug("kNN cache hit (k=%d)", n_neighbors)
            return _NEIGHBORS_CACHE[key]

    # ── Backend dispatch ───────────────────────────────────────────
    if backend == "auto":
        use_jax = _JAX_AVAILABLE and _GPU_AVAILABLE
    elif backend == "jax":
        if not _JAX_AVAILABLE:
            raise RuntimeError("JAX backend requested but JAX is not installed.")
        use_jax = True
    elif backend == "pynndescent":
        use_jax = False
    else:
        raise ValueError(
            f"Unknown backend '{backend}'. Use 'auto', 'jax', or 'pynndescent'."
        )

    if use_jax:
        # ── Select GPU path based on input size ─────────────────
        data_gb = (X_arr.shape[0] * X_arr.shape[1] * 4) / (1024**3)
        if data_gb <= 10.0:
            result = _jax_exact_knn(X_arr, n_neighbors)
        else:
            # Large atlas — streaming GPU kNN (query×ref chunked)
            logger.debug(
                "Streaming GPU kNN (%.1f GB input, chunks q=%d r=%d)",
                data_gb,
                15000,
                10000,
            )
            result = _jax_streaming_knn(X_arr, n_neighbors, tri_memmap=tri_memmap)
    else:
        result = _pynndescent_knn(X_arr, n_neighbors, n_jobs=n_jobs)

    # ── Store in cache ─────────────────────────────────────────────
    if use_cache:
        _NEIGHBORS_CACHE[key] = result

    return result
