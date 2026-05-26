"""
Unified kNN abstraction for CheckAtlas metrics.

Replaces four independent, duplicate kNN implementations
(LISI, kBET, DBCVI, cal_dimred) with a single shared module.

Provides:
- ``NeighborResults`` dataclass — container for kNN indices + distances
- ``compute_neighbors()``  — auto-dispatches GPU (JAX) vs CPU (pynndescent)
"""

from dataclasses import dataclass, field
from functools import cached_property
import logging
from typing import Optional, Tuple

import numpy as np
from scipy.sparse import csr_matrix, coo_matrix

from ._jax_utils import _JAX_AVAILABLE, _GPU_AVAILABLE, pdist_squareform, _get_ndarray

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

        nbrs = NearestNeighbors(
            n_neighbors=n_neighbors, algorithm="auto", n_jobs=n_jobs
        )
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

    For N ≤ 15 000 cells the full N×N distance matrix fits in GPU
    memory in a single pass.  For larger atlases the computation is
    chunked in the query axis.
    """
    import jax
    import jax.numpy as jnp
    import functools

    n = X.shape[0]

    # ── Small atlas: one-shot full matrix ──────────────────────────
    # 10k×10k float32 = 400 MB — easily fits in 40 GB A100
    if n <= 15000:
        db = jnp.asarray(X, dtype=jnp.float32)
        D = pdist_squareform(X)  # returns numpy via our jax_utils path
        D_jax = jnp.asarray(D, dtype=jnp.float32)
        knn_dists, knn_indices = jax.lax.approx_min_k(D_jax, k=n_neighbors)
        return NeighborResults(
            indices=_get_ndarray(knn_indices),
            distances=_get_ndarray(knn_dists),
        )

    # ── Large atlas: chunked ───────────────────────────────────────
    db = jnp.asarray(X, dtype=jnp.float32)
    all_dists = []
    all_idx = []

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        qy = db[start:end]
        # Compute (chunk × n) distance matrix on GPU
        qy_sq = jnp.sum(qy**2, axis=1)           # (chunk,)
        db_sq = jnp.sum(db**2, axis=1)            # (n,)
        dot = jnp.dot(qy, db.T)                    # (chunk, n) ← GPU matmul
        D_chunk_sq = qy_sq[:, None] + db_sq[None, :] - 2 * dot
        D_chunk = jnp.sqrt(jnp.maximum(D_chunk_sq, 0.0))

        dist, idx = jax.lax.approx_min_k(D_chunk, k=n_neighbors)
        all_dists.append(_get_ndarray(dist))
        all_idx.append(_get_ndarray(idx))

    return NeighborResults(
        indices=np.concatenate(all_idx, axis=0),
        distances=np.concatenate(all_dists, axis=0),
    )


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
        flat[max(0, n // 2 - 1000): n // 2 + 1000].tobytes(),
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
        raise ValueError(f"Unknown backend '{backend}'. Use 'auto', 'jax', or 'pynndescent'.")

    if use_jax:
        result = _jax_exact_knn(X_arr, n_neighbors)
    else:
        result = _pynndescent_knn(X_arr, n_neighbors, n_jobs=n_jobs)

    # ── Store in cache ─────────────────────────────────────────────
    if use_cache:
        _NEIGHBORS_CACHE[key] = result

    return result
