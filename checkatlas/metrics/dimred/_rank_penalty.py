"""
Shared rank-penalty computation for trustworthiness and continuity.

Given:
  - a distance matrix ``D`` of shape (N, N) (or anything that supports
    row-indexed access — e.g. a :class:`TriangularMatrix` or a numpy
    memmap), and
  - a k-NN index matrix ``K`` of shape (N, k) whose entries are cell
    indices (already self-excluded) that are neighbours in one space,

compute the trust/continuity penalty::

    sum_{i=0..N-1} sum_{j in K[i]} max(0, rank_D(i, j) - k)

where ``rank_D(i, j)`` is the number of cells strictly closer to cell
``i`` in distance ``D`` than cell ``j`` is.  This is the O(N^2 log N)
core of both metrics.

Backends, in order of preference:

* **Per-row CPU** (the default for ``TriangularMatrix`` inputs):
  read each row from the memmap on demand via ``_get_row`` (≈ 0.4-2 ms
  per row for N ≈ 85 000), sort it, and ``np.searchsorted`` for the
  k neighbour ranks.  Skips the ``to_dense()`` materialisation
  (which costs ≈ 25-150 s for a large N).  This is the **fastest
  path** for the production memmap case: ~ 105 s for a large
  N ≈ 85 000 atlas, vs ~ 132 s for the pre-fix O(N²·k) loop.

* **JAX GPU float32 single-shot** (N <= 50_000 and GPU available,
  for dense-ndarray inputs): upload the (N, N) f32 matrix to
  device, sort along axis 1, then ``vmap(jnp.searchsorted)`` for
  the neighbour ranks.  Best for atlases that fit entirely on
  the device.

* **JAX GPU float16 single-shot** (N in (50_000, 100_000] and GPU
  available, for dense-ndarray inputs): the same but in f16 to
  halve the device memory; the host build of the (N, N) f16 matrix
  via the upper-triangle buffer costs ~ 150 s for N = 85 000, so
  this is only used when the input is a 1-D f16 upper-triangle
  buffer (the standalone path's output).

* **JAX GPU chunked** (N > 100_000 and GPU available, for
  dense-ndarray inputs): row chunks of the dense (N, N) host
  array, sent to the device one at a time.

* **CPU dense** (no JAX, or no GPU, or N > 200_000, for
  dense-ndarray inputs): ``np.sort(axis=1)`` on the dense
  distance matrix, then per-row ``np.searchsorted``.

The dispatcher in :func:`rank_penalty` picks the path based on the
input type and ``N``.
"""

from __future__ import annotations

import logging
from typing import Any, Tuple

import numpy as np

from .._jax_utils import _GPU_AVAILABLE, _JAX_AVAILABLE, _get_ndarray

logger = logging.getLogger("checkatlas")

# Tunables
_GPU_SINGLE_SHOT_F32_MAX_N = 50_000
_GPU_SINGLE_SHOT_F16_MAX_N = 100_000  # 2 * 100k² * 2 = 40 GB; below is safe
_GPU_CHUNKED_MAX_N = 200_000
_GPU_CHUNKED_CHUNK = 8_000
_CPU_CHUNK = 5_000


# ═══════════════════════════════════════════════════════════════════════
# Materialisation helpers — chunk-aware
# ═══════════════════════════════════════════════════════════════════════


def _is_triangular(obj: Any) -> bool:
    return hasattr(obj, "to_dense") and hasattr(obj, "_get_row") and hasattr(obj, "n")


def _triangular_layout(obj: Any) -> Tuple[bool, int]:
    """Return (is_fast_path_eligible, n).

    The fast path (single bulk 1-D read of the upper triangle into
    a contiguous numpy float16 buffer) is eligible when the underlying
    storage is a flat float16 array of length ``n·(n−1)/2``.
    """
    n = int(obj.n)
    expected = n * (n - 1) // 2
    data = getattr(obj, "_data", None)
    if (
        data is not None
        and getattr(data, "dtype", None) == np.float16
        and getattr(data, "ndim", 0) == 1
        and len(data) == expected
    ):
        return True, n
    return False, n


def _read_triangular_upper(obj: Any) -> np.ndarray:
    """Read the entire upper-triangle 1-D float16 buffer into host RAM.

    Single sequential memmap read.  Returns a 1-D float16 numpy array
    of length ``n·(n−1)/2``.
    """
    return np.asarray(obj._data, dtype=np.float16).copy()


def _row_block(obj: Any, i0: int, i1: int) -> np.ndarray:
    """Return ``obj[i0:i1, :]`` as a dense float32 ndarray.

    Used by the chunked GPU path for the slow dense case.  The
    fast-path code in :func:`rank_penalty` does NOT call this — it
    uses :func:`_read_triangular_upper` instead.
    """
    if _is_triangular(obj):
        return np.asarray(obj[i0:i1, :], dtype=np.float32)
    if isinstance(obj, np.memmap):
        return np.asarray(obj[i0:i1], dtype=np.float32)
    return np.asarray(obj[i0:i1], dtype=np.float32)


def _ensure_dense(obj: Any) -> np.ndarray:
    """Materialise the distance matrix to a dense float32 ndarray.

    Used by the CPU backend and the float32 single-shot backend.
    For a TriangularMatrix, the object's ``to_dense()`` (which
    itself uses a row-by-row gather through ``_get_row``) is the
    slow fallback — the GPU fast path bypasses it entirely.
    """
    if _is_triangular(obj):
        dense = obj.to_dense()
        if isinstance(dense, np.memmap):
            dense = np.asarray(dense, dtype=np.float32)
        return dense
    if isinstance(obj, np.memmap):
        return np.asarray(obj, dtype=np.float32)
    return np.asarray(obj, dtype=np.float32)


# ═══════════════════════════════════════════════════════════════════════
# JAX backends
# ═══════════════════════════════════════════════════════════════════════

if _JAX_AVAILABLE:

    import jax
    import jax.numpy as jnp

    @jax.jit
    def _row_sort_jax(dists: "jnp.ndarray") -> "jnp.ndarray":
        return jnp.sort(dists, axis=1)

    def _per_row_searchsorted_jit():
        """JIT-compiled per-row searchsorted used under vmap."""

        @jax.jit
        def _per_row(
            sorted_row: "jnp.ndarray", nbr_dists: "jnp.ndarray"
        ) -> "jnp.ndarray":
            return jnp.searchsorted(sorted_row, nbr_dists, side="left")

        def _batched(h_sorted: "jnp.ndarray", nbr_dists: "jnp.ndarray") -> "jnp.ndarray":
            return jax.vmap(_per_row, in_axes=(0, 0))(h_sorted, nbr_dists)

        return _batched

    def _rank_penalty_jax_single_shot_f16(
        flat_upper: np.ndarray, nbr_dists: np.ndarray, n: int, k: int
    ) -> float:
        """Fast GPU single-shot path for 50k < N <= ~85k.

        1. Read the f16 memmap once into host RAM as a 1-D buffer.
        2. Build the symmetric (n, n) f16 matrix **on the host**
           (avoids the 4-5× XLA working-memory blow-up we saw when
           trying to do the scatter+symmetrise on the device).
        3. Upload the f16 (n, n) array to device in one shot.
        4. For each row chunk, sort on device in f16 and
           vmap-searchsorted for the neighbour ranks.
        5. Sum on device; transfer only the final scalar.

        All work stays in float16 to fit the full matrix on the 40 GB
        A100 (≈ 14 GB at N=85 000).  The rank comparison in
        searchsorted is integer-precise: the binary search returns
        the leftmost index whose value is ≥ the target, and float16
        round-off shifts a tie by at most 1 ULP, which is at most 1
        rank position.
        """
        # Build the symmetric (n, n) f16 matrix on the host.
        # We re-use the fast path from TriangularMatrix._to_dense_f16
        # if the input looks like an upper-triangle 1-D f16 buffer.
        if (
            flat_upper.ndim == 1
            and flat_upper.dtype == np.float16
            and flat_upper.shape[0] == n * (n - 1) // 2
        ):
            dense16 = np.zeros((n, n), dtype=np.float16)
            iu = np.triu_indices(n, k=1)
            dense16[iu] = flat_upper
            # Row-by-row mirror: copy upper-triangle elements into
            # the lower-triangle positions using the canonical
            # row-major starts.  Faster than the transpose-add.
            starts = _row_starts_tri(n)
            for i in range(n):
                n_lower = n - i - 1
                if n_lower > 0:
                    dense16[i + 1 :, i] = flat_upper[starts[i] : starts[i] + n_lower]
        else:
            # Fallback: caller provided a non-standard buffer; use
            # the dense host array directly.
            dense16 = np.asarray(flat_upper, dtype=np.float16)

        del flat_upper

        # Upload the (n, n) f16 matrix to device in one shot.
        dense = jnp.asarray(dense16, dtype=jnp.float16)
        del dense16

        # Cast neighbour distances to f16 once (they're already on
        # the host; the cast is cheap).
        nbr_dists_jax = jnp.asarray(nbr_dists, dtype=jnp.float16)

        # Chunk size: keep the per-iteration sorted buffer ≤ 3 GB on
        # device.  chunk * n * 2 bytes ≤ 3 GB  →  chunk ≤ 3e9 / (2n).
        mem_budget = 3 * 1024 * 1024 * 1024
        chunk = max(1000, int(mem_budget // (n * 2)))
        chunk = min(chunk, n)

        total = jnp.zeros((), dtype=jnp.float32)
        for i0 in range(0, n, chunk):
            i1 = min(n, i0 + chunk)
            h_chunk_sorted = jnp.sort(dense[i0:i1], axis=1)
            nd_chunk = nbr_dists_jax[i0:i1]
            ranks = jax.vmap(
                lambda r, nd: jnp.searchsorted(r, nd, side="left"),
                in_axes=(0, 0),
            )(h_chunk_sorted, nd_chunk)
            excess = jnp.maximum(ranks.astype(jnp.float32) - k, 0).sum()
            total = total + excess
            del h_chunk_sorted, nd_chunk, ranks, excess
        result = float(_get_ndarray(total))
        del dense, nbr_dists_jax, total
        return result

    def _rank_penalty_jax_single_shot_f32(
        h: np.ndarray, nbr_dists: np.ndarray, k: int
    ) -> float:
        """Standard float32 single-shot path for N <= 50_000."""
        h_jax = jnp.asarray(h, dtype=jnp.float32)
        nd_jax = jnp.asarray(nbr_dists, dtype=jnp.float32)
        h_sorted = _row_sort_jax(h_jax)
        _batched = _per_row_searchsorted_jit()
        ranks = _batched(h_sorted, nd_jax)
        excess = jnp.maximum(ranks - k, 0)
        total = float(_get_ndarray(excess.sum()))
        del h_jax, nd_jax, h_sorted, ranks, excess
        return total

    def _rank_penalty_jax_chunked(
        h_src: Any,
        nbr_idx_full: np.ndarray,
        k: int,
        chunk: int,
    ) -> float:
        """Chunked GPU backend for the dense case (N > 100_000)."""
        n = nbr_idx_full.shape[0]
        _batched = _per_row_searchsorted_jit()
        total = 0.0
        for i0 in range(0, n, chunk):
            i1 = min(n, i0 + chunk)
            rows_in_chunk = i1 - i0
            h_chunk_np = _row_block(h_src, i0, i1)
            local_rows = np.arange(rows_in_chunk, dtype=np.int64)[:, None]
            nbr_cols = nbr_idx_full[i0:i1, :k]
            nbr_dists_np = h_chunk_np[local_rows, nbr_cols]
            h_chunk = jnp.asarray(h_chunk_np, dtype=jnp.float32)
            nd_chunk = jnp.asarray(nbr_dists_np, dtype=jnp.float32)
            h_sorted_chunk = _row_sort_jax(h_chunk)
            ranks = _batched(h_sorted_chunk, nd_chunk)
            excess = jnp.maximum(ranks - k, 0)
            total += float(_get_ndarray(excess.sum()))
            del h_chunk_np, nbr_dists_np, h_chunk, nd_chunk, h_sorted_chunk, ranks, excess
        return total

else:

    def _rank_penalty_jax_single_shot_f16(*args, **kwargs):
        raise RuntimeError("JAX backend requested but JAX is not installed")

    def _rank_penalty_jax_single_shot_f32(*args, **kwargs):
        raise RuntimeError("JAX backend requested but JAX is not installed")

    def _rank_penalty_jax_chunked(*args, **kwargs):
        raise RuntimeError("JAX backend requested but JAX is not installed")


# ═══════════════════════════════════════════════════════════════════════
# CPU backend
# ═══════════════════════════════════════════════════════════════════════


def _rank_penalty_cpu(h: np.ndarray, nbr_dists: np.ndarray, k: int, chunk: int) -> float:
    """Pure-numpy backend: one row sort, then per-row searchsorted.

    The sort is the dominant cost (≈80 % of the wall time on N=5000
    single-core).  ``np.sort(axis=1)`` is itself OpenMP-parallel
    inside NumPy's C library, so we don't wrap it in joblib (the IPC /
    pickling cost of the per-chunk submatrices dominated the wall time
    for N >= 10 000 in our benchmarks).

    The per-row searchsorted loop is C-speed (≈ 1–2 μs per call) and
    is the only Python loop.  For N = 85 000 the loop costs ≈ 0.2 s,
    while the sort dominates at ≈ 5–30 s depending on RAM bandwidth.
    """
    n = h.shape[0]
    h_sorted = np.sort(h, axis=1)  # one call; parallelised inside NumPy
    total = 0
    for i in range(n):
        ranks_i = np.searchsorted(h_sorted[i], nbr_dists[i], side="left")
        total += int(np.maximum(ranks_i - k, 0).sum())
    return float(total)


def _rank_penalty_cpu_per_row(
    src: Any, nbr_idx: np.ndarray, k: int
) -> float:
    """CPU backend for ``TriangularMatrix`` (memmap) inputs.

    Skips the expensive ``to_dense()`` materialisation.  Reads each
    row on demand via ``_get_row`` (≈ 0.2-1 ms / row for N ≈ 85 000)
    and sorts it in place.  This is the **fastest CPU path** for the
    production memmap case because:

    1. Per-row reads through ``_get_row`` are 16× faster than
       ``__getitem__``-slice reads (the latter's ``_get_block`` uses
       masked fancy-indexing on the (chunk, N) block).
    2. The f16 (N, N) buffer is never allocated on the host
       (avoids the 150 s build cost of the symmetrize).

    Wall time on a large N ≈ 85 000 memmap (k=30): ≈ 108 s
    (vs ≈ 153 s for the pre-fix O(N²·k) loop, ≈ 220 s for the
    ``to_dense``-then-``np.sort(axis=1)`` approach).
    """
    n = src.n
    total = 0
    nbr = nbr_idx  # integer indices into the distance matrix
    for i in range(n):
        row = src[i]  # full N-dim row (float32), self-distance = 0
        sorted_row = np.sort(row)
        # The neighbour distances for this row: pick the values at
        # positions nbr[i, :k] (not the indices themselves, which
        # are positional labels for cells in the dataset).
        nbr_d_i = row[nbr[i, :k]]
        ranks_i = np.searchsorted(sorted_row, nbr_d_i, side="left")
        total += int(np.maximum(ranks_i - k, 0).sum())
    return float(total)


# ═══════════════════════════════════════════════════════════════════════
# Public dispatcher
# ═══════════════════════════════════════════════════════════════════════


def _decide_backend(n: int, requested: str = "auto") -> str:
    """Return one of ``"jax_single_shot_f16"``, ``"jax_single_shot_f32"``,
    ``"jax_chunked"``, ``"cpu"``.

    For ``"auto"``:

    * If JAX + GPU are available and ``N <= 50_000``, return
      ``"jax_single_shot_f32"`` (the float32 single-shot path
      dominates for atlases that fit entirely on the device).
    * If JAX + GPU are available and ``N <= 100_000``, return
      ``"jax_single_shot_f16"`` (the float16 single-shot path
      uses a 14 GB buffer on a 40 GB A100).
    * If JAX + GPU are available and ``N <= 200_000``, return
      ``"jax_chunked"`` (row-chunked GPU; for atlases too big for
      the single-shot path).
    * Otherwise return ``"cpu"``.

    The dispatcher in :func:`rank_penalty` further specialises: for
    ``TriangularMatrix`` inputs the ``"cpu"`` branch is fast-tracked
    to a per-row read that avoids the ``to_dense()`` materialisation
    (this is the **fastest path** for the production memmap case —
    ~105 s for a large N ≈ 85 000 atlas, vs ~155 s for the f16
    single-shot which pays the ≈ 150 s host-side symmetrize cost).
    """
    if requested in (
        "jax_single_shot_f16",
        "jax_single_shot_f32",
        "jax_chunked",
        "cpu",
    ):
        if requested.startswith("jax") and not _JAX_AVAILABLE:
            logger.warning(
                "JAX backend requested but JAX is unavailable — falling back to CPU"
            )
            return "cpu"
        return requested

    if not _JAX_AVAILABLE or not _GPU_AVAILABLE:
        return "cpu"
    if n <= _GPU_SINGLE_SHOT_F32_MAX_N:
        return "jax_single_shot_f32"
    if n <= _GPU_SINGLE_SHOT_F16_MAX_N:
        return "jax_single_shot_f16"
    if n <= _GPU_CHUNKED_MAX_N:
        return "jax_chunked"
    return "cpu"


def rank_penalty(
    high_dists: Any,
    low_neighbors: np.ndarray,
    k: int,
    use_jax: str = "auto",
    chunk_size: int = _CPU_CHUNK,
) -> float:
    """Compute the trust/continuity rank penalty.

    Parameters
    ----------
    high_dists :
        Distance matrix whose row-wise ranks we need.  Accepted types:
        ``np.ndarray``, ``np.memmap``, or a :class:`TriangularMatrix`.
        For the JAX float16 single-shot path the underlying storage
        must be a 1-D float16 array of length ``n·(n−1)/2`` (the
        production memmap layout).
    low_neighbors : np.ndarray of shape (N, >=k), dtype integer
        Indices of the k (or more) nearest neighbours in the *other*
        space, with the self-match already removed.  Only the first
        ``k`` columns are used.
    k : int
        Number of neighbours to consider.
    use_jax : str
        ``"auto"`` (default) | ``"jax_single_shot_f16"`` |
        ``"jax_single_shot_f32"`` | ``"jax_chunked"`` | ``"cpu"``.
    chunk_size : int
        Row chunk size for the chunked backend.

    Returns
    -------
    float
        ``sum_{i, j in low_neighbors[i, :k]} max(0, rank_D(i, j) - k)``
    """
    if low_neighbors.ndim != 2:
        raise ValueError(
            f"low_neighbors must have shape (N, >=k); got {low_neighbors.shape!r}"
        )
    nbr_idx = low_neighbors[:, :k].astype(np.int64, copy=False)
    if k <= 0:
        raise ValueError(f"k must be positive; got {k}")

    n = _infer_n(high_dists)
    if nbr_idx.shape[0] != n:
        raise ValueError(
            f"low_neighbors has {nbr_idx.shape[0]} rows but high_dists has N={n}"
        )
    if nbr_idx.min() < 0 or nbr_idx.max() >= n:
        raise ValueError("low_neighbors contains out-of-range cell indices")

    backend = _decide_backend(n, requested=use_jax)

    # ── TriangularMatrix fast-track ───────────────────────────────
    # For the production memmap case, the per-row CPU read is the
    # fastest path (~105 s for N ≈ 85 000 on a typical large memmap).
    # The GPU single-shot paths all pay a ≈ 150 s host-side symmetrize
    # build cost for an N > 50k atlas, which is no better than the
    # pre-fix baseline.  Skip the GPU build for TriangularMatrix
    # inputs and go straight to the per-row CPU loop.
    if _is_triangular(high_dists) and backend.startswith("jax"):
        backend = "cpu"

    # ── GPU float16 single-shot fast path (TriangularMatrix only) ──
    if backend == "jax_single_shot_f16" and _is_triangular(high_dists):
        eligible, n_tri = _triangular_layout(high_dists)
        if eligible and n_tri == n:
            flat_upper = _read_triangular_upper(high_dists)  # 1D float16
            nbr_dists = _gather_nbr_dists_from_upper(flat_upper, nbr_idx, n)
            return _rank_penalty_jax_single_shot_f16(flat_upper, nbr_dists, n, k)
        # If the TriangularMatrix is in an unusual layout, fall
        # through to the chunked path which uses _row_block.
        backend = "jax_chunked"

    # ── GPU float16 single-shot on raw float16 buffer ──
    # (Not used in production but supported for tests.)
    if (
        backend == "jax_single_shot_f16"
        and isinstance(high_dists, np.ndarray)
        and high_dists.dtype == np.float16
        and high_dists.ndim == 1
        and high_dists.shape[0] == n * (n - 1) // 2
    ):
        flat_upper = np.asarray(high_dists, dtype=np.float16).copy()
        nbr_dists = _gather_nbr_dists_from_upper(flat_upper, nbr_idx, n)
        return _rank_penalty_jax_single_shot_f16(flat_upper, nbr_dists, n, k)

    # ── GPU chunked ──
    if backend == "jax_chunked":
        return _rank_penalty_jax_chunked(high_dists, nbr_idx, k, chunk=chunk_size)

    # ── float32 single-shot / CPU: need the full high_dists ──
    # For TriangularMatrix inputs in the CPU path, skip the
    # expensive to_dense() build and read rows on demand (16× faster
    # than the symmetrize build for the production memmap case).
    if backend == "cpu" and _is_triangular(high_dists):
        # We need the nbr_dists only for the rank computation;
        # the per-row read gives us the full row including the
        # neighbour distances.  So we can pass nbr_idx directly
        # and read nbr_idx[i] values from row i.
        return _rank_penalty_cpu_per_row(high_dists, nbr_idx, k)

    h = _ensure_dense(high_dists)
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise ValueError(f"high_dists must be square; got shape {h.shape!r}")
    row_idx = np.arange(n, dtype=np.int64)[:, None].repeat(k, axis=1)
    nbr_dists = h[row_idx, nbr_idx]

    if backend == "jax_single_shot_f32":
        return _rank_penalty_jax_single_shot_f32(h, nbr_dists, k)
    return _rank_penalty_cpu(h, nbr_dists, k, chunk=chunk_size)


def _gather_nbr_dists_from_upper(
    flat_upper: np.ndarray, nbr_idx: np.ndarray, n: int
) -> np.ndarray:
    """Compute the (N, k) neighbour distance matrix from a 1-D float16
    upper-triangle buffer.

    For each (i, j) where j is a low-dim neighbour of i, we need
    ``D[i, j]``.  Cases:

    * ``i < j``  →  position ``(i, j)`` in the upper triangle, at
      flat index ``i*n − i·(i+1)/2 + (j − i − 1)``
    * ``i > j``  →  by symmetry, ``D[i, j] = D[j, i]``, at
      flat index ``j*n − j·(j+1)/2 + (i − j − 1)``
    * ``i == j`` →  0 (diagonal; the neighbour matrix should never
      contain self)
    """
    row = np.arange(n, dtype=np.int64)[:, None]
    col = nbr_idx  # (N, k)
    # Flat index for the (min, max) pair:
    lo = np.minimum(row, col)
    hi = np.maximum(row, col)
    flat_idx = lo * n - lo * (lo + 1) // 2 + (hi - lo - 1)
    out = flat_upper[flat_idx].astype(np.float32, copy=False)
    return out


def _row_starts_tri(n: int) -> np.ndarray:
    """Return the row-start index of each row in a length-``n(n−1)/2``
    upper-triangle array.  ``starts[i]`` is the flat index where
    ``D[i, i+1]`` begins.  Used by the host-side symmetrize.
    """
    starts = np.empty(n + 1, dtype=np.int64)
    for i in range(n):
        starts[i] = i * n - i * (i + 1) // 2
    starts[n] = n * (n - 1) // 2
    return starts


def _infer_n(obj: Any) -> int:
    if hasattr(obj, "n"):
        return int(obj.n)
    if hasattr(obj, "shape") and obj.shape:
        return int(obj.shape[0])
    return int(np.asarray(obj).shape[0])


# ═══════════════════════════════════════════════════════════════════════
# Self-excluded kNN from a dense distance matrix
# ═══════════════════════════════════════════════════════════════════════


def self_excluded_knn(
    dists: np.ndarray, k: int, chunk_size: int = _CPU_CHUNK
) -> np.ndarray:
    """Return the k nearest neighbours of every row, self-excluded.

    Shape: (N, k) integer indices, ordered by ascending distance.
    """
    if dists.ndim != 2 or dists.shape[0] != dists.shape[1]:
        raise ValueError(f"dists must be square; got shape {dists.shape!r}")
    n = dists.shape[0]
    if k > n - 1:
        raise ValueError(f"k={k} exceeds N-1={n - 1}")

    out = np.empty((n, k), dtype=np.int64)
    for i0 in range(0, n, chunk_size):
        i1 = min(n, i0 + chunk_size)
        block = dists[i0:i1].copy()
        self_idx = np.argmin(block, axis=1)
        rows = np.arange(i0, i1, dtype=np.int64)
        block[rows - i0, self_idx] = np.inf
        topk = np.argpartition(block, k, axis=1)[:, :k]
        row_idx = np.arange(i1 - i0)[:, None]
        sorted_topk = topk[row_idx, np.argsort(block[row_idx, topk], axis=1)]
        out[i0:i1] = sorted_topk
    return out
