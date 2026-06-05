"""
Upper-triangle float16 symetric matrix storage for distance matrices.

Provides numpy-array-like read access to a full N×N symmetric matrix
that is physically stored as only the upper triangle (excluding diagonal)
at float16 precision — **75 % less space** than a full float32 matrix.

Space comparison for N = 85 233 cells:
    float32 full:     85 233² × 4  = 27.1 GB
    float16 full:     85 233² × 2  = 13.5 GB
    float16 tri (this):  N·(N−1)÷2 × 2 =  6.8 GB  ←

Usage::

    tri = TriangularMatrix(n_cells, filepath)  # or existing memmap
    row_i = tri[i]          # full row as float32 (reconstructed on read)
    upper = tri[i, i+1:]    # upper triangle of row i
    block = tri[100:200, :] # block reconstruction
    print(tri.shape)        # (85233, 85233)
"""

from typing import Literal

import numpy as np


class TriangularMatrix:
    """Upper-triangle float16 symetric matrix with full N×N read access.

    The matrix ``D`` is stored as a 1‑D ``float16`` memmap holding
    the **strict** upper‑triangle elements ``D[i, j]`` for ``i < j``
    in row‑major order.  The diagonal (``D[i, i] = 0``) and lower
    triangle are reconstructed on‑the‑fly via symmetry.

    Parameters
    ----------
    n : int
        Number of cells (matrix width).
    data : np.ndarray | np.memmap | None
        Pre‑existing ``float16`` array of length ``n·(n−1)//2``.
        If ``None``, **writable** memmap storage is created at *filepath*.
    filepath : str | None
        Path for memmap storage.  Required when ``data is None``.
    mode : str
        Memmap mode (``"w+"`` for write, ``"r"`` for read).  Ignored
        when *data* is supplied directly.
    """

    def __init__(
        self,
        n: int,
        data: np.ndarray | None = None,
        filepath: str | None = None,
        mode: Literal["r", "r+", "w+", "c"] = "w+",
    ):
        self.n = int(n)
        self._tri_size = self.n * (self.n - 1) // 2

        if data is not None:
            if len(data) != self._tri_size:
                raise ValueError(
                    f"data length {len(data):,} != expected "
                    f"{self._tri_size:,} for N={self.n}"
                )
            self._data = np.asarray(data, dtype=np.float16)
            self._filepath = None
        else:
            if filepath is None:
                raise ValueError("filepath required when data is None")
            self._data = np.memmap(
                filepath, dtype=np.dtype(np.float16), mode=mode, shape=(self._tri_size,)
            )
            self._filepath = filepath

        # Cache for row index computation
        self._row_starts = self._precompute_row_starts()

    @property
    def shape(self) -> tuple:
        return (self.n, self.n)

    @property
    def dtype(self):
        return np.float32  # public dtype — values converted on read

    # ── Row-start precomputation ────────────────────────────────────
    def _precompute_row_starts(self) -> np.ndarray:
        """Cache start indices of each row in the 1‑D upper‑triangle array.

        ``row_starts[i]`` is the index in ``_data`` where ``D[i, i+1]``
        begins.  The row ends at ``row_starts[i+1]`` (or ``_tri_size``).
        """
        row_starts = np.zeros(self.n + 1, dtype=np.int64)
        for i in range(self.n):
            row_starts[i] = i * self.n - i * (i + 1) // 2
        row_starts[self.n] = self._tri_size
        return row_starts

    def _upper_idx(self, i: int, j: int) -> int:
        """Flat index of element ``(i, j)`` with ``i < j``."""
        return int(self._row_starts[i] + (j - i - 1))

    # ── Full-row reconstruction ─────────────────────────────────────
    def _get_row(self, i: int) -> np.ndarray:
        """Reconstruct full row *i* as float32."""
        n = self.n
        out = np.zeros(n, dtype=np.float32)

        # Lower triangle: D[i, 0:i] = D[0:i, i] (symmetry)
        if i > 0:
            k = np.arange(i, dtype=np.int64)
            lower_idx = k.astype(np.int64) * n + i - k * (k + 1) // 2 - k - 1
            out[:i] = self._data[lower_idx].astype(np.float32)

        # Upper triangle: D[i, i+1:n] — contiguous in _data
        if i < n - 1:
            start = int(self._row_starts[i])
            end = int(self._row_starts[i + 1])
            out[i + 1 :] = self._data[start:end].astype(np.float32)

        # D[i, i] = 0 (already zero-initialised)
        return out

    # ── Block reconstruction ────────────────────────────────────────
    def _get_block(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """Reconstruct a block as float32.

        Parameters
        ----------
        rows : np.ndarray
            Integer array of row indices.
        cols : np.ndarray
            Integer array of column indices.

        Returns
        -------
        shape (len(rows), len(cols))
        """
        n_rows = len(rows)
        n_cols = len(cols)
        out = np.zeros((n_rows, n_cols), dtype=np.float32)

        row_arr = rows[:, None].astype(np.int64)  # (n_rows, 1)
        col_arr = cols[None, :].astype(np.int64)  # (1, n_cols)

        # ── Diagonal: where row == col → 0 (already zero) ─────────
        # ── Upper triangle: D[row, col] for row < col ────────────
        upper_mask = row_arr < col_arr
        if np.any(upper_mask):
            ri = np.broadcast_to(row_arr, (n_rows, n_cols))[upper_mask].astype(
                np.int64, copy=False
            )
            cj = np.broadcast_to(col_arr, (n_rows, n_cols))[upper_mask].astype(
                np.int64, copy=False
            )
            idx = self._row_starts[ri] + (cj - ri - 1)
            out[upper_mask] = self._data[idx].astype(np.float32)

        # ── Lower triangle: D[row, col] = D[col, row] ────────────
        lower_mask = row_arr > col_arr
        if np.any(lower_mask):
            ri = np.broadcast_to(row_arr, (n_rows, n_cols))[lower_mask].astype(
                np.int64, copy=False
            )
            cj = np.broadcast_to(col_arr, (n_rows, n_cols))[lower_mask].astype(
                np.int64, copy=False
            )
            # Swap: D[ri, cj] = D[cj, ri] since cj < ri
            idx = self._row_starts[cj] + (ri - cj - 1)
            out[lower_mask] = self._data[idx].astype(np.float32)

        return out

    # ── Element-wise pair extraction ──────────────────────────────────
    def _get_elements(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """Element-wise pair extraction (like ``arr[rows, cols]`` for 1‑D arrays).

        Returns a 1‑D float32 array of length ``len(rows)`` where
        ``out[i] = self[rows[i], cols[i]]``.
        """
        n = len(rows)
        out = np.empty(n, dtype=np.float32)

        ri = np.asarray(rows, dtype=np.int64)
        cj = np.asarray(cols, dtype=np.int64)

        diag_mask = ri == cj
        out[diag_mask] = np.float32(0.0)

        upper_mask = ri < cj
        if np.any(upper_mask):
            idx = self._row_starts[ri[upper_mask]] + (cj[upper_mask] - ri[upper_mask] - 1)
            out[upper_mask] = self._data[idx].astype(np.float32)

        lower_mask = ri > cj
        if np.any(lower_mask):
            idx = self._row_starts[cj[lower_mask]] + (ri[lower_mask] - cj[lower_mask] - 1)
            out[lower_mask] = self._data[idx].astype(np.float32)

        return out

    # ── Numpy indexing interface ─────────────────────────────────────
    def __getitem__(self, key):
        """Numpy‑compatible read‑only access.

        Supports:
        - ``tri[i]``      → full row *i* as float32 array
        - ``tri[i:j]``    → block of rows [i:j)
        - ``tri[i, j]``   → single element
        - ``tri[i, j:k]`` → slice of row *i* over columns [j:k)
        - ``tri[i:j, k:l]`` → block
        """
        if isinstance(key, tuple):
            # 2‑D access
            row_k, col_k = key

            # ── Scalar × scalar → single element ────────────────
            if isinstance(row_k, (int, np.integer)) and isinstance(
                col_k, (int, np.integer)
            ):
                i = int(row_k)
                j = int(col_k)
                if i == j:
                    return np.float32(0.0)
                elif i < j:
                    return np.float32(self._data[self._upper_idx(i, j)])
                else:
                    return np.float32(self._data[self._upper_idx(j, i)])

            # ── Scalar row × column slice ────────────────────────
            if isinstance(row_k, (int, np.integer)):
                i = int(row_k)
                cols = np.arange(self.n)[col_k]
                return self._get_block(np.array([i]), cols)[0]

            # ── Element‑wise pair indexing (1‑D arrays) ──────────
            if (
                isinstance(row_k, np.ndarray)
                and row_k.ndim == 1
                and isinstance(col_k, np.ndarray)
                and col_k.ndim == 1
                and row_k.dtype.kind not in ("b",)
            ):
                return self._get_elements(row_k, col_k)

            # ── Row slice × column slice → block ────────────────
            rows = np.arange(self.n)[row_k]
            cols = np.arange(self.n)[col_k]
            return self._get_block(rows, cols)

        # ── 1‑D access (row or row slice) ───────────────────────────
        if isinstance(key, (int, np.integer)):
            return self._get_row(int(key))

        # Slice of rows
        rows = np.arange(self.n)[key]
        return self._get_block(rows, np.arange(self.n))

    # ── Repr ─────────────────────────────────────────────────────────
    def __repr__(self) -> str:
        size_mb = self._tri_size * 2 / (1024**2)
        return (
            f"TriangularMatrix(n={self.n:,}, "
            f"tri_size={self._tri_size:,}, "
            f"dtype=float16, "
            f"disk={size_mb:.1f} MB)"
        )

    def to_dense(self) -> "np.ndarray":
        """Materialize full N×N matrix as float32 numpy array.

        This is expensive (O(N²) memory) — only use when a function
        requires a true numpy array (e.g. joblib pickling).
        The result is cached so subsequent calls are free.

        Implementation
        --------------
        The slow path (legacy, row-by-row gather through ``_get_row``)
        is kept as a fallback for any object that wasn't built with a
        contiguous memmap of the strict upper triangle (e.g. a unit
        test that constructed the object from an in-memory array).

        The fast path is taken when the underlying storage is a 1-D
        float16 array of length ``n·(n−1)/2`` (the normal case for
        the blood-atlas memmap).  It:

        1. Reads the entire memmap in one sequential pass;
        2. Allocates a dense float16 (N, N) buffer and fills the
           upper triangle by reshaping the 1-D read with strides;
        3. Symmetrises via transpose (lower triangle mirrored from
           upper);
        4. Converts to float32 in-place.

        For N = 85 000 this is ~ 25-40 s on a 40 GB machine,
        replacing the ~ 85-150 s of the row-by-row gather loop.
        """
        if hasattr(self, "_dense_cache"):
            return self._dense_cache

        dense16 = self._to_dense_f16()
        # 4. Promote to float32 (cached).
        self._dense_cache = dense16.astype(np.float32, copy=True)
        return self._dense_cache

    def to_dense_f16(self) -> "np.ndarray":
        """Materialize full N×N matrix as float16 numpy array.

        Same as :meth:`to_dense` but skips the float32 promotion.  Use
        this when the consumer is the GPU fast path
        (:func:`checkatlas.metrics.dimred._rank_penalty.rank_penalty`)
        which can upload a float16 (N, N) array directly to device —
        halving the host and device memory cost compared to the
        float32 version.

        Returns the full symmetric float16 (N, N) array (zeros on
        the diagonal, the lower triangle mirrored from the upper).
        """
        if hasattr(self, "_dense_f16_cache"):
            return self._dense_f16_cache
        self._dense_f16_cache = self._to_dense_f16()
        return self._dense_f16_cache

    def _to_dense_f16(self) -> "np.ndarray":
        """Build the symmetric float16 (N, N) array from the upper
        triangle.  Used by both ``to_dense`` and ``to_dense_f16``.

        Fast path (the memmap case): single sequential read of the
        1-D upper-triangle buffer, scatter to upper triangle, then
        row-by-row mirror into the lower triangle using the
        pre-computed ``_row_starts``.  This is ~2.5× faster than the
        ``dense += dense.T`` (transpose + add) approach because the
        transpose materialises a full copy.

        Slow path: per-row gather via ``_get_row`` for any object that
        isn't backed by a contiguous 1-D float16 buffer of the right
        length.
        """
        expected = self.n * (self.n - 1) // 2
        data = self._data
        if (
            data is not None
            and getattr(data, "dtype", None) == np.float16
            and getattr(data, "ndim", 0) == 1
            and len(data) == expected
        ):
            try:
                # 1. Single sequential read of the entire memmap.
                flat = np.asarray(data, dtype=np.float16).copy()
                # 2. Dense float16 (N, N) buffer, zero-initialised.
                dense16 = np.zeros((self.n, self.n), dtype=np.float16)
                # 3. Scatter the flat read into the upper triangle.
                iu = np.triu_indices(self.n, k=1)
                if iu[0].shape[0] == flat.shape[0]:
                    dense16[iu] = flat
                else:
                    raise ValueError("flat length mismatch")
                # 4. Symmetrise via row-by-row mirror: for each row i,
                #    copy the upper-triangle elements (already in
                #    ``dense16[i, i+1:]``) into the lower-triangle
                #    positions (``dense16[i+1:, i]``).  This is just
                #    a contiguous slice assignment per row, so it
                #    bypasses the costly transpose + add.
                starts = self._row_starts
                for i in range(self.n):
                    n_lower = self.n - i - 1
                    if n_lower > 0:
                        dense16[i + 1 :, i] = flat[starts[i] : starts[i] + n_lower]
                return dense16
            except Exception:
                pass  # fall through to slow path

        # ── Slow path: row-by-row gather (legacy) ──
        dense16 = np.zeros((self.n, self.n), dtype=np.float16)
        for i in range(self.n):
            dense16[i] = self._get_row(i)
        return dense16

    # ── Memmap flush / close ─────────────────────────────────────────
    def flush(self):
        """Flush memmap to disk (no‑op for in‑memory arrays)."""
        if hasattr(self._data, "flush"):
            self._data.flush()

    def _close_mmap(self):
        """Close the underlying mmap file handle."""
        if self._filepath is not None:
            import mmap as _mmap_mod

            attr = getattr(self._data, "_mmap", None)
            if attr is not None and isinstance(attr, _mmap_mod.mmap):
                try:
                    attr.close()
                except Exception:
                    pass


def save_triangular(
    filepath: str, n: int, upper_triangle: np.ndarray
) -> TriangularMatrix:
    """Save an upper‑triangle float16 array as a memmap and return wrapper.

    Parameters
    ----------
    filepath : str
    n : int
    upper_triangle : shape (n*(n-1)//2,)
    """
    tri = TriangularMatrix(n=n, filepath=filepath, mode="w+")
    tri._data[:] = np.asarray(upper_triangle, dtype=np.float16)
    tri.flush()
    return tri


def load_triangular(filepath: str, n: int) -> TriangularMatrix:
    """Open an existing triangular memmap for read‑only access."""
    return TriangularMatrix(n=n, filepath=filepath, mode="r")


def store_upper_triangle(tri_data, block, row_start, col_start, n_total):
    """Store the upper‑triangle portion of a distance block into a 1‑D triangular array.

    Parameters
    ----------
    tri_data : np.ndarray
        1‑D float16 memmap, length ``n_total·(n_total−1)//2``.
    block : np.ndarray
        Sub‑matrix of pairwise distances, shape ``(qsize, rsize)``.
    row_start : int
        Global row offset of *block*.
    col_start : int
        Global column offset of *block*.
    n_total : int
        Full matrix dimension.
    """
    import numpy as _np

    qsize, rsize = block.shape
    block_f16 = _np.asarray(block, dtype=_np.float16)

    rows = _np.arange(qsize, dtype=_np.int64) + row_start
    cols_all = _np.arange(rsize, dtype=_np.int64) + col_start
    upper_mask = (rows[:, None] < cols_all[None, :]) & (cols_all[None, :] < n_total)
    if not _np.any(upper_mask):
        return
    ri = _np.broadcast_to(rows[:, None], (qsize, rsize))[upper_mask]
    cj = _np.broadcast_to(cols_all[None, :], (qsize, rsize))[upper_mask]
    global_rows = ri.astype(_np.int64, copy=False)
    flat_idx = (
        global_rows * _np.int64(n_total)
        - global_rows * (global_rows + 1) // 2
        + cj
        - global_rows
        - 1
    )
    tri_data[flat_idx] = block_f16[upper_mask]
