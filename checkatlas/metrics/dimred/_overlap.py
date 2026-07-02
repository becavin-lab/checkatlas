import numpy as np


def _overlap_count(low_idx: np.ndarray, high_idx: np.ndarray) -> np.ndarray:
    """Count the intersection size between each pair of k-NN index rows.

    For two ``(N, k)`` integer index matrices, returns a length-``N``
    array of the number of indices that appear in *both* the low-dim and
    high-dim kNN lists for each cell.  Set semantics: a value that
    appears twice in a single row counts once in the intersection, which
    matches the prior ``set(low) & set(high)`` semantic and is robust to
    kNN indices that contain ties.

    Implementation
    --------------
    Sort each row's concatenated 2k indices and count unique values
    in the union via the standard "diff+1" trick on a sorted array.
    Then ``|intersection| = |low_set| + |high_set| - |union|`` and
    both ``|low_set|`` and ``|high_set|`` are obtained by the same
    sorted-diff trick applied to the per-row k arrays.

    For k=30 the temporary sort is over an (N, 2k=60) int64 array
    (~5 MB for N=85 000) and runs in O(N k log k) — well under a
    second on any modern CPU.  This replaces the previous
    ``joblib``+``set`` per-row implementation, which cost ~5-15 µs
    of pickling + Python overhead per cell.

    Parameters
    ----------
    low_idx, high_idx : np.ndarray
        Integer index matrices of identical shape ``(N, k)``.

    Returns
    -------
    np.ndarray of shape ``(N,)``, dtype ``int64``.
    """
    if low_idx.shape != high_idx.shape:
        raise ValueError(
            f"shape mismatch: low_idx {low_idx.shape} vs high_idx {high_idx.shape}"
        )
    if low_idx.ndim != 2:
        raise ValueError(f"expected 2-D index matrices; got ndim={low_idx.ndim}")

    n, k = low_idx.shape

    low_set_size = _row_unique_count(low_idx)
    high_set_size = _row_unique_count(high_idx)
    union = _row_unique_count(np.concatenate([low_idx, high_idx], axis=1))
    intersection = low_set_size + high_set_size - union
    return intersection.astype(np.int64)


def _row_unique_count(idx: np.ndarray) -> np.ndarray:
    """Count unique values per row of a 2-D array.

    Uses the standard sorted-difference trick: sort each row, then
    count ``np.sum(diff != 0) + 1``.  O(N k log k).
    """
    sorted_idx = np.sort(idx, axis=1, kind="stable")
    diffs = np.diff(sorted_idx, axis=1)
    return (diffs != 0).sum(axis=1) + 1
