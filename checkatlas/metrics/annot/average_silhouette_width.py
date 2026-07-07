from typing import Callable, Optional

import numpy as np
from scipy.sparse import issparse
from sklearn.metrics import pairwise_distances

# ── Centroid-variance approx (shared with cluster silhouette) ──────
from ..cluster.silhouette import _centroid_variance_silhouette

# ── GPU cdist import (optional, falls back to sklearn) ────────────
_gpu_cdist: Optional[Callable[..., np.ndarray]] = None
_HAS_JAX: bool = False
try:
    from .._jax_utils import _JAX_AVAILABLE as _HAS_JAX
    from .._jax_utils import cdist as _gpu_cdist
except ImportError:
    pass

# ── GPU silhouette imports (optional) ─────────────────────────────
_HAS_SILHOUETTE_GPU: bool = False
try:
    import jax  # noqa: F401
    import jax.numpy as jnp  # noqa: F401

    from .._jax_utils import _GPU_AVAILABLE, _gpu_free_memory_bytes  # noqa: F811

    _HAS_SILHOUETTE_GPU = _HAS_JAX and _GPU_AVAILABLE
except ImportError:
    pass

_TRIANGULAR_DENSE_THRESHOLD = 30_000
"""At or below this number of cells, :func:`TriangularMatrix.to_dense` is
called to produce a full float32 matrix for sklearn's vectorised
``silhouette_score(metric="precomputed")``.  Above the threshold the
chunked :func:`_silhouette_from_triangular` path is taken to avoid
materialising the full N×N dense array."""


def run(X, labels, n_jobs=-1, verbose=True, sample_size=None, precomputed_dists=None, method="approx"):
    """
    Calculate the Average Silhouette Width for clustering quality.

    Silhouette measures how well each point fits within its assigned
    cluster relative to its nearest neighbouring cluster.
    Values range from -1 (misassigned) to +1 (well assigned).

    Three computation paths, resolved in priority order:

    1. **Centroid-Variance Approximation** (default, ``method="approx"``) --
       O(N·K·d) using the bluster formula:
       ``D_tilde(i, C_X) = sqrt(||x_i - mu_X||**2 + sigma_X**2)``
       Takes priority over precomputed distances.  Runs in seconds.

    2. **Precomputed distances** (``precomputed_dists``, ``method="exact"``) --
       exact result from a precomputed N×N matrix or TriangularMatrix.

    3. **Exact raw** (``method="exact"``, no precomputed) -- full intra-cluster
       and inter-cluster mean Euclidean distances via chunked blocks.  O(N^2·d).

    `Average Silhouette Width readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/clustering/silhouette/>`__

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix.  Sparse inputs are densified.
    labels : array-like of shape (n_samples,)
        Cluster labels.
    n_jobs : int, default=-1
        Parallel jobs (-1 = all cores).  Exact path only.
    verbose : bool, default=True
        Print progress.
    sample_size : int, optional
        **Ignored** -- always processes every sample.
    precomputed_dists : ndarray or TriangularMatrix, optional
        Precomputed N×N pairwise distance matrix.  Only used when
        ``method="exact"``.
    method : str, optional (default="approx")
        - ``"approx"`` -- Centroid-Variance Approximation (bluster method), fast
        - ``"exact"`` -- exact pairwise distances (precomputed or raw)

    Returns
    -------
    float
        Mean silhouette coefficient.  Range [-1, 1].
    """
    # Path 1: Centroid-Variance Approximation (default, O(N·K·d))
    if method == "approx":
        if issparse(X):
            X = X.toarray()
        labels_arr = np.asarray(labels)
        if len(np.unique(labels_arr)) < 2:
            return 0.0
        return _centroid_variance_silhouette(X, labels_arr, verbose=verbose)

    # Path 2: Exact — precomputed distances
    if precomputed_dists is not None:
        if issparse(precomputed_dists):
            precomputed_dists = precomputed_dists.toarray()
        labels_arr = np.asarray(labels)
        n_total = precomputed_dists.shape[0]
        if n_total != len(labels_arr):
            raise ValueError(
                f"precomputed_dists and labels must have same length. "
                f"Got dists: {n_total}, labels: {len(labels_arr)}"
            )
        if len(np.unique(labels_arr)) < 2:
            return 0.0
        if hasattr(precomputed_dists, "to_dense"):
            if precomputed_dists.n <= _TRIANGULAR_DENSE_THRESHOLD:
                if verbose:
                    print("TriangularMatrix — densifying for Silhouette (sklearn)...")
                dense = precomputed_dists.to_dense()
                from sklearn.metrics import silhouette_score

                return float(silhouette_score(dense, labels_arr, metric="precomputed"))
            else:
                if verbose:
                    print("TriangularMatrix — computing Silhouette from precomputed distances...")
                return _silhouette_from_triangular(
                    precomputed_dists, labels_arr, verbose=verbose
                )
        else:
            if verbose:
                print("Using precomputed distances for Silhouette...")
            from sklearn.metrics import silhouette_score

            return float(silhouette_score(precomputed_dists, labels_arr, metric="precomputed"))

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    labels = np.asarray(labels)
    n_total = len(X)

    if n_total != len(labels):
        raise ValueError(
            f"X and labels must have same length. "
            f"Got X: {n_total}, labels: {len(labels)}"
        )

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    if n_clusters < 2:
        return 0.0

    if verbose:
        print(
            f"Computing Average Silhouette Width "
            f"({n_total:,} samples, {n_clusters} clusters, n_jobs={n_jobs})..."
        )

    # ── Step 1: a(i) — exact intra-cluster mean Euclidean distance ─
    a_values = np.zeros(n_total)
    cluster_sizes = {}

    for lbl in unique_labels:
        idx = np.where(labels == lbl)[0]
        nk = len(idx)
        cluster_sizes[lbl] = nk
        if nk < 2:
            continue
        Xk = X[idx]
        if verbose and nk > 5000:
            print(f"  Intra-cluster '{lbl}' ({nk:,} samples)…")

        _exact_intra_mean(Xk, a_values, idx, n_jobs=n_jobs, chunk=2000)

    # ── Step 2: b(i) via exact chunked mean Euclidean distance ────
    # For each pair (Ci, Cj), compute the mean Euclidean distance from
    # each point in Ci to all points in Cj.  This is the *exact* mean
    # distance required by the silhouette formula, not the centroid
    # RMS approximation which gives wrong values.
    #
    # Chunking over Cj keeps peak memory at O(|Ci| · chunk_j) instead
    # of O(|Ci| · |Cj|).  Total work is O(N² · d) across all pairs.
    b_values = np.full(n_total, np.inf)

    for lbl_i in unique_labels:
        idx_i = np.where(labels == lbl_i)[0]
        Xi = X[idx_i]
        ni = len(idx_i)

        b_i = np.full(ni, np.inf)
        norms_i = (Xi ** 2).sum(axis=1)  # (ni,)

        for lbl_j in unique_labels:
            if lbl_i == lbl_j:
                continue
            idx_j = np.where(labels == lbl_j)[0]
            Xj = X[idx_j]
            nj = len(idx_j)

            # Compute mean Euclidean distance from each point in Xi to
            # all points in Xj.  Chunk over Xj to keep memory bounded.
            sum_dist = np.zeros(ni)
            norms_j = (Xj ** 2).sum(axis=1)
            _chunk_j = 4000
            for cs in range(0, nj, _chunk_j):
                ce = min(cs + _chunk_j, nj)
                Xjc = Xj[cs:ce]
                norms_jc = norms_j[cs:ce]
                # Expansion: ||Xi - Xjc||² = ||Xi||² + ||Xjc||² - 2·Xi·Xjc.T
                dot = Xi @ Xjc.T
                dists_sq = norms_i[:, None] + norms_jc[None, :] - 2.0 * dot
                dists_sq = np.maximum(dists_sq, 0.0)
                sum_dist += np.sqrt(dists_sq).sum(axis=1)

            mean_dist = sum_dist / nj
            b_i = np.minimum(b_i, mean_dist)

        b_values[idx_i] = b_i

    if verbose:
        print("  Completed b(i) via exact chunked mean distance.")

    # ── Step 3: per-point silhouette ──────────────────────────────
    s = _silhouette_from_apoint(a_values, b_values)

    if verbose:
        print(f"  ASW = {s:.6f}")

    return s


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════


# ── GPU silhouette (single‑shot + chunked) ─────────────────────────
_GPU_SILHOUETTE_SAFETY = 0.80  # fraction of free GPU memory to use (80 % cap)
_GPU_SILHOUETTE_MAX_N = 100_000  # hard cap, even on large GPUs
_GPU_JAX_OVERHEAD = 3.0  # JAX matmul working memory + autotune multiplier
_GPU_MIN_CHUNK = 1_000  # minimum rows per GPU chunk


def _compute_chunk_size(free_gpu_bytes: int, n: int) -> int:
    """Return row chunk size so the working set fits on GPU.

    The block on GPU is ``chunk_size × n`` float32 elements.
    With JAX matmul overhead the working set is roughly
    ``chunk_size × n × 4 × (1 + JAX_OVERHEAD)`` bytes.
    """
    if free_gpu_bytes <= 0:
        return 0
    bytes_per_element = 4  # float32
    budget = free_gpu_bytes * _GPU_SILHOUETTE_SAFETY
    chunk = int(budget / (n * bytes_per_element * (1 + _GPU_JAX_OVERHEAD)))
    return max(_GPU_MIN_CHUNK, min(n, chunk))


def _silhouette_gpu_chunked(tri, labels, n: int, verbose: bool = True):
    """GPU‑chunked silhouette — batched matmul per row chunk.

    When the full N×N matrix does not fit on the GPU this function
    reads one row‑chunk at a time (``chunk_size × N`` float32),
    uploads it to the GPU, and computes **all** per‑cluster means
    for that chunk via a single JAX matmul
    (``block @ label_mask.T / cluster_sizes``).

    Returns the mean silhouette coefficient **or None** when the
    GPU path cannot be used.
    """
    if not _HAS_SILHOUETTE_GPU:
        return None

    labels_np = np.asarray(labels)
    unique_labels = np.unique(labels_np)
    n_clusters = len(unique_labels)
    if n_clusters < 2:
        return 0.0

    # ── Memory gate ────────────────────────────────────────────────
    try:
        free = _gpu_free_memory_bytes()
    except Exception:
        free = 0
    if free <= 0:
        if verbose:
            print("GPU memory query failed — falling back to CPU silhouette")
        return None

    chunk_size = _compute_chunk_size(free, n)
    if chunk_size >= n:
        # Chunk size equals N — use single-shot instead (one matmul)
        return None  # let caller try _silhouette_gpu
    if chunk_size < _GPU_MIN_CHUNK:
        if verbose:
            print(
                f"GPU chunk size too small ({chunk_size:,} rows) — "
                f"falling back to CPU"
            )
        return None

    if verbose:
        n_chunks = (n + chunk_size - 1) // chunk_size
        print(
            f"GPU chunked silhouette: "
            f"{n_chunks} chunks × {chunk_size:,} rows, "
            f"{n_clusters} clusters"
        )

    # ── Build label mask on CPU (small: n_clusters × n) ────────────
    label_indices = np.searchsorted(unique_labels, labels_np).astype(np.int32)
    label_mask = np.zeros((n_clusters, n), dtype=np.float32)
    for k in range(n_clusters):
        label_mask[k, label_indices == k] = 1.0
    cluster_sizes_arr = label_mask.sum(axis=1)

    # ── Prepare accumulators ───────────────────────────────────────
    a_values = np.zeros(n, dtype=np.float64)
    b_values = np.full(n, np.inf, dtype=np.float64)

    import jax.numpy as jnp
    from tqdm import tqdm

    n_chunks = (n + chunk_size - 1) // chunk_size
    pbar = tqdm(
        total=n_chunks,
        desc="  Silhouette (GPU chunked)",
        disable=not verbose,
    )

    # Warm-up JIT
    _warm = jnp.mean(jnp.zeros((min(100, chunk_size), 2), dtype=jnp.float32), axis=1)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        cur_chunk = end - start

        # ── Read all columns for this chunk ─────────────────────────
        block = tri[start:end, :]  # (cur_chunk, n) float32

        # ── Upload to GPU ───────────────────────────────────────────
        block_gpu = jnp.asarray(block)
        label_mask_gpu = jnp.asarray(label_mask)
        cluster_sizes_gpu = jnp.asarray(cluster_sizes_arr)

        # ── Compute all cluster means at once ───────────────────────
        # (cur_chunk, n) @ (n, n_clusters) / cluster_sizes
        #   → (cur_chunk, n_clusters)
        cluster_means = (
            block_gpu @ label_mask_gpu.T
        ) / cluster_sizes_gpu[None, :]

        # ── Download to CPU ─────────────────────────────────────────
        cluster_means_np = np.asarray(cluster_means)

        # ── a(i): rows that belong to each cluster ──────────────────
        chunk_labels = labels_np[start:end]
        chunk_label_idx = label_indices[start:end]
        for k in range(n_clusters):
            nk = int(cluster_sizes_arr[k])
            if nk < 2:
                continue
            mask_in = chunk_labels == unique_labels[k]
            if np.any(mask_in):
                a_global_idx = np.arange(start, end)[mask_in]
                a_values[a_global_idx] = (
                    cluster_means_np[mask_in, k] * nk / (nk - 1)
                )

        # ── b(i): min over all other clusters ───────────────────────
        cmeans = cluster_means_np.copy()
        cmeans[np.arange(cur_chunk), chunk_label_idx] = np.inf
        b_chunk = np.min(cmeans, axis=1)
        b_values[start:end] = np.minimum(b_values[start:end], b_chunk)

        pbar.update(1)

    pbar.close()

    denom = np.maximum(a_values, b_values)
    mask = denom > 0
    s = np.zeros(n)
    s[mask] = (b_values[mask] - a_values[mask]) / denom[mask]
    return float(np.mean(s))


def _silhouette_gpu(tri, labels, n: int, verbose: bool = True):
    """GPU‑accelerated silhouette from a TriangularMatrix (JAX single‑shot).

    Returns the mean silhouette coefficient **or None** if the GPU
    path cannot be used (insufficient memory, JAX unavailable, etc.).
    The caller must fall back to the CPU path when ``None`` is returned.

    Algorithm
    ---------
    1. Densify the TriangularMatrix to a host‑side float16 matrix via
       :meth:`~TriangularMatrix._to_dense_f16`.
    2. Upload to GPU and compute per‑cluster mean distances via a single
       JAX matrix‑multiply (``dists @ label_mask.T``).
    3. Compute per‑point *a(i)* and *b(i)* entirely on the GPU and
       transfer only the final scalar back to the host.

    Memory budget (conservative)
    ----------------------------
    * float16 dense matrix on GPU:      ``N² × 2`` bytes
    * JAX float32 matmul working set:   ``N² × 4`` bytes
    * Misc (label_mask, cluster_means): ``N × K × 8`` bytes  (negligible)
    * Total ≈ ``N² × 6`` bytes

    With the :data:`_GPU_SILHOUETTE_SAFETY` factor (45 % of free GPU
    memory) this limits the single‑shot path on common hardware:

    ===========  ===========
    GPU          max N ≈
    ===========  ===========
    A100‑40 GB   ~54 000
    A6000‑48 GB  ~60 000
    A100‑80 GB   ~77 000
    ===========  ===========
    """
    if not _HAS_SILHOUETTE_GPU:
        return None

    n_clusters = len(np.unique(labels))
    if n_clusters < 2:
        return 0.0

    # ── Memory gate ────────────────────────────────────────────────
    try:
        free = _gpu_free_memory_bytes()
    except Exception:
        free = 0
    if free <= 0:
        if verbose:
            print("GPU memory query failed — falling back to CPU silhouette")
        return None

    required = n * n * 8  # f16 dense + f32 JAX matmul + autotune
    if required > free * _GPU_SILHOUETTE_SAFETY or n > _GPU_SILHOUETTE_MAX_N:
        if verbose:
            print(
                f"GPU memory insufficient for single-shot silhouette "
                f"(N={n:,}, req={required/1024**3:.1f} GB, "
                f"free={free/1024**3:.1f} GB) — falling back to CPU"
            )
        return None

    # ── Densify (host side, float16) ───────────────────────────────
    if verbose:
        import time as _time

        _t0 = _time.time()
        print("Densifying TriangularMatrix for GPU silhouette (float16)...")
    dense_f16 = tri._to_dense_f16()
    if verbose:
        print(f"  Densified in {_time.time() - _t0:.1f}s, uploading to GPU...")

    # ── GPU computation ────────────────────────────────────────────
    import jax.numpy as jnp

    _t1 = _time.time() if verbose else None

    dists = jnp.asarray(dense_f16)
    del dense_f16  # free host memory early

    labels_np = np.asarray(labels)
    unique_labels = np.unique(labels_np)
    label_indices = np.searchsorted(unique_labels, labels_np).astype(np.int32)

    # Build label_mask: (n_clusters, n) float32 indicator matrix
    label_idx_jax = jnp.asarray(label_indices)
    label_mask = jax.nn.one_hot(label_idx_jax, n_clusters, dtype=jnp.float32).T
    cluster_sizes = label_mask.sum(axis=1)  # (n_clusters,)

    # dists @ label_mask.T: (n, n) f16 × (n, n_clusters) f32 → (n, n_clusters) f32
    cluster_means = (dists @ label_mask.T) / cluster_sizes[None, :]

    # a(i): mean distance to own cluster (with nk/(nk-1) correction)
    a_values = cluster_means[jnp.arange(n), label_idx_jax]
    nk = cluster_sizes[label_idx_jax]
    a_values = jnp.where(nk > 1, a_values * nk / (nk - 1), 0.0)

    # b(i): min mean distance to any OTHER cluster
    cluster_means = cluster_means.at[jnp.arange(n), label_idx_jax].set(jnp.inf)
    b_values = jnp.min(cluster_means, axis=1)

    denom = jnp.maximum(a_values, b_values)
    s = jnp.where(denom > 0, (b_values - a_values) / denom, 0.0)
    result = float(jnp.mean(s))

    if verbose:
        print(
            f"  GPU silhouette computed in {_time.time() - _t1:.1f}s"
            f"  →  {result:.6f}"
        )

    return result


def _silhouette_from_triangular(tri, labels, verbose=True):
    """Compute silhouette score from a TriangularMatrix distance matrix.

    Processes one cluster at a time for all rows, reading only the
    columns belonging to that cluster from the TriangularMatrix via
    column‑subsetted ``_get_block``.  This *avoids* the intermediate
    read‑all‑columns‑then‑subset double‑copy pattern of the previous
    implementation (`tri[start:end, :]` → `block[:, idx]`) by
    reading the distance sub‑matrix ``tri[start:end, idx]`` directly.
    The per‑cluster boolean masks in ``_get_block`` are therefore
    ``chunk_size × |cluster|`` instead of ``chunk_size × N`` — a
    **50× reduction** for the typical ~1 700‑cell cluster — which
    cuts the total temporary allocation from ~263 GB to ~60 GB per
    call at N=85 000 and eliminates the second pass through the data.

    Total work is still O(N²) element reads (inherent to silhouette)
    but the constant‑factor overhead is roughly halved.
    """
    n = tri.shape[0]
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    # ── Try GPU single-shot path ───────────────────────────────────
    if _HAS_SILHOUETTE_GPU:
        gpu_result = _silhouette_gpu(tri, labels, n=n, verbose=verbose)
        if gpu_result is not None:
            return gpu_result

        # ── Try GPU chunked path ─────────────────────────────────────
        gpu_chunked_result = _silhouette_gpu_chunked(
            tri, labels, n=n, verbose=verbose
        )
        if gpu_chunked_result is not None:
            return gpu_chunked_result

    cluster_idx = {lbl: np.where(labels == lbl)[0] for lbl in unique_labels}
    cluster_sizes = {lbl: len(idx) for lbl, idx in cluster_idx.items()}

    a_values = np.zeros(n, dtype=np.float64)
    b_values = np.full(n, np.inf, dtype=np.float64)

    _ROW_CHUNK = max(500, min(5000, 50_000_000 // max(n, 1)))

    from tqdm import tqdm

    pbar = tqdm(
        total=n_clusters,
        desc="  Silhouette (precomputed)",
        disable=not verbose,
    )

    for lbl, idx in cluster_idx.items():
        nk = cluster_sizes[lbl]

        for start in range(0, n, _ROW_CHUNK):
            end = min(start + _ROW_CHUNK, n)

            # Read D[start:end, cluster_cols] directly — no intermediate
            # full‑column block, no subsequent fancy‑index copy.
            block = tri[start:end, idx]
            col_means = block.mean(axis=1, dtype=np.float64)

            row_labels = labels[start:end]

            # a(i): rows that belong to this cluster
            mask_in = row_labels == lbl
            if np.any(mask_in) and nk > 1:
                a_global_idx = np.arange(start, end)[mask_in]
                a_values[a_global_idx] = col_means[mask_in] * nk / (nk - 1)

            # b(i): rows that belong to OTHER clusters — keep minimum
            mask_out = row_labels != lbl
            if np.any(mask_out):
                b_global_idx = np.arange(start, end)[mask_out]
                b_values[b_global_idx] = np.minimum(
                    b_values[b_global_idx], col_means[mask_out]
                )

        pbar.update(1)

    pbar.close()

    denom = np.maximum(a_values, b_values)
    mask = denom > 0
    s = np.zeros(n)
    s[mask] = (b_values[mask] - a_values[mask]) / denom[mask]
    return float(np.mean(s))


def _exact_intra_mean(Xk, a_out, global_idx, n_jobs=1, chunk=2000):
    """
    Compute a(i) = mean Euclidean distance from each point to all
    *other* points in the same cluster.  Uses chunked upper-triangular
    pairwise blocks so the full |Ck|×|Ck| matrix is never stored.

    For each off-diagonal block (r < c), the column sums are also
    accumulated into a_out so that points in Xc (columns) get their
    cross-distances counted too — without this, the last chunk of
    each cluster row misses half its distances.

    When JAX is available, GPU cdist replaces sklearn pairwise_distances
    for ~100× faster intra-cluster distance computation.
    """
    n = len(Xk)
    counts = np.zeros(n)  # per-point count of other points seen
    _use_gpu = _HAS_JAX and _gpu_cdist is not None

    for r in range(0, n, chunk):
        re = min(r + chunk, n)
        Xr = Xk[r:re]
        for c in range(r, n, chunk):
            ce = min(c + chunk, n)
            Xc = Xk[c:ce]
            if _use_gpu:
                block = np.array(_gpu_cdist(Xr, Xc), copy=True)  # JAX returns read-only
            else:
                block = pairwise_distances(Xr, Xc, n_jobs=n_jobs)
            # block shape: (re-r) × (ce-c)

            if r == c:
                # Diagonal block — zero out self-distances, sum rows.
                dlen = block.shape[0]
                for d in range(dlen):
                    if d < block.shape[1]:
                        block[d, d] = np.nan
                with np.errstate(all="ignore"):
                    row_sum = np.nansum(block, axis=1)
                row_cnt = (ce - c) - 1  # exclude self
                a_out[global_idx[r:re]] += row_sum
                counts[r:re] += row_cnt
            else:
                # Off-diagonal block — accumulate BOTH row sums (for Xr
                # points) and column sums (for Xc points).  Each
                # unordered pair is counted exactly once from each
                # direction, giving each point its full n−1 count.
                row_sum = block.sum(axis=1)
                col_sum = block.sum(axis=0)
                row_cnt = ce - c
                col_cnt = re - r
                a_out[global_idx[r:re]] += row_sum
                a_out[global_idx[c:ce]] += col_sum
                counts[r:re] += row_cnt
                counts[c:ce] += col_cnt

    # Normalise: a(i) = sum / (n-1)
    mask = counts > 0
    a_out[global_idx[mask]] /= counts[mask]


def _silhouette_from_apoint(a, b):
    """Per-point silhouette from precomputed a(i) and b(i)."""
    denom = np.maximum(a, b)
    mask = denom > 0
    s = np.zeros(len(a))
    s[mask] = (b[mask] - a[mask]) / denom[mask]
    return float(np.mean(s))


# Keep old name for backward-compat references
_silhouette_per_point = _silhouette_from_apoint
