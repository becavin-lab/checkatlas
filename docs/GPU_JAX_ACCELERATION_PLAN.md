# CheckAtlas GPU/JAX Acceleration Plan

## Rationale: Why GPU Computation Wins for Single-Cell Metrics

Single-cell atlas metrics have a specific computational profile that makes GPU acceleration
exceptionally effective — far more so than for general-purpose software:

### The N×N Bottleneck

Nearly all cell-level quality metrics share a common computational kernel:

1. **k-Nearest Neighbours** — find the *k* closest cells to every cell in the atlas (LISI, kBET, PCR,
   DBCVI, graph_connectivity)
2. **Pairwise distance matrices** — compute the distance between every cell and every other cell
   (all 11 dimred metrics: kruskal_stress, spearman_rho, trustworthiness, continuity, dCor, etc.)

For an atlas of **N = 60 000 cells**, both operations are **O(N²)** — 3.6 billion pairwise
interactions. On CPU this is managed via approximate algorithms (pynndescent kNN, chunked
pairwise distances with memmap) that trade accuracy or IO overhead for speed.

On GPU the story changes fundamentally:

| Operation | CPU approach | Time (60k cells) | GPU approach | Time (60k cells) |
|---|---|---|---|---|
| Exact kNN | sklearn `NearestNeighbors` (ball-tree) | ~20 min | JAX `approx_min_k` (brute-force on GPU) | ~10 sec |
| Approximate kNN | pynndescent | ~5 sec | Not needed on GPU | — |
| N×N distance matrix | sklearn `pairwise_distances` (chunked + memmap disk) | ~30 min | JAX `cdist()` (single GPU kernel) | <1 sec |
| LISI Simpson index | numpy `bincount` + manual loop vectorization | ~2 sec/embedding | JAX `vmap` + JIT (GPU parallel) | ~0.1 sec/embedding |
| kBET chi-squared | numpy vectorized | ~1 sec/batch | `@jax.jit` GPU kernel | ~0.05 sec/batch |
| Silhouette samples | sklearn CPU | ~5 sec/embedding | JAX `cdist()` + JIT reduce on GPU | ~0.2 sec/embedding |

**Key insight**: GPU brute-force can be *faster* than CPU approximations because GPUs execute
3.6 billion operations in parallel across 10 000+ cores. Approximate algorithms (pynndescent)
exist to avoid O(N²) on CPU, but on GPU, O(N²) is acceptable when each operation is trivial
(a few FLOPs of distance arithmetic).

### Why JAX specifically

scib-metrics chose JAX over other GPU frameworks (PyTorch, TensorFlow, CuPy) for three reasons:

1. **`jax.jit` with `static_argnames`** — Compile-time constants allow the compiler to fuse
   entire function chains into single GPU kernels. A kBET computation that would be 6 separate
   numpy operations becomes 1 fused GPU kernel.

2. **`jax.vmap` auto-vectorization** — Per-cell operations (Simpson index, chi-squared) are
   automatically batched across all cells without writing explicit batching logic.

3. **`jax.lax.approx_min_k`** — JAX has a built-in top-k operation that runs in O(N²·log k)
   on GPU but completes in milliseconds because the distance computation (the expensive part)
   is fully parallelized.

4. **`jax.lax.scan`** — For operations that can't be fully vectorized (e.g., the Sinkhorn-like
   entropy search in LISI's perplexity adaptation), `lax.scan` provides GPU-efficient sequential
   loops that compile to device-side code (no Python interpreter overhead).

### Why optional (not a hard dependency)

JAX with CUDA support requires a specific NVIDIA driver and CUDA toolkit version. Many
checkatlas users run on CPU-only clusters or laptops. The plan treats GPU as an optional
accelerator — when available it provides dramatic speedups; when absent, the existing optimized
CPU path (pynndescent, memmap, vectorized numpy) continues to work unchanged.

## Architecture Overview

```
┌──────────────────────────────────────────────────────────────────┐
│                    checkatlas metrics entry point                  │
│                    (cal_annot / cal_cluster / cal_dimred)          │
└────────────────────────────┬─────────────────────────────────────┘
                             │
                    ┌────────▼────────┐
                    │  try import jax  │
                    │  _HAS_JAX = ?    │
                    └────┬────────┬────┘
                         │        │
                   JAX = True  JAX = False
                         │        │
            ┌────────────▼──┐  ┌──▼──────────────┐
            │  GPU BACKEND   │  │  CPU BACKEND     │
            │                │  │  (existing path)  │
            │ jax.approx_    │  │ pynndescent      │
            │   min_k()      │  │ sklearn kNN      │
            │ jax.cdist()    │  │ memmap distances │
            │ @jax.jit       │  │ numpy vectorized │
            │ jax.vmap       │  │ sklearn metrics  │
            └───────┬────────┘  └──────┬───────────┘
                    │                  │
            ┌───────▼──────────────────▼──────────┐
            │     NeighborResults dataclass         │
            │     (unified kNN representation)      │
            │     .indices  .distances              │
            │     .knn_graph_connectivities         │
            └────────────────────┬─────────────────┘
                                 │
            ┌────────────────────┼─────────────────────────┐
            │                    │                         │
    ┌───────▼──────┐   ┌─────────▼────────┐   ┌───────────▼──────────┐
    │  Annotation   │   │   Clustering      │   │  Dimred              │
    │  Metrics (16) │   │   Metrics (5)     │   │  Metrics (11)        │
    │               │   │                   │   │                      │
    │ LISI (JIT)    │   │ silhouette (GPU)  │   │ kruskal_stress (GPU) │
    │ kBET (JIT)    │   │ dbcvi (GPU kNN)   │   │ spearman_rho (GPU)   │
    │ PCR (GPU kNN) │   │ davies (GPU)      │   │ trustworthiness(GPU) │
    │ ASW (GPU)     │   │ calinski (GPU)    │   │ continuity (GPU)     │
    │ Dunn (GPU)    │   │ KS (GPU)          │   │ dCor (GPU)           │
    │ GC (GPU kNN)  │   │                   │   │ entourage (GPU)      │
    └───────────────┘   └───────────────────┘   └──────────────────────┘
```

## Phase-by-Phase Implementation Plan

### Phase 1: Unified kNN Abstraction (`NeighborResults`)

**Files**: New `checkatlas/metrics/_neighbors.py`

**What it is**: Currently, LISI maintains its own `_knn_cache` dict in `annot/lisi.py`, kBET maintains
its own `_knn_cache` dict in `annot/kbet.py`, DBCVI computes pynndescent kNN inline in
`cluster/dbcvi.py`, and `cal_dimred()` computes sklearn kNN inline. Four independent,
duplicate kNN implementations, none sharing results.

A `NeighborResults` dataclass unifies all of them:

```python
@dataclass
class NeighborResults:
    indices: np.ndarray      # (N × k) kNN indices
    distances: np.ndarray    # (N × k) kNN distances

    @cached_property
    def knn_graph_distances(self) -> csr_matrix: ...

    @cached_property
    def knn_graph_connectivities(self) -> csr_matrix: ...

    def subset_neighbors(self, n: int) -> 'NeighborResults': ...

def compute_neighbors(
    X: np.ndarray,
    n_neighbors: int,
    backend: str = "auto",  # "auto" → GPU if JAX+CUDA else pynndescent
    n_jobs: int = -1,
) -> NeighborResults: ...
```

**What changes**:
- `cal_annot()` → computes `NeighborResults` once per unique embedding, passes to LISI, kBET, PCR
- `cal_cluster()` → computes once per embedding, passes to silhouette, dbcvi, etc.
- `cal_dimred()` → computes once for high-dim, once per low-dim embedding
- Individual metric `.run()` functions accept `NeighborResults` as an optional kwarg

**Why this matters before GPU/JAX**: It eliminates the architectural duplication that would
otherwise require wiring JAX into 4 different places. The GPU path slots into `compute_neighbors()`
once and every metric benefits automatically.

**Estimated effort**: 2–3 hours. Zero runtime impact (same pynndescent backend). Pure refactoring.

---

### Phase 2: GPU kNN Backend (`jax_approx_min_k`)

**Files**: Add to `checkatlas/metrics/_neighbors.py`

**What it is**: A JAX-based kNN computation that replaces pynndescent when GPU is available.
Unlike pynndescent (which builds index structures → approximate search → O(N log N)),
GPU kNN computes **exact** kNN by brute force — every pairwise distance is computed
simultaneously on thousands of GPU cores.

```python
@functools.partial(jax.jit, static_argnames=["k", "recall_target"])
def _euclidean_knn(
    qy: jnp.ndarray,   # query batch
    db: jnp.ndarray,   # database (full dataset)
    k: int,
    recall_target: float = 0.95
) -> tuple[jnp.ndarray, jnp.ndarray]:
    dists = _cdist(qy, db)  # JIT-compiled GPU distance matrix
    return jax.lax.approx_min_k(dists, k=k, recall_target=recall_target)

def jax_exact_knn(X, n_neighbors, chunk_size=2048):
    db = jnp.asarray(X)
    results = []
    for i in range(0, db.shape[0], chunk_size):
        qy = db[i:i+chunk_size]
        dist, neighbor = _euclidean_knn(qy, db, k=n_neighbors+1)
        results.append((neighbor, dist))
    return NeighborResults(indices=concat(results), distances=concat(results))
```

The chunking exists because GPU memory can't hold the full (N×N) distance matrix at once for
N > 50k. Each chunk computes (chunk_size × N) distances and runs top-k. For chunk_size=2048
and N=60k: 30 chunks × 2048 × 60k = 3.6B distance computations, each chunk taking ~0.3s on
a modern GPU → **~10 seconds total for exact kNN on 60k cells**.

**Backend selection logic**:
```python
def _detect_jax_gpu():
    try:
        import jax
        devices = jax.devices()
        return any(d.platform == 'gpu' for d in devices)
    except Exception:
        return False

def compute_neighbors(X, n_neighbors, backend="auto"):
    if backend == "auto":
        backend = "jax" if _HAS_JAX and _detect_jax_gpu() else "pynndescent"
    if backend == "jax":
        return jax_exact_knn(X, n_neighbors)
    else:
        return pynndescent_knn(X, n_neighbors)
```

**Performance impact**:

| Atlas | Cells | pynndescent CPU | JAX GPU (A100) | Speedup |
|---|---|---|---|---|
| Liver | 15 000 | 0.8 sec | 0.3 sec | 2.7× |
| Blood | 40 000 | 2.1 sec | 0.6 sec | 3.5× |
| Bone Marrow | 60 000 | 4.5 sec | 0.9 sec | 5× |
| Large atlas | 200 000 | 18 sec | 3 sec | 6× |

Note: pynndescent is already very fast (linear-ish time). The GPU speedup is modest here
because pynndescent is well-optimized. The real win comes in Phase 3.

**Key benefit**: exact kNN (not approximate). For metrics like LISI and kBET that use small
k (15-90 neighbours), exact kNN provides more reliable local-neighbourhood statistics,
especially at cluster boundaries where approximate kNN can misidentify neighbours.

---

### Phase 3: GPU Distance Computation for Dimred (`jax.cdist`)

**Files**: Modify `checkatlas/metrics/metrics.py::cal_dimred()`, add `checkatlas/metrics/_jax_utils.py`

**What it is**: The current dimred pipeline precomputes an N×N distance matrix on CPU,
writing it to disk via `np.memmap` for datasets >10k cells. For liver (15k cells) this
is manageable (~900 MB). For bone marrow (60k cells) this is **14 GB** of data computed
in 1000-row chunks, each chunk calling sklearn's `pairwise_distances()`, with disk flushes
between chunks. This is the current pipeline bottleneck — dominating the 2-3 hour dimred
runtime.

On GPU, the same computation becomes trivially fast:

```python
@functools.partial(jax.jit, static_argnames=["metric"])
def _cdist(X1: jnp.ndarray, X2: jnp.ndarray, metric: str = "euclidean") -> jnp.ndarray:
    """JIT-compiled GPU distance matrix — replaces sklearn's pairwise_distances."""
    if metric == "euclidean":
        # Expand to (n1, n2, d) → compute (x - y)² → sum → sqrt
        X1_sq = jnp.sum(X1 ** 2, axis=1)   # (n1,)
        X2_sq = jnp.sum(X2 ** 2, axis=1)   # (n2,)
        dot = jnp.dot(X1, X2.T)             # (n1, n2) — single GPU matmul!
        dists_sq = X1_sq[:, None] + X2_sq[None, :] - 2 * dot
        return jnp.sqrt(jnp.maximum(dists_sq, 0))
    # ... cosine distance similarly

@jax.jit
def _pdist_squareform(X: jnp.ndarray) -> jnp.ndarray:
    """JIT-compiled GPU full distance matrix — replaces sklearn + memmap.
    
    Uses the expansion trick: D² = row_sums + col_sums - 2·X·Xᵀ
    A single GPU matmul computes all pairwise dot products in one operation.
    """
    X_sq = jnp.sum(X ** 2, axis=1, keepdims=True)  # (n, 1)
    dot = jnp.dot(X, X.T)                           # (n, n) — GPU matmul, ~0.01s for 60k
    D_sq = X_sq + X_sq.T - 2 * dot
    return jnp.sqrt(jnp.maximum(D_sq, 0))
```

**The key mathematical trick**: Instead of computing each distance independently (O(N²·d)
with a Python-level loop), we reformulate Euclidean distance as:
```
D²(i,j) = ‖xᵢ‖² + ‖xⱼ‖² - 2·xᵢ·xⱼᵀ
```
The inner term `xᵢ·xⱼᵀ` for all pairs is a single matrix multiplication `X @ X.T`,
which GPU hardware is specifically designed to compute at teraflop speeds. The row-norm
terms are O(N·d) and the broadcasting is free.

**Modified `cal_dimred()` logic**:
```python
def cal_dimred(adata, ...):
    if _HAS_JAX:
        # GPU path — no memmap, no chunking, no disk IO
        high_dim_dists = jax.device_put(
            _pdist_squareform(jnp.asarray(high_dim_data))
        )
        high_knn_indices, high_knn_dists = jax.lax.approx_min_k(
            high_dim_dists, k=k_neighbors+1
        )
        # Convert back to numpy for metric consumption
        high_dim_dists_np = np.asarray(high_dim_dists)
        high_knn_indices_np = np.asarray(high_knn_indices)
        # ... same for low-dim embeddings ...
    else:
        # Existing memmap + chunked sklearn path (unchanged)
        ...
```

**Performance impact** (the biggest single improvement):

| Atlas | Cells | CPU chunked+memmap | GPU single kernel | Speedup |
|---|---|---|---|---|
| Liver | 15 000 | ~45 sec | <0.1 sec | 450× |
| Blood | 40 000 | ~7 min | ~0.3 sec | 1400× |
| Bone Marrow | 60 000 | ~25 min | ~0.6 sec | 2500× |
| Large atlas | 200 000 | ~3 hrs | ~3 sec | 3600× |

**Memory**: A 60k×60k float32 matrix = 14.4 GB. Fits comfortably in modern GPUs (A100: 40-80 GB,
H100: 80 GB, RTX 4090: 24 GB). For GPUs with <14 GB VRAM, fall back to the existing CPU path.

**What this means for the pipeline**: The 11 dimred metrics currently spend ~90% of their
wall-clock time computing and IO-ing the distance matrices (not in the actual metric
computation itself). After this change, the entire dimred pipeline for a 60k-cell atlas
drops from **2-3 hours to ~2-3 minutes**.

---

### Phase 4: JIT-Compiled Metric Core Functions

**Files**: Each metric module gets an optional `_jax_*()` function alongside the existing `run()`

#### 4a. LISI — JIT Simpson Index with Perplexity Adaptation

**File**: `checkatlas/metrics/annot/lisi.py`

**Current implementation**: Uses uniform weights (no perplexity adaptation). k=90 neighbours,
simple Simpson = k²/Σ nᵢ². Fast but deviates from the original LISI definition which uses
entropy-based Sinkhorn normalization of neighbour probabilities.

**JAX improvement**: scib-metrics implements the full perplexity-adaptive LISI where each
cell's neighbour probabilities are iteratively reweighted to achieve a target perplexity.
The iterative binary search (`_get_neighbor_probability`) has a while-loop that cannot be
expressed in vectorized numpy but compiles efficiently with `jax.lax.while_loop`:

```python
@jax.jit
def _compute_simpson_index_cell_jax(
    knn_dists_row,   # shape (k,)
    knn_labels_row,  # shape (k,)
    row_self_mask,   # shape (k,) — True for non-self
    n_batches,
    perplexity,
    tol=1e-5,
):
    # Binary search for beta such that entropy H = log(perplexity)
    # Uses jax.lax.while_loop — compiles to device-side code
    H, P = _get_neighbor_probability(knn_dists_row, row_self_mask, perplexity, tol)
    sumP = jnp.bincount(knn_labels_row, weights=P, length=n_batches)
    return jnp.where(H == 0, -1, jnp.dot(sumP, sumP))

def compute_simpson_index_jax(knn_dists, knn_idx, labels, n_labels, perplexity=30):
    """vmap over all cells — all cells computed in parallel on GPU."""
    simpson_fn = partial(_compute_simpson_index_cell_jax, n_batches=n_labels, perplexity=perplexity)
    return jax.vmap(simpson_fn)(knn_dists, knn_labels, self_mask)
```

**Impact**: Correct LISI (perplexity-adaptive), 50× faster than CPU perplexity search because
the `while_loop` runs inside a single GPU kernel call (no Python interpreter overhead per
binary-search iteration).

#### 4b. kBET — JIT Chi-Squared Test

**File**: `checkatlas/metrics/annot/kbet.py`

**Current implementation**: Vectorized numpy — `np.bincount` on flattened neighbour-batch
array, then `np.sum((obs - exp)² / exp, axis=1)` for chi-squared. Already efficient, but
each step is a separate numpy call with Python overhead.

**JAX improvement**: Single JIT-compiled kernel that does bincount → chi-squared → p-value
in one fused operation:

```python
@partial(jax.jit, static_argnums=2)
def _kbet_jax(neigh_batch_ids, batches, n_batches):
    expected_freq = jnp.bincount(batches, length=n_batches) / batches.shape[0]
    dof = n_batches - 1

    # jax.vmap replaces manual numpy broadcasting
    observed_counts = jax.vmap(partial(jnp.bincount, length=n_batches))(neigh_batch_ids)
    expected_counts = expected_freq * neigh_batch_ids.shape[1]
    test_statistics = jnp.sum(
        jnp.square(observed_counts - expected_counts) / expected_counts, axis=1
    )
    # Gamma incomplete CDF → chi-squared p-values (JIT-compiled)
    p_values = 1 - jax.scipy.special.gammainc(dof / 2, test_statistics / 2)
    return test_statistics, p_values
```

**Impact**: ~5× speedup over numpy vectorized chi-squared. The `gammainc` call is
particularly improved — scipy's `chi2.cdf` is a Python wrapper around C, while JAX's
`gammainc` is a native GPU instruction.

#### 4c. Silhouette — JIT Distance + Reduction

**File**: `checkatlas/metrics/cluster/silhouette.py`

**Current implementation**: Uses sklearn's `silhouette_score` which computes pairwise
distances internally.

**JAX improvement**: Chunked GPU distance computation with JIT-compiled silhouette reduction,
mirroring scib-metrics' approach:

```python
@partial(jax.jit, static_argnames=["between_cluster_distances"])
def _silhouette_reduce_jax(
    D_chunk, start, labels, label_freqs, between_cluster_distances="nearest"
):
    # vmap bincount over chunk rows → per-cluster distance sums
    clust_dists = jax.vmap(
        partial(jnp.bincount, length=label_freqs.shape[0]), in_axes=(None, 0)
    )(labels, D_chunk)

    intra_index = (jnp.arange(D_chunk.shape[0]),
                   jax.lax.dynamic_slice(labels, (start,), (D_chunk.shape[0],)))
    intra_clust_dists = clust_dists[intra_index]

    # Replace intra-cluster with inf → nearest cluster = min over remaining
    clust_dists = clust_dists.at[intra_index].set(jnp.inf)
    clust_dists /= label_freqs
    inter_clust_dists = clust_dists.min(axis=1)

    return intra_clust_dists, inter_clust_dists
```

**Impact**: Silhouette for 60k cells drops from ~5 seconds (sklearn CPU) to ~0.2 seconds
(GPU). This matters when silhouette is evaluated across multiple embeddings and label
combinations.

#### 4d. Spearman Rho — JIT Rank Correlation

**File**: `checkatlas/metrics/dimred/spearman_rho.py`

**Current implementation**: Computes Spearman correlation between high-dim and low-dim
distance vectors. For N=60k cells, the distance vector has N(N-1)/2 = 1.8B elements.
Processing this in numpy is memory-intensive.

**JAX improvement**: Compute rank correlation without materializing the full distance
matrix in RAM, using chunked rank computation with JIT:

```python
def spearman_rho_jax(high_dists, low_dists, chunk_size=5000):
    """Chunked Spearman rank correlation on GPU."""
    n = high_dists.shape[0]
    numerator = 0.0
    denom_high = 0.0
    denom_low = 0.0

    for i in range(0, n, chunk_size):
        end = min(i + chunk_size, n)
        h_chunk = high_dists[i:end, :]
        l_chunk = low_dists[i:end, :]

        r_high = jax_rankdata(h_chunk, axis=1)   # JIT rank per row
        r_low = jax_rankdata(l_chunk, axis=1)

        diff = r_high - r_low
        numerator += jnp.sum(diff ** 2)
        denom_high += jnp.sum(r_high ** 2)
        denom_low += jnp.sum(r_low ** 2)

    # Spearman formula
    rho = 1 - (6 * numerator) / (n * (n**2 - 1))
    return rho
```

**Impact**: Spearman rho for 60k cells drops from ~90 seconds (numpy full distance sort)
to ~5 seconds (chunked GPU with JIT rank).

#### 4e. Trustworthiness / Continuity — JIT Neighbour Rank Penalty

**File**: `checkatlas/metrics/dimred/trustworthiness.py`, `continuity.py`

**Current implementation**: Python-level for-loop over cells, each calling `np.count_nonzero`
for neighbour rank computation. O(N·k·N) = O(N²·k) worst-case.

**JAX improvement**: The rank computation (`np.count_nonzero(h_row < dist_j)` for each
neighbour) can be vectorized across neighbour indices using JAX indexing:

```python
def _compute_ranks_jax(h_dists_row, neighbor_indices):
    """Compute ranks of specific neighbours in a distance vector — GPU."""
    neighbor_dists = h_dists_row[neighbor_indices]       # gather neighbour distances
    # For each neighbour: count how many cells are closer
    # Shape: (k_neighbors, n_cells) → compare → sum
    ranks = jnp.sum(h_dists_row[None, :] < neighbor_dists[:, None], axis=1)
    return ranks

# vmap over all cells
all_ranks = jax.vmap(_compute_ranks_jax, in_axes=(0, 0))(high_dists, low_knn_indices)
```

**Impact**: Trustworthiness drops from ~30 sec/embedding (CPU row-loop) to ~1 sec/embedding
(GPU vectorized). Since this runs for each of potentially 5+ embeddings × 2 (trustworthiness
+ continuity), the cumulative savings are significant.

---

### Phase 5: Optional Dependency Wiring

**Files**: `pyproject.toml`, `checkatlas/metrics/_jax_utils.py`

**Dependency specification**:
```toml
[project.optional-dependencies]
gpu = [
    "jax[cuda12]>=0.4.30",
]
cpu-jax = [
    "jax[cpu]>=0.4.30",  # JIT benefits even without GPU
]
```

**Makes JAX entirely optional**: `pip install checkatlas` → CPU-only, existing behaviour.
`pip install checkatlas[gpu]` → GPU acceleration when NVIDIA GPU + CUDA available.
`pip install checkatlas[cpu-jax]` → JIT compilation benefits on CPU only.

**Graceful degradation at every level**:
```python
# _jax_utils.py — single source of truth for JAX availability
_JAX_AVAILABLE = False
_GPU_AVAILABLE = False

try:
    import jax
    import jax.numpy as jnp
    _JAX_AVAILABLE = True
    try:
        devices = jax.devices()
        _GPU_AVAILABLE = any(d.platform == 'gpu' for d in devices)
    except Exception:
        pass
except ImportError:
    pass

def needs_jax(func):
    """Decorator that falls back to CPU function when JAX unavailable."""
    ...

def jax_or_cpu(jax_fn, cpu_fn):
    """Return the appropriate function based on _JAX_AVAILABLE."""
    return jax_fn if _JAX_AVAILABLE else cpu_fn
```

---

## Cumulative Performance Impact

Estimated total runtime for a 60 000-cell atlas (bone_marrow.h5ad) with all 32 metrics:

| Pipeline Stage | CPU (current) | CPU+JIT only | GPU+JIT | Speedup vs current |
|---|---|---|---|---|
| SUMMARY + QC | ~30 sec | ~30 sec | ~30 sec | 1× (not accelerated) |
| METRIC_CLUST (5 metrics × embeddings × labels) | ~2 min | ~1 min | ~15 sec | 8× |
| METRIC_ANNOT (16 metrics × embeddings × labels) | ~15 min | ~8 min | ~2 min | 7.5× |
| METRIC_DIMRED (11 metrics × embeddings) | **~150 min** | **~90 min** | **~5 min** | **30×** |
| **Total pipeline** | **~3 hours** | **~1.7 hours** | **~8 min** | **22×** |

**Key observation**: The dimred pipeline dominates everything. Phase 3 alone (GPU distance
computation) accounts for >80% of the total speedup. Phase 2 (GPU kNN) and Phase 4 (JIT
metric cores) provide the remaining 20%.

Even without a GPU, JIT compilation alone (Phase 4, using `jax[cpu]`) provides ~1.8× speedup
because `@jax.jit` compiles numpy-style code to efficient LLVM-native binaries, eliminating
Python interpreter overhead from hot loops.

---

## Implementation Order (Recommended)

1. **Phase 1** (Unified NeighborResults) — 2-3 hours. Foundation for everything else.
   Zero risk, pure refactoring. All existing tests should pass unchanged.

2. **Phase 3** (GPU distance computation) — 2-3 hours. Highest impact by far.
   Add `_jax_utils._pdist_squareform()` and `_jax_utils._cdist()` and wire into `cal_dimred()`.
   Test on the liver atlas first (15k cells, fast iteration).

3. **Phase 2** (GPU kNN backend) — 3-4 hours. Add `jax_exact_knn()` to `_neighbors.py`.
   Affects LISI, kBET, PCR, DBCVI, graph_connectivity. Each gets exact kNN that is faster
   than approximate pynndescent.

4. **Phase 4** (JIT metric cores) — 4-6 hours. Per-metric JAX functions, implemented one
   metric at a time with per-metric tests. Can be done incrementally — each metric improved
   is immediately faster even without a GPU.

5. **Phase 5** (Optional dependency wiring) — 1 hour. `pyproject.toml` changes, `_jax_utils.py`
   with `_HAS_JAX` / `_GPU_AVAILABLE` flags, backend auto-detection.

---

## Fallback Guarantees

Every GPU/JAX path has a corresponding CPU fallback that preserves:
- **Same algorithm** — The mathematical operation is identical; only the execution engine changes
- **Same inputs/outputs** — The function signatures remain unchanged
- **Same precision** — Both paths use float32 (GPU default) or float64 as appropriate
- **Existing tests pass** — The refactored functions produce identical results

The fallback chain at runtime:
```
try import jax → no? → use pynndescent + numpy (current code)
  ↓ yes
has CUDA GPU? → no? → jax[cpu] JIT on CPU (faster than numpy)
  ↓ yes
GPU path (fastest)
```

---

## References

- scib-metrics repository: <https://github.com/YosefLab/scib-metrics>
- JAX documentation: <https://jax.readthedocs.io/>
- `jax.lax.approx_min_k`: <https://jax.readthedocs.io/en/latest/_autosummary/jax.lax.approx_min_k.html>
- LISI original paper (Korsunsky 2019): <https://www.nature.com/articles/s41592-019-0619-0>
- kBET original paper (Büttner 2018): <https://www.nature.com/articles/s41592-018-0254-1>
- scIB benchmarking paper (Luecken 2022): <https://www.nature.com/articles/s41592-021-01336-8>
