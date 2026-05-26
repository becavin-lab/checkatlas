import numpy as np
from scipy.stats import chi2
from scipy.sparse import issparse
from anndata import AnnData

# Module-level cache: (hash_key, k) -> neighbour indices
_knn_cache = {}


def _knn(X, k, n_jobs):
    """Fast approximate kNN using pynndescent, with caching per embedding."""
    flat = X.ravel()
    n = len(flat)
    chunks = [
        flat[:2000].tobytes(),
        flat[max(0, n // 2 - 1000):n // 2 + 1000].tobytes(),
        flat[-2000:].tobytes() if n > 2000 else b'',
    ]
    key = (X.shape, hash(b''.join(chunks)), k)
    if key in _knn_cache:
        return _knn_cache[key]
    try:
        from pynndescent import NNDescent
        index = NNDescent(X, n_neighbors=k + 1, metric='euclidean',
                          n_jobs=n_jobs, random_state=42)
        idx, _ = index.neighbor_graph
    except ImportError:
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm='auto',
                                n_jobs=n_jobs)
        nbrs.fit(X)
        _, idx = nbrs.kneighbors(X)
    result = idx[:, 1:]  # drop self
    _knn_cache[key] = result
    return result


def _clear_knn_cache():
    """Clear the kNN cache (call at end of pipeline to free memory)."""
    _knn_cache.clear()


def run(X, labels, k=25, alpha=0.05, n_jobs=-1, verbose=True):
    """
    Calculate kBET rejection rate.

    kBET tests whether the batch composition in local neighbourhoods
    matches the global batch distribution via a chi-squared test.  Lower
    rejection rates indicate better batch mixing.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix (e.g. PCA embedding, UMAP, scVI latent space).
    labels : array-like of shape (n_samples,)
        Batch labels to evaluate mixing.
    k : int, default=25
        Number of nearest neighbours to consider.
    alpha : float, default=0.05
        Significance level for the chi-squared test.
    n_jobs : int, default=-1
        Number of parallel jobs for kNN.  -1 uses all cores.
    verbose : bool, default=True
        Print progress information.

    Returns
    -------
    float
        kBET rejection rate.  Range [0, 1]; lower = better mixing.
    """
    # ── backward compat: if first arg is AnnData, extract X_pca ──
    if isinstance(X, AnnData):
        adata = X
        batch_label = labels if isinstance(labels, str) else "batch"
        if 'X_pca' in adata.obsm:
            X_arr = np.asarray(adata.obsm['X_pca'], dtype=np.float64)
        elif issparse(adata.X):
            X_arr = adata.X.toarray().astype(np.float64)
        else:
            X_arr = np.asarray(adata.X, dtype=np.float64)
        labels_arr = np.asarray(adata.obs[batch_label])
        return run(X_arr, labels_arr, k=k, alpha=alpha,
                   n_jobs=n_jobs, verbose=verbose)

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X, dtype=np.float64)

    labels = np.asarray(labels)

    n_samples = X.shape[0]
    if n_samples < 2:
        return 0.0

    unique_batches, batch_counts = np.unique(labels, return_counts=True)
    n_batches = len(unique_batches)
    if n_batches <= 1:
        return 0.0

    # ── fast approximate kNN ──────────────────────────────────────
    k = min(k, n_samples - 1)

    if verbose:
        print(f"Computing kBET ({n_samples:,} samples, "
              f"k={k}, n_jobs={n_jobs})...")

    neighbour_idx = _knn(X, k, n_jobs)   # (n_samples, k)
    k_actual = neighbour_idx.shape[1]

    # ── vectorised chi-squared ─────────────────────────────────────
    expected_freqs = batch_counts / n_samples           # (n_batches,)
    expected_counts = expected_freqs * k_actual
    expected_counts = np.maximum(expected_counts, 1e-10)

    # Encode batches → integer codes
    batch_to_code = {b: i for i, b in enumerate(unique_batches)}
    batch_codes = np.array([batch_to_code[b] for b in labels],
                           dtype=np.intp)

    neighbour_batch = batch_codes[neighbour_idx]        # (n_samples, k)

    offsets = np.arange(n_samples, dtype=np.int64) * n_batches
    flat = (neighbour_batch + offsets[:, None]).ravel()
    observed_flat = np.bincount(flat, minlength=n_samples * n_batches)
    observed = observed_flat.reshape(n_samples, n_batches).astype(np.float64)

    # Chi-squared statistic per cell
    chi2_stats = np.sum((observed - expected_counts) ** 2 / expected_counts,
                        axis=1)

    df = n_batches - 1
    p_values = 1.0 - chi2.cdf(chi2_stats, df)

    rejections = np.sum(p_values < alpha)
    rejection_rate = rejections / n_samples

    if verbose:
        print(f"  kBET rejection rate = {rejection_rate:.4f} "
              f"({rejections:,}/{n_samples:,} cells rejected)")

    return float(rejection_rate)


# ═══════════════════════════════════════════════════════════════════════
# JAX GPU-accelerated kBET (Phase 4b)
# ═══════════════════════════════════════════════════════════════════════

try:
    import functools
    import jax
    import jax.numpy as jnp
    _HAS_JAX_KBET = True
except ImportError:
    _HAS_JAX_KBET = False


if _HAS_JAX_KBET:

    @functools.partial(jax.jit, static_argnums=(2,))
    def _kbet_chi2_jax(
        neigh_batch_ids: jnp.ndarray,
        batches: jnp.ndarray,
        n_batches: int,
    ) -> jnp.ndarray:
        """JIT-compiled chi-squared test per cell (GPU kernel).

        Fuses bincount → expected_counts → chi² → p-value into one
        GPU kernel, avoiding Python interpreter overhead between steps.
        """
        expected_freq = jnp.bincount(batches, length=n_batches) / batches.shape[0]
        dof = n_batches - 1

        observed_counts = jax.vmap(
            functools.partial(jnp.bincount, length=n_batches)
        )(neigh_batch_ids)
        expected_counts = expected_freq * neigh_batch_ids.shape[1]
        test_stats = jnp.sum(
            jnp.square(observed_counts - expected_counts) / expected_counts,
            axis=1,
        )
        p_values = 1.0 - jax.scipy.special.gammainc(dof / 2.0, test_stats / 2.0)
        return p_values


def run_with_neighbors(neighbor_results, labels, alpha=0.05, verbose=True):
    """kBET using precomputed kNN (NeighborResults).

    Uses JAX GPU acceleration when available; falls back to numpy.

    Parameters
    ----------
    neighbor_results : NeighborResults
        Precomputed kNN from :func:`checkatlas.metrics._neighbors.compute_neighbors`.
    labels : array-like
        Batch labels.
    alpha : float
        Significance level.
    verbose : bool

    Returns
    -------
    float
        kBET rejection rate [0, 1]; lower = better mixing.
    """
    from .._neighbors import NeighborResults

    labels = np.asarray(labels)

    n_samples = neighbor_results.n_samples
    if n_samples < 2:
        return 0.0

    unique_batches, batch_counts = np.unique(labels, return_counts=True)
    n_batches = len(unique_batches)
    if n_batches <= 1:
        return 0.0

    neighbour_idx = neighbor_results.indices
    k_actual = neighbour_idx.shape[1]

    batch_to_code = {b: i for i, b in enumerate(unique_batches)}
    batch_codes = np.array([batch_to_code[b] for b in labels], dtype=np.intp)
    neighbour_batch = batch_codes[neighbour_idx]

    if _HAS_JAX_KBET:
        jax_neigh = jnp.asarray(neighbour_batch, dtype=jnp.int32)
        jax_batches = jnp.asarray(batch_codes, dtype=jnp.int32)
        p_values = _kbet_chi2_jax(jax_neigh, jax_batches, n_batches)
        p_values = np.asarray(p_values)
    else:
        expected_freqs = batch_counts / n_samples
        expected_counts = expected_freqs * k_actual
        expected_counts = np.maximum(expected_counts, 1e-10)

        offsets = np.arange(n_samples, dtype=np.int64) * n_batches
        flat = (neighbour_batch + offsets[:, None]).ravel()
        observed_flat = np.bincount(flat, minlength=n_samples * n_batches)
        observed = observed_flat.reshape(n_samples, n_batches).astype(np.float64)

        chi2_stats = np.sum(
            (observed - expected_counts) ** 2 / expected_counts, axis=1
        )
        p_values = 1.0 - chi2.cdf(chi2_stats, n_batches - 1)

    rejections = np.sum(p_values < alpha)
    rejection_rate = rejections / n_samples

    if verbose:
        print(f"  kBET rejection rate = {rejection_rate:.4f} "
              f"({rejections:,}/{n_samples:,} cells rejected)")

    return float(rejection_rate)
