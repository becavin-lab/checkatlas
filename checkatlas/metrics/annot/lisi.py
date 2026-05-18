import numpy as np
import pandas as pd
from scipy.sparse import issparse

# Module-level cache: (hash_of_X, n_neighbors) -> neighbour_idx
_knn_cache = {}


def _knn(X, k, n_jobs):
    """Fast approximate kNN using pynndescent, with caching per embedding."""
    flat = X.ravel()
    # Robust cache key: shape + hash of first/mid/last chunks
    n = len(flat)
    chunks = [
        flat[:2000].tobytes(),
        flat[max(0, n//2 - 1000):n//2 + 1000].tobytes(),
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


def run(X, label, perplexity=30, n_jobs=-1, verbose=True):
    """
    Calculate Local Inverse Simpson's Index (LISI) for batch mixing evaluation.

    Measures how well batches are mixed in local neighbourhoods.  For each
    cell we find its *k* nearest neighbours, compute the Simpson index of
    the batch labels in that neighbourhood, and take the inverse:

        LISI_cell = k² / Σ_b count(b)²

    The final score is scaled to [0, 1] where 1 indicates perfect mixing.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix (e.g. PCA embedding, integrated space).
    label : pandas Series, Index or array-like
        Batch or cell-type labels to evaluate mixing.
    perplexity : int, default=30
        **Kept for API compatibility — ignored internally.**
    n_jobs : int, default=-1
        Number of parallel jobs for kNN computation.  -1 uses all cores.
    verbose : bool, default=True
        Print progress information.

    Returns
    -------
    float
        Scaled LISI score.  Range [0, 1]; higher = better mixing.
    """
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X, dtype=np.float64)

    if isinstance(label, (pd.Series, pd.Index)):
        label = label.values
    else:
        label = np.asarray(label)

    n_samples = X.shape[0]
    if n_samples < 2:
        return 0.0

    unique_batches, batch_codes = np.unique(label, return_inverse=True)
    n_batches = len(unique_batches)
    if n_batches <= 1:
        return 0.0

    n_neighbors = min(90, n_samples - 1)

    if verbose:
        print(f"Computing LISI ({n_samples:,} samples, "
              f"{n_batches} labels, k={n_neighbors}, n_jobs={n_jobs})...")

    neighbour_idx = _knn(X, n_neighbors, n_jobs)   # (n_samples, n_neighbors)
    k_actual = neighbour_idx.shape[1]

    # ── vectorised batch counts ─────────────────────────────────
    neighbour_batch = batch_codes[neighbour_idx]      # (n_samples, k)
    offsets = np.arange(n_samples, dtype=np.int64) * n_batches
    flat = (neighbour_batch + offsets[:, None]).ravel()
    counts = np.bincount(flat, minlength=n_samples * n_batches)
    counts = counts.reshape(n_samples, n_batches).astype(np.float64)

    # Simpson index per cell:  Σ_b count_b² / k²
    sum_sq = np.sum(counts * counts, axis=1)
    sum_sq = np.maximum(sum_sq, 1e-10)

    lisi_scores = (k_actual * k_actual) / sum_sq

    # Scale to [0, 1]
    mean_lisi = np.mean(lisi_scores)
    scaled_lisi = (mean_lisi - 1.0) / (n_batches - 1.0)

    if verbose:
        print(f"  Scaled LISI = {scaled_lisi:.6f}")

    return float(scaled_lisi)


def ilisi_graph(adata, batch_key, k0=90, perplexity=None, scale=True,
                n_jobs=-1, verbose=False):
    """Integration LISI (iLISI) score — convenience wrapper."""
    X = adata.obsm.get('X_pca', adata.X)
    return run(X, adata.obs[batch_key], perplexity=perplexity,
               n_jobs=n_jobs, verbose=verbose)


def clisi_graph(adata, label_key, k0=90, perplexity=None, scale=True,
                n_jobs=-1, verbose=False):
    """Cell-type LISI (cLISI) score — convenience wrapper."""
    X = adata.obsm.get('X_pca', adata.X)
    return run(X, adata.obs[label_key], perplexity=perplexity,
               n_jobs=n_jobs, verbose=verbose)
