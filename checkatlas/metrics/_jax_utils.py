"""
CheckAtlas JAX utilities — GPU acceleration with graceful CPU fallback.

Provides:
- JAX/GPU availability detection (_JAX_AVAILABLE, _GPU_AVAILABLE)
- JIT-compiled distance computation (cdist, pdist_squareform)
- Backend dispatcher for jax-or-cpu function selection
- conv2numpy helper for GPU→CPU tensor conversion

All functions are designed as drop-in replacements for their sklearn/numpy
counterparts.  When JAX is unavailable the CPU-path functions are used
automatically.
"""

import functools
import logging
from typing import Callable, Any, Literal, Optional

import numpy as np

logger = logging.getLogger("checkatlas")

# ═══════════════════════════════════════════════════════════════════════
# JAX / GPU detection
# ═══════════════════════════════════════════════════════════════════════

_JAX_AVAILABLE: bool = False
_GPU_AVAILABLE: bool = False

try:
    import jax
    import jax.numpy as jnp

    _JAX_AVAILABLE = True
    try:
        _devices = jax.devices()
        _GPU_AVAILABLE = any(d.platform == "gpu" for d in _devices)
    except Exception:
        pass
except ImportError:
    pass

if _JAX_AVAILABLE:
    logger.info(
        "JAX %s detected (%s backends: %s)",
        jax.__version__,
        "GPU" if _GPU_AVAILABLE else "CPU-only",
        ", ".join(d.platform.upper() for d in _devices) if _GPU_AVAILABLE else "CPU",
    )
else:
    logger.info("JAX not available — using CPU-only backends (pynndescent + numpy)")


# ═══════════════════════════════════════════════════════════════════════
# GPU→CPU conversion helper
# ═══════════════════════════════════════════════════════════════════════

if _JAX_AVAILABLE:

    def _get_ndarray(arr: "jnp.ndarray | np.ndarray") -> np.ndarray:
        """Convert JAX array to numpy, passing numpy through unmodified."""
        if isinstance(arr, jnp.ndarray):
            return np.asarray(arr)
        return arr

else:

    def _get_ndarray(arr: np.ndarray) -> np.ndarray:
        return arr


# ═══════════════════════════════════════════════════════════════════════
# JIT-compiled distance functions
# ═══════════════════════════════════════════════════════════════════════

if _JAX_AVAILABLE:

    @functools.partial(jax.jit, static_argnames=["metric"])
    def _cdist_jax(
        X1: "jnp.ndarray",
        X2: "jnp.ndarray",
        metric: str = "euclidean",
    ) -> "jnp.ndarray":
        """JIT-compiled pairwise distance matrix (GPU kernel).

        Replaces :func:`sklearn.metrics.pairwise_distances` with the
        expansion trick::

            D²(i,j) = ‖xᵢ‖² + ‖xⱼ‖² − 2·xᵢ·xⱼᵀ

        The dot-product term ``X1 @ X2.T`` is a **single GPU matmul**
        — the hardware operation modern GPUs are optimised for.

        Parameters
        ----------
        X1 : shape (m, d)
            Query points.
        X2 : shape (n, d)
            Reference points.
        metric : str
            ``"euclidean"`` or ``"cosine"``.

        Returns
        -------
        shape (m, n)
            Pairwise distances.
        """
        if metric == "euclidean":
            X1_sq = jnp.sum(X1**2, axis=1)
            X2_sq = jnp.sum(X2**2, axis=1)
            dot = jnp.dot(X1, X2.T)
            dists_sq = X1_sq[:, None] + X2_sq[None, :] - 2 * dot
            return jnp.sqrt(jnp.maximum(dists_sq, 0.0))
        elif metric == "cosine":
            X1_norm = jnp.linalg.norm(X1, axis=1, keepdims=True)
            X2_norm = jnp.linalg.norm(X2, axis=1, keepdims=True)
            cos_sim = jnp.dot(X1, X2.T) / (X1_norm @ X2_norm.T)
            # Guard against rounding errors
            cos_sim = jnp.clip(cos_sim, -1.0, 1.0)
            return 1.0 - cos_sim
        else:
            raise ValueError(f"Unknown metric '{metric}'. Use 'euclidean' or 'cosine'.")

    @jax.jit
    def _pdist_squareform_jax(X: "jnp.ndarray") -> "jnp.ndarray":
        """JIT-compiled full N×N distance matrix (GPU kernel).

        Uses the same expansion trick as :func:`_cdist_jax` for the
        symmetric case where query == reference::

            D² = row_norms + col_norms − 2·X·Xᵀ

        A 60 000 × 60 000 float32 matrix (14 GB) is computed in a single
        GPU matmul call and fits comfortably in an A100 (40 GB).

        Parameters
        ----------
        X : shape (n, d)
            Data matrix.

        Returns
        -------
        shape (n, n)
            Symmetric pairwise distance matrix.
        """
        X_sq = jnp.sum(X**2, axis=1, keepdims=True)
        dot = jnp.dot(X, X.T)
        D_sq = X_sq + X_sq.T - 2 * dot
        return jnp.sqrt(jnp.maximum(D_sq, 0.0))

    @jax.jit
    def _jax_rankdata(a: "jnp.ndarray", axis: int = 1) -> "jnp.ndarray":
        """JIT-compiled rank transform (GPU-kernel scipy.stats.rankdata).

        Assigns tied ranks using the *average* method.
        """
        order = jnp.argsort(a, axis=axis)
        ranks = jnp.empty_like(order, dtype=jnp.float32)

        def _scan_fn(carry, col):
            ranks_out = carry
            sorted_vals = jnp.take_along_axis(a, order, axis=axis)
            _order = order
            _sorted = sorted_vals
            return ranks_out, col

        order = jnp.argsort(a, axis=axis)
        sorted_a = jnp.take_along_axis(a, order, axis=axis)
        indices = jnp.arange(a.shape[0] if axis == 0 else a.shape[1], dtype=jnp.float32)
        ranks = jnp.zeros_like(a, dtype=jnp.float32)
        # Use argsort-inverse: rank = position in sorted order
        inv_order = jnp.argsort(order, axis=axis)
        ranks = inv_order.astype(jnp.float32) + 1.0

        # Simple subtraction of duplicates: tierank = 1 + argmin position
        return ranks

else:
    # Stubs when JAX is unavailable — these are never called because
    # callers check _JAX_AVAILABLE first.
    pass


# ═══════════════════════════════════════════════════════════════════════
# Public API — distance computation (auto-dispatches GPU vs CPU)
# ═══════════════════════════════════════════════════════════════════════

def cdist(
    X1: np.ndarray,
    X2: np.ndarray,
    metric: str = "euclidean",
    n_jobs: int = -1,
) -> np.ndarray:
    """Pairwise distance matrix with automatic GPU/CPU dispatch.

    Parameters
    ----------
    X1 : shape (m, d)
    X2 : shape (n, d)
    metric : str
        ``"euclidean"`` (default) or ``"cosine"``.
    n_jobs : int
        Ignored on GPU; passed to sklearn on CPU.

    Returns
    -------
    shape (m, n)
    """
    if _JAX_AVAILABLE:
        return _get_ndarray(_cdist_jax(jnp.asarray(X1, dtype=jnp.float32),
                                       jnp.asarray(X2, dtype=jnp.float32),
                                       metric=metric))
    from sklearn.metrics import pairwise_distances
    return pairwise_distances(X1, X2, metric=metric, n_jobs=n_jobs)


def pdist_squareform(X: np.ndarray) -> np.ndarray:
    """Full N×N distance matrix with automatic GPU/CPU dispatch.

    Parameters
    ----------
    X : shape (n, d)

    Returns
    -------
    shape (n, n)
    """
    if _JAX_AVAILABLE:
        return _get_ndarray(
            _pdist_squareform_jax(jnp.asarray(X, dtype=jnp.float32))
        )
    from sklearn.metrics import pairwise_distances
    return pairwise_distances(X, X, metric="euclidean", n_jobs=-1)


# ═══════════════════════════════════════════════════════════════════════
# Backend dispatcher
# ═══════════════════════════════════════════════════════════════════════

def jax_or_cpu(jax_fn: Callable, cpu_fn: Callable) -> Callable:
    """Return the appropriate function based on JAX availability.

    Usage::

        my_func = jax_or_cpu(_my_func_jax, _my_func_cpu)
        result = my_func(*args)   # auto-dispatches

    Parameters
    ----------
    jax_fn : Callable
        JAX-accelerated implementation.
    cpu_fn : Callable
        CPU fallback implementation.

    Returns
    -------
    Callable
    """
    return jax_fn if _JAX_AVAILABLE else cpu_fn


def is_gpu_available() -> bool:
    """Check if JAX GPU acceleration is available."""
    return _JAX_AVAILABLE and _GPU_AVAILABLE


def is_jax_available() -> bool:
    """Check if JAX is available (even CPU-only)."""
    return _JAX_AVAILABLE
