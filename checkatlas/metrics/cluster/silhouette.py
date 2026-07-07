import numpy as np
from sklearn.metrics import pairwise_distances, silhouette_score
from tqdm import tqdm


_TRIANGULAR_DENSE_THRESHOLD = 30_000


def run(
    count_repr,
    annotations,
    sample_size=None,
    n_jobs=-1,
    random_state=42,
    verbose=True,
    batch_size=1000,
    precomputed_dists=None,
    method="approx",
):
    """
    Calculate the Average Silhouette Width (Silhouette Score).

    Four computation paths, resolved in priority order:

    1. **Sampling** (``sample_size``) -- sklearn's built-in random subsample
       path.  Fast approximate, quality depends on subsample size.

    2. **Centroid-Variance Approximation** (default, ``method="approx"``) --
       O(N·K·d) using the bluster formula:
       ``D_tilde(i, C_X) = sqrt(||x_i - mu_X||**2 + sigma_X**2)``
       where ``mu_X`` is the cluster centroid and ``sigma_X**2`` is the
       summed feature variance.  Runs in seconds even on 100k+ cells.
       *Takes priority over precomputed distances* -- avoids O(N^2)
       reads and GPU OOM crashes on large atlases.

    3. **Precomputed distances** (``precomputed_dists``, ``method="exact"``) --
       exact result from a precomputed N×N matrix.  TriangularMatrix support
       via sklearn (small) or ``_silhouette_from_triangular`` (large).

    4. **Exact raw** (``method="exact"``, no precomputed) -- full N×N
       distance matrix via chunked ``pairwise_distances``.  O(N^2·d).

    :param count_repr: array-like of shape (n_samples, n_features)
        Feature matrix containing the data points
    :param annotations: array-like of shape (n_samples,)
        Cluster labels for each sample
    :param sample_size: int, optional (default=None)
        The size of the sample to use when computing the Silhouette Coefficient
        on a random subset of the data. If None, no sampling is used.
    :param n_jobs: int, optional (default=-1)
        Number of parallel jobs for distance computation.
        -1 means using all processors.
    :param random_state: int, optional (default=42)
        Random state for reproducibility when sampling.
    :param verbose: bool, optional (default=True)
        Whether to show progress during computation.
    :param batch_size: int, optional (default=1000)
        Batch size for chunked distance computation (exact path only).
    :param precomputed_dists: np.ndarray or TriangularMatrix, optional
        Precomputed pairwise distance matrix.  Only used when
        ``method="exact"``.
    :param method: str, optional (default="approx")
        - ``"approx"`` -- Centroid-Variance Approximation (bluster method), fast
        - ``"exact"`` -- full N×N distance matrix (precomputed or raw)
    :return: float
        The mean Silhouette Coefficient over all samples.
        Range: [-1, 1], where higher values indicate better clustering.
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)
    n_samples = count_repr.shape[0]

    # Path 1: Sampling (fast approximate via sklearn)
    if sample_size is not None:
        return silhouette_score(
            count_repr,
            annotations,
            sample_size=sample_size,
            random_state=random_state,
        )

    # Path 2: Centroid-Variance Approximation (default, O(N·K·d))
    # Takes priority over precomputed distances — avoids O(N²) reads
    # and GPU OOM crashes on large atlases.
    if method == "approx":
        return _centroid_variance_silhouette(count_repr, annotations, verbose=verbose)

    # Path 3: Exact — use precomputed distances if available
    if precomputed_dists is not None:
        if hasattr(precomputed_dists, "to_dense"):
            if precomputed_dists.n <= _TRIANGULAR_DENSE_THRESHOLD:
                if verbose:
                    print("TriangularMatrix -- densifying for Silhouette (sklearn)...")
                dense = precomputed_dists.to_dense()
                return float(silhouette_score(dense, annotations, metric="precomputed"))
            else:
                if verbose:
                    print("TriangularMatrix -- computing Silhouette from precomputed distances...")
                from ..annot.average_silhouette_width import _silhouette_from_triangular

                return _silhouette_from_triangular(precomputed_dists, annotations, verbose=verbose)
        if verbose:
            print("Using precomputed distances for Silhouette...")
        return silhouette_score(precomputed_dists, annotations, metric="precomputed")

    # Path 4: Exact — full pairwise distance matrix from raw data
    if verbose:
        n_batches = (n_samples + batch_size - 1) // batch_size

        distances = np.zeros((n_samples, n_samples), dtype=np.float32)

        with tqdm(
            total=n_batches,
            desc="Computing pairwise distances",
            disable=not verbose,
        ) as pbar:
            for i in range(0, n_samples, batch_size):
                end_i = min(i + batch_size, n_samples)
                distances[i:end_i, :] = pairwise_distances(
                    count_repr[i:end_i], count_repr, n_jobs=n_jobs
                )
                pbar.update(1)

        np.fill_diagonal(distances, 0)
        if verbose:
            print("Computing silhouette scores...")
        return silhouette_score(distances, annotations, metric="precomputed")

    else:
        distances = pairwise_distances(count_repr, n_jobs=n_jobs)
        np.fill_diagonal(distances, 0)
        return silhouette_score(distances, annotations, metric="precomputed")


def _centroid_variance_silhouette(count_repr, annotations, verbose=True):
    """Compute silhouette via Centroid-Variance Approximation (bluster method).

    D_tilde(i, C_X) = sqrt(||x_i - mu_X||**2 + sigma_X**2)

    where ``mu_X`` is the cluster centroid and ``sigma_X**2`` is the
    summed feature variance (mean squared distance from points to
    their own centroid, multiplied by the number of features).

    This replaces the full O(N^2·d) pairwise distance matrix with an
    O(N·K·d) computation that is mathematically equivalent to the
    root-mean-squared distance from a point to all points in a cluster.

    Parameters
    ----------
    count_repr : ndarray of shape (n_samples, n_features)
    annotations : ndarray of shape (n_samples,)
    verbose : bool

    Returns
    -------
    float
        Mean silhouette coefficient in [-1, 1].
    """
    count_repr = np.asarray(count_repr)
    annotations = np.asarray(annotations)

    unique_labels = np.unique(annotations)
    n_clusters = len(unique_labels)

    if n_clusters < 2:
        return 0.0

    n_samples, n_features = count_repr.shape

    # ── Step 1: Cluster centroids and summed feature variances ──
    centroids = np.zeros((n_clusters, n_features), dtype=np.float64)
    variances = np.zeros(n_clusters, dtype=np.float64)
    cluster_sizes = np.zeros(n_clusters, dtype=np.int64)

    label_to_idx = {lbl: k for k, lbl in enumerate(unique_labels)}

    for lbl in unique_labels:
        k = label_to_idx[lbl]
        mask = annotations == lbl
        X_k = count_repr[mask]
        n_k = X_k.shape[0]
        cluster_sizes[k] = n_k

        centroid = X_k.mean(axis=0)
        centroids[k] = centroid

        if n_k > 1:
            # sigma^2 = sum_f Var_f(C) = (1 / n_k) * sum_j ||x_j - mu||^2
            variances[k] = np.sum((X_k - centroid) ** 2) / n_k
        else:
            variances[k] = 0.0

    # ── Step 2: Squared distances from all cells to all centroids ──
    # ||x_i - c_k||^2 = ||x_i||^2 + ||c_k||^2 - 2 x_i·c_k   (N, K)
    norms_x = np.sum(count_repr ** 2, axis=1, dtype=np.float64)
    norms_c = np.sum(centroids ** 2, axis=1, dtype=np.float64)

    dists_sq = norms_x[:, None] + norms_c[None, :] - 2.0 * count_repr @ centroids.T
    dists_sq = np.maximum(dists_sq, 0.0)

    # ── Step 3: Apply centroid-variance approximation ───────────
    # D_tilde(i, C_k) = sqrt(dist_to_centroid_k^2 + sigma_k^2)
    D_tilde = np.sqrt(dists_sq + variances[None, :])
    # (N, K)

    # ── Step 4: a(i) and b(i) per point ─────────────────────────
    own_idx = np.array([label_to_idx[lbl] for lbl in annotations])

    # a(i): D_tilde to own cluster
    a_values = D_tilde[np.arange(n_samples), own_idx]

    # b(i): minimum D_tilde to any OTHER cluster
    D_tilde[np.arange(n_samples), own_idx] = np.inf
    b_values = np.min(D_tilde, axis=1)

    # ── Step 5: per-point silhouette → mean ─────────────────────
    denom = np.maximum(a_values, b_values)
    mask = denom > 0
    s = np.zeros(n_samples)
    s[mask] = (b_values[mask] - a_values[mask]) / denom[mask]
    result = float(np.mean(s))

    if verbose:
        print(
            f"  Centroid-Variance Silhouette: "
            f"{n_samples:,} cells x {n_clusters} clusters "
            f"-> {result:.6f}"
        )

    return result
