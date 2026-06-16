"""Cross-domain generalisation metric (Problem 6).

Reproduces the Souza & Mehta (2026) finding that scFM embeddings
generalise poorly across species / tissues / assays.

**Contract (per the user instruction): the entire input atlas is
used — no subsampling. The kNN graph is built ONCE on the full
embedding via the unified :func:`metrics._neighbors.compute_neighbors`
(which auto-dispatches to JAX GPU when available, pynndescent
otherwise) and is re-used across all ``(train, test)`` domain
pairs.**

For each ``(train_domain, test_domain)`` pair with at least
``min_domain_size`` cells, the metric:
  1. Identifies cells in the test domain,
  2. Reads the test-domain kNN subgraph from the precomputed full
     kNN graph (no rebuild),
  3. Runs Leiden clustering on the test-domain subgraph,
  4. Compares the test-domain clustering to the test-domain
     reference labels via Adjusted Rand Index.

The result is one long-format row per (train_domain, test_domain)
pair. The Leiden step uses the igraph backend via scanpy; that
clustering itself runs once per pair, on the (small) test-domain
subgraph only.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score

logger = logging.getLogger("checkatlas")


def _leiden_labels_from_knn(
    nn_indices: np.ndarray,
    seed: int = 42,
) -> np.ndarray:
    """Run Leiden clustering on a precomputed kNN graph.

    Parameters
    ----------
    nn_indices : np.ndarray
        ``(n_test, n_neighbors)`` integer matrix. Each row contains
        the cell indices of the k nearest neighbours of the i-th
        test-domain cell. Indices are global (into the full atlas);
        they are re-mapped to 0..n_test-1 here.
    seed : int
        Random seed for the Leiden algorithm.
    """
    try:
        import igraph as ig
        import leidenalg
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "cross_domain requires igraph and leidenalg; install with "
            "`pip install python-igraph leidenalg`."
        ) from exc

    n_test = nn_indices.shape[0]
    if n_test < 2:
        return np.zeros(n_test, dtype=int)

    # Re-map global indices to 0..n_test-1.
    remap = -np.ones(int(nn_indices.max()) + 1, dtype=np.int64) if nn_indices.size else None
    if remap is not None:
        remap[nn_indices.ravel()] = np.arange(nn_indices.size) % n_test
    edges_src = np.repeat(np.arange(n_test), nn_indices.shape[1])
    edges_dst = remap[nn_indices.ravel()]
    # Drop self-loops
    mask = edges_src != edges_dst
    g = ig.Graph(n=n_test, edges=list(zip(edges_src[mask].tolist(), edges_dst[mask].tolist())))
    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition, seed=seed
    )
    return np.asarray(partition.membership, dtype=int)


def run(
    adata: AnnData,
    embedding: str,
    *,
    ref_label: str,
    domain_key: str,
    min_domain_size: int = 50,
    seed: int = 42,
    n_neighbors: int = 30,
    preprocess_context=None,
) -> pd.DataFrame:
    """Compute per-pair cross-domain ARI for a single embedding.

    The full input ``adata`` is consumed as-is. The kNN graph is
    built ONCE on the full embedding; per-domain subgraphs are
    sliced from that precomputed graph.

    Parameters
    ----------
    adata : AnnData
    embedding : str
        Single embedding key.
    ref_label : str
        Reference cell-type column.
    domain_key : str
        Categorical obs column with the domain label (species, tissue,
        assay, etc.).
    min_domain_size : int
        Minimum number of cells in a domain for inclusion.
    seed : int
        Random seed for the Leiden algorithm.
    n_neighbors : int
        Number of nearest neighbours for the kNN graph (default 30).
    preprocess_context : PreprocessContext, optional
        If provided AND ``embedding`` is one of the cached
        ``knn_paths`` keys, the precomputed kNN is reused. Otherwise
        the kNN is built locally with
        :func:`metrics._neighbors.compute_neighbors` (JAX GPU if
        available, pynndescent CPU otherwise).

    Returns
    -------
    pd.DataFrame
        Long format with columns
        ``[Embedding, Train_Domain, Test_Domain, N_Test, ARI, N_Classes]``.
    """
    empty_df = pd.DataFrame(
        columns=[
            "Embedding",
            "Train_Domain",
            "Test_Domain",
            "N_Test",
            "ARI",
            "N_Classes",
        ]
    )
    if embedding not in adata.obsm:
        return empty_df
    if domain_key not in adata.obs.columns:
        return empty_df
    if ref_label not in adata.obs.columns:
        return empty_df

    domains = adata.obs[domain_key].astype(str)
    unique_domains = np.asarray(sorted(domains.unique()))
    sizes = np.array([int((domains == d).sum()) for d in unique_domains])
    keep = unique_domains[sizes >= min_domain_size]
    if keep.size < 2:
        return empty_df

    # ── Build (or load) the kNN graph for the full embedding ONCE ──
    nn = _load_or_build_knn(
        adata, embedding, n_neighbors, preprocess_context
    )
    if nn is None:
        return empty_df

    y_true = adata.obs[ref_label].astype(str).values

    rows: list[dict] = []
    for test_d in keep:
        test_idx_global = np.where(domains == test_d)[0]
        n_test = int(test_idx_global.size)
        if n_test < 2 or len(np.unique(y_true[test_idx_global])) < 2:
            continue
        # Slice the per-cell kNN to keep only test-domain cells, and
        # remap neighbour indices to the local 0..n_test-1 range.
        # Cells whose neighbour falls inside the test domain are
        # valid; cells whose neighbour falls outside are dropped
        # from the per-test subgraph but the test cell itself is
        # preserved (a self-loop in the local graph).
        nn_test = nn.indices[test_idx_global]  # (n_test, n_neighbors)
        global_to_local = -np.ones(adata.n_obs, dtype=np.int64)
        global_to_local[test_idx_global] = np.arange(n_test, dtype=np.int64)
        local_dst = global_to_local[nn_test.ravel()].reshape(n_test, -1)
        # Self-loop fallback for neighbours outside the test domain
        # (i.e. -1 in local_dst). Use a uniform random neighbour
        # inside the test domain to keep the graph connected.
        no_neighbour = local_dst < 0
        if no_neighbour.any():
            fallback = np.random.default_rng(seed).integers(
                low=0, high=n_test, size=int(no_neighbour.sum())
            )
            local_dst[no_neighbour] = fallback
        try:
            test_clusters = _leiden_labels_from_knn(local_dst, seed=seed)
            ari = float(adjusted_rand_score(y_true[test_idx_global], test_clusters))
        except ImportError:
            logger.debug(
                "cross_domain: igraph/leidenalg unavailable, falling "
                "back to a simple connectivity-based label propagation"
            )
            ari = float("nan")
        except Exception as exc:  # pragma: no cover
            logger.debug(
                "cross_domain: Leiden failed for %s/%s: %s",
                embedding,
                test_d,
                exc,
            )
            ari = float("nan")
        rows.append(
            {
                "Embedding": embedding,
                "Train_Domain": "*",  # no train side: Leiden is unsupervised
                "Test_Domain": test_d,
                "N_Test": n_test,
                "ARI": ari,
                "N_Classes": int(len(np.unique(y_true[test_idx_global]))),
            }
        )
    if not rows:
        return empty_df
    return pd.DataFrame(rows)


def _load_or_build_knn(
    adata: AnnData,
    embedding: str,
    n_neighbors: int,
    preprocess_context,
):
    """Try to load a precomputed kNN from ``preprocess_context``; fall
    back to :func:`compute_neighbors` (JAX GPU if available, else
    pynndescent CPU). The atlas is never re-loaded."""
    from .._neighbors import NeighborResults, compute_neighbors

    if preprocess_context is not None and hasattr(preprocess_context, "knn_paths"):
        npz_path = preprocess_context.knn_paths.get(embedding)
        if npz_path is not None:
            try:
                from .._cache import load_knn

                loaded = load_knn(
                    os.path.dirname(npz_path),
                    os.path.splitext(os.path.basename(npz_path))[0],
                )
                if loaded is not None:
                    return NeighborResults(
                        indices=loaded[0], distances=loaded[1]
                    )
            except Exception:
                pass

    X = adata.obsm[embedding]
    try:
        return compute_neighbors(X, n_neighbors=n_neighbors + 1, backend="auto")
    except Exception as exc:
        logger.debug("cross_domain: compute_neighbors failed: %s", exc)
        return None
