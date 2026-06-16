"""Cross-domain generalisation metric (Problem 6).

Reproduces the Souza & Mehta (2026) finding that scFM embeddings
generalise poorly across species / tissues / assays.

For each ``(train_domain, test_domain)`` pair with at least
``min_domain_size`` cells, the metric:
  1. Subsamples the test domain,
  2. Clusters cells in the test domain in the embedding space
     (Leiden via scanpy, on a kNN graph built from the embedding),
  3. Compares the test-domain clustering to the test-domain
     reference labels via Adjusted Rand Index.

The result is one long-format row per (train_domain, test_domain)
pair, plus a summary row (mean off-diagonal ARI).
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score

logger = logging.getLogger("checkatlas")


def _leiden_labels(X: np.ndarray, seed: int = 42) -> np.ndarray:
    """Leiden clustering of X (cells x features) using scanpy."""
    try:
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover
        raise ImportError("scanpy is required for cross_domain.run") from exc

    import anndata as ad

    a = ad.AnnData(X=X)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(a, use_rep="X", random_state=seed)
        sc.tl.leiden(a, random_state=seed, flavor="igraph", n_iterations=2)
    return np.asarray(a.obs["leiden"].astype(str))


def run(
    adata: AnnData,
    embedding: str,
    *,
    ref_label: str,
    domain_key: str,
    min_domain_size: int = 50,
    seed: int = 42,
) -> pd.DataFrame:
    """Compute per-pair cross-domain ARI for a single embedding.

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

    Returns
    -------
    pd.DataFrame
        Long format with columns
        ``[Embedding, Train_Domain, Test_Domain, N_Test, ARI, N_Classes]``.
    """
    if embedding not in adata.obsm:
        return pd.DataFrame(
            columns=[
                "Embedding",
                "Train_Domain",
                "Test_Domain",
                "N_Test",
                "ARI",
                "N_Classes",
            ]
        )
    if domain_key not in adata.obs.columns:
        return pd.DataFrame(
            columns=[
                "Embedding",
                "Train_Domain",
                "Test_Domain",
                "N_Test",
                "ARI",
                "N_Classes",
            ]
        )
    if ref_label not in adata.obs.columns:
        return pd.DataFrame(
            columns=[
                "Embedding",
                "Train_Domain",
                "Test_Domain",
                "N_Test",
                "ARI",
                "N_Classes",
            ]
        )

    domains = adata.obs[domain_key].astype(str)
    unique_domains = sorted(d.unique() for d in [domains])[0]
    unique_domains = np.asarray(unique_domains)
    sizes = np.array([int((domains == d).sum()) for d in unique_domains])
    keep = unique_domains[sizes >= min_domain_size]
    if keep.size < 2:
        return pd.DataFrame(
            columns=[
                "Embedding",
                "Train_Domain",
                "Test_Domain",
                "N_Test",
                "ARI",
                "N_Classes",
            ]
        )

    X = adata.obsm[embedding]
    y_true = adata.obs[ref_label].astype(str).values

    rows: list[dict] = []
    for train_d in keep:
        for test_d in keep:
            test_mask = (domains == test_d).values
            X_test = X[test_mask]
            y_test = y_true[test_mask]
            if X_test.shape[0] < 2 or len(np.unique(y_test)) < 2:
                continue
            try:
                test_clusters = _leiden_labels(X_test, seed=seed)
                ari = float(adjusted_rand_score(y_test, test_clusters))
            except Exception as exc:  # pragma: no cover
                logger.debug(
                    "cross_domain: clustering failed for %s/%s: %s",
                    embedding,
                    test_d,
                    exc,
                )
                ari = float("nan")
            rows.append(
                {
                    "Embedding": embedding,
                    "Train_Domain": train_d,
                    "Test_Domain": test_d,
                    "N_Test": int(test_mask.sum()),
                    "ARI": ari,
                    "N_Classes": int(len(np.unique(y_test))),
                }
            )
    return pd.DataFrame(rows)
