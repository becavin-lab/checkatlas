"""scFM-specific metrics.

Each module exposes a ``run(...)`` function with a uniform signature:

    run(adata, embeddings, **kwargs) -> pd.DataFrame

The DataFrame is in long format with at least
``Embedding, Metric, Value`` columns so it can be stacked onto the
existing CheckAtlas metric table.

All modules accept an optional ``preprocess_context`` (Becavin
comment 2) and never re-load adata, re-detect columns, or re-compute
kNN / distance matrices that the context already provides.
"""

from . import cross_domain, rare_types, scaling, stability

__all__ = ["scaling", "stability", "rare_types", "cross_domain"]
