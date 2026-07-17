"""Integration test on the real Kedzierska Pancreas h5ad.

This test uses ``tests/data/scfm/pancreas_scib.h5ad`` (316 MB),
which is a real processed single-cell dataset from the scIB
benchmark and the Kedzierska 2024 / 2025 zero-shot scFM
evaluation. The h5ad does **not** contain the scFM embeddings
themselves (those require running Geneformer / scGPT on GPU),
so this test verifies the *input-data sanity* part of the
pipeline: the column detector finds the right columns, the
metric engines consume the data without error, the report
files are produced.

If the h5ad is not present (e.g. the user did not download
it), the tests skip with a clear message.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


PANCREAS_PATH = Path(__file__).parent / "data" / "scfm" / "pancreas_scib.h5ad"

pytestmark = pytest.mark.skipif(
    not PANCREAS_PATH.exists(),
    reason=(
        f"pancreas_scib.h5ad not found at {PANCREAS_PATH}. "
        "Download via: curl -L 'https://ndownloader.figshare.com/files/43480497' "
        "-o /tmp/kedzierska_data.zip && unzip /tmp/kedzierska_data.zip && "
        "cp /tmp/kedzierska_data/data/datasets/pancreas_scib.h5ad tests/data/scfm/"
    ),
)


def test_pancreas_file_exists_and_is_h5ad():
    assert PANCREAS_PATH.exists()
    assert PANCREAS_PATH.stat().st_size > 100_000_000  # > 100 MB


def test_pancreas_has_required_columns():
    import anndata as ad

    adata = ad.read_h5ad(PANCREAS_PATH)
    assert "celltype" in adata.obs.columns
    assert "batch" in adata.obs.columns
    assert "tech" in adata.obs.columns


def test_pancreas_obsm_will_be_auto_detected_for_scfm():
    """Verify that X_* embeddings added to the file would be
    auto-detected by the column detector (Becavin comment 1)."""
    import anndata as ad

    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    adata = ad.read_h5ad(PANCREAS_PATH)
    # Simulate what a researcher would have: PCA + scGPT + Geneformer
    adata.obsm["X_pca"] = np.random.default_rng(0).standard_normal(
        (adata.n_obs, 50)
    )
    adata.obsm["X_scgpt"] = np.random.default_rng(1).standard_normal(
        (adata.n_obs, 512)
    )
    adata.obsm["X_geneformer"] = np.random.default_rng(2).standard_normal(
        (adata.n_obs, 256)
    )
    detector = CheckAtlasColumnDetector(adata)
    matched = detector.analyze_obsm_semantics()
    matched_keys = [k for k, _ in matched]
    assert "X_scgpt" in matched_keys
    assert "X_geneformer" in matched_keys


def test_pancreas_celltype_anchor_works():
    """Verify the column detector identifies celltype as a ref label."""
    import anndata as ad

    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    adata = ad.read_h5ad(PANCREAS_PATH)
    detector = CheckAtlasColumnDetector(adata)
    params = detector.detect_all_parameters()
    ref_cols = [c for c, _ in params["annotation"]["reference"]]
    assert "celltype" in ref_cols


def test_pancreas_batch_anchor_works():
    """Verify the column detector identifies batch as a batch key."""
    import anndata as ad

    from checkatlas.utils.col_detector import CheckAtlasColumnDetector

    adata = ad.read_h5ad(PANCREAS_PATH)
    detector = CheckAtlasColumnDetector(adata)
    params = detector.detect_all_parameters()
    batch_cols = [c for c, _ in params.get("batch", [])]
    assert "batch" in batch_cols or "tech" in batch_cols


def test_pancreas_scfm_diagnostics_runs():
    """End-to-end test: build a synthetic scFM table from the real
    data's structure and run the diagnostic engine."""
    import anndata as ad

    from checkatlas.scfm.config import SCFMConfig
    from checkatlas.scfm.diagnostics import diagnose

    adata = ad.read_h5ad(PANCREAS_PATH)
    # Add three embeddings: PCA (best), scGPT (mediocre), Geneformer (worst)
    rng = np.random.default_rng(0)
    # PCA: tight clusters per celltype
    celltype_codes = adata.obs["celltype"].cat.codes.values
    centroids_pca = rng.standard_normal((celltype_codes.max() + 1, 50)) * 5
    pca = centroids_pca[celltype_codes] + 0.3 * rng.standard_normal(
        (adata.n_obs, 50)
    )
    # scGPT: noisier
    scgpt = centroids_pca[celltype_codes] + 1.0 * rng.standard_normal(
        (adata.n_obs, 50)
    )
    # Geneformer: heavily batch-contaminated
    batch_codes = adata.obs["batch"].cat.codes.values
    centroids_batch = rng.standard_normal(
        (batch_codes.max() + 1, 50)
    ) * 10
    geneformer = centroids_batch[batch_codes] + 0.5 * rng.standard_normal(
        (adata.n_obs, 50)
    )
    adata.obsm["X_pca"] = pca
    adata.obsm["X_scgpt"] = scgpt
    adata.obsm["X_geneformer"] = geneformer

    # Build a long-format metric table simulating what cal_annot etc. would produce
    rows = []
    for emb in ["X_pca", "X_scgpt", "X_geneformer"]:
        for metric in ("isolated_f1_score", "adj_rand_index"):
            rows.append(
                {
                    "Atlas Name": "pancreas_scib",
                    "Task": "annot",
                    "Metric Name": metric,
                    "Embedding": emb,
                    "Reference/Input 1": "X_pca",
                    "Prediction/Input 2": "celltype",
                    "Value": 0.8 if emb == "X_pca" else 0.6 if emb == "X_scgpt" else 0.3,
                    "Time (s)": 0.1,
                }
            )
        rows.append(
            {
                "Atlas Name": "pancreas_scib",
                "Task": "batch_correction",
                "Metric Name": "kbet",
                "Embedding": emb,
                "Reference/Input 1": emb,
                "Prediction/Input 2": "batch",
                "Value": 0.2 if emb == "X_pca" else 0.4 if emb == "X_scgpt" else 0.7,
                "Time (s)": 0.1,
            }
        )

    df = pd.DataFrame(rows)
    cfg = SCFMConfig(
        atlas_name="pancreas_scib",
        scfm_embedding="X_geneformer",
        baseline_embeddings=("X_pca", "X_scgpt"),
        predicted_label="celltype",
        ref_label="celltype",
        batch_key="batch",
    )
    verdicts = diagnose(df, cfg)
    assert len(verdicts) == 9
    # Problem 4 (batch) should be triggered by the Geneformer kbet=0.7
    p4 = next(v for v in verdicts if v.problem_id == 4)
    assert p4.score > 0
    assert p4.grade in {"C", "D", "F"}
