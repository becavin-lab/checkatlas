"""Tests for the scfm metric-collection layer.

These tests build a small synthetic AnnData and exercise the four
new metric modules (``scaling``, ``stability``, ``rare_types``,
``cross_domain``) and the ``cal_scfm`` orchestrator.

The tests require scanpy, anndata and scikit-learn to be installed.
If any of them are missing, the tests skip gracefully so the
package can still be imported in minimal environments.
"""


import numpy as np
import pandas as pd
import pytest


_scanpy_available = True
try:
    import scanpy as sc  # noqa: F401
    import anndata  # noqa: F401
    from sklearn.metrics import adjusted_rand_score  # noqa: F401
except ImportError:
    _scanpy_available = False


pytestmark = pytest.mark.skipif(
    not _scanpy_available,
    reason="scanpy / anndata / scikit-learn not installed",
)


@pytest.fixture
def small_adata() -> "anndata.AnnData":
    """Build a tiny AnnData with two embeddings and three obs columns.

    The data is engineered so that:

    * ``X_pca`` clusters the cells well (low intra-cluster distance).
    * ``X_scfm`` is a noisy version of ``X_pca`` (lower ARI).
    * ``celltype_author`` is imbalanced: Tcell is the common class
      (~200 cells), Bcell and NK are the rare classes (~50 each).
    * ``celltype_scfm_pred`` overpredicts the common class Tcell.
    * ``donor_id`` is a 2-batch factor.
    * ``tissue`` is a 2-domain factor.
    """
    import anndata

    rng = np.random.default_rng(42)
    n = 300
    n_genes = 50

    # Imbalanced labels: 200 Tcell, 50 Bcell, 50 NK.
    n_t, n_b, n_nk = 200, 50, 50
    y = np.concatenate([
        np.zeros(n_t, dtype=int),
        np.ones(n_b, dtype=int),
        np.full(n_nk, 2, dtype=int),
    ])
    rng.shuffle(y)

    # Three well-separated clusters in 5D.
    centers = np.array([[0, 0, 0, 0, 0], [5, 5, 5, 0, 0], [-5, -5, -5, 0, 0]])
    X = centers[y] + 0.5 * rng.standard_normal((n, 5))
    X_full = np.zeros((n, n_genes))
    X_full[:, :5] = X
    X_full[:, 5:] = rng.standard_normal((n, n_genes - 5)) * 0.1

    # Prediction: correctly classify Tcell, but always predict Tcell
    # for Bcell and NK (the rare classes). This reproduces the
    # scGPT CD14 overprediction finding.
    pred = np.where(y == 0, 0, 0)  # always Tcell except for half the Tcells
    pred[: n_t // 2] = 0

    obs = pd.DataFrame(
        {
            "celltype_author": pd.Categorical(
                ["Tcell" if c == 0 else "Bcell" if c == 1 else "NK" for c in y]
            ),
            "celltype_scfm_pred": pd.Categorical(
                ["Tcell" if p == 0 else "Bcell" if p == 1 else "NK" for p in pred]
            ),
            "donor_id": pd.Categorical(
                ["A" if i < n // 2 else "B" for i in range(n)]
            ),
            "tissue": pd.Categorical(
                ["PBMC" if i < n // 2 else "BM" for i in range(n)]
            ),
        }
    )

    adata = anndata.AnnData(X=X_full, obs=obs)
    adata.obsm["X_pca"] = X
    adata.obsm["X_scfm"] = X + 0.5 * rng.standard_normal(X.shape)
    adata.obsm["X_noisy"] = X + 2.0 * rng.standard_normal(X.shape)
    return adata


# ──────────────────────────────────────────────────────────────────
# scaling
# ──────────────────────────────────────────────────────────────────


def test_scaling_returns_long_dataframe(small_adata):
    from checkatlas.metrics.scfm.scaling import run

    df = run(
        small_adata,
        ["X_pca", "X_scfm"],
        ref_label="celltype_author",
        pred_label="celltype_scfm_pred",
        batch_key="donor_id",
        fractions=(0.5, 1.0),
    )
    assert not df.empty
    assert "Embedding" in df.columns
    assert "Metric" in df.columns
    assert "Fraction" in df.columns
    assert "Value" in df.columns
    assert set(df["Embedding"].unique()) == {"X_pca", "X_scfm"}


def test_scaling_contains_plateau_ratio(small_adata):
    from checkatlas.metrics.scfm.scaling import run

    df = run(
        small_adata,
        ["X_pca"],
        ref_label="celltype_author",
        pred_label="celltype_scfm_pred",
        fractions=(0.01, 0.10, 0.50, 1.00),
    )
    plateau = df[df["Fraction"] < 0]
    assert not plateau.empty
    assert "plateau_ratio" in str(plateau["Metric"].iloc[0]) or "ari" in str(
        plateau["Metric"].iloc[0]
    )


def test_scaling_skips_missing_embedding(small_adata):
    from checkatlas.metrics.scfm.scaling import run

    df = run(
        small_adata,
        ["X_does_not_exist", "X_pca"],
        ref_label="celltype_author",
        fractions=(1.0,),
    )
    assert "X_does_not_exist" not in df["Embedding"].unique()


# ──────────────────────────────────────────────────────────────────
# stability
# ──────────────────────────────────────────────────────────────────


def test_stability_returns_cv_per_metric(small_adata):
    from checkatlas.metrics.scfm.stability import run

    df = run(
        small_adata,
        ["X_pca"],
        ref_label="celltype_author",
        pred_label="celltype_scfm_pred",
        batch_key="donor_id",
        n_seeds=3,
    )
    assert not df.empty
    stats = set(df["Statistic"].unique())
    assert {"mean", "std", "cv"}.issubset(stats)


def test_stability_handles_missing_columns(small_adata):
    from checkatlas.metrics.scfm.stability import run

    df = run(
        small_adata,
        ["X_pca"],
        ref_label=None,
        pred_label=None,
        batch_key=None,
        n_seeds=2,
    )
    # No metrics evaluable -> empty DataFrame
    assert df.empty


# ──────────────────────────────────────────────────────────────────
# rare_types
# ──────────────────────────────────────────────────────────────────


def test_rare_types_detects_overprediction(small_adata):
    from checkatlas.metrics.scfm.rare_types import run

    df = run(
        small_adata,
        ["X_pca"],
        pred_label="celltype_scfm_pred",
        ref_label="celltype_author",
        rare_quantile=0.34,  # Bcell and NK are the rarer classes
    )
    assert not df.empty
    metrics = set(df["Metric"].unique())
    assert {"rare_f1", "common_f1", "rare_common_gap"}.issubset(metrics)
    # The overprediction of Tcell for NK should create a gap.
    gap_row = df[df["Metric"] == "rare_common_gap"].iloc[0]
    assert gap_row["Value"] > 0


def test_rare_types_handles_missing_columns(small_adata):
    from checkatlas.metrics.scfm.rare_types import run

    df = run(
        small_adata,
        ["X_pca"],
        pred_label="nonexistent",
        ref_label="celltype_author",
    )
    assert df.empty


# ──────────────────────────────────────────────────────────────────
# cross_domain
# ──────────────────────────────────────────────────────────────────


def test_cross_domain_returns_per_pair_ari(small_adata):
    from checkatlas.metrics.scfm.cross_domain import run

    df = run(
        small_adata,
        "X_pca",
        ref_label="celltype_author",
        domain_key="tissue",
        min_domain_size=20,
    )
    assert not df.empty
    assert "Train_Domain" in df.columns
    assert "Test_Domain" in df.columns
    assert "ARI" in df.columns
    # Self-pair ARI should be reasonable
    diag = df[df["Train_Domain"] == df["Test_Domain"]]
    assert not diag.empty


def test_cross_domain_handles_too_few_domains(small_adata):
    from checkatlas.metrics.scfm.cross_domain import run

    # Remove one domain so only one is left

    ad = small_adata.copy()
    ad.obs["tissue"] = pd.Categorical(["PBMC"] * ad.n_obs)
    df = run(
        ad,
        "X_pca",
        ref_label="celltype_author",
        domain_key="tissue",
        min_domain_size=10,
    )
    assert df.empty


# ──────────────────────────────────────────────────────────────────
# cal_scfm orchestrator
# ──────────────────────────────────────────────────────────────────


def test_cal_scfm_returns_long_table(small_adata):
    from checkatlas.metrics.scfm.run import cal_scfm

    df = cal_scfm(
        small_adata,
        scfm_embedding="X_scfm",
        baseline_embeddings=("X_pca",),
        ref_label="celltype_author",
        pred_label="celltype_scfm_pred",
        batch_key="donor_id",
        domain_key="tissue",
        atlas_name="test",
        scaling_fractions=(0.5, 1.0),
        n_seeds=2,
    )
    assert not df.empty
    assert "Atlas Name" in df.columns
    assert "Task" in df.columns
    tasks = set(df["Task"].unique())
    assert {"scfm_scaling", "scfm_stability", "scfm_rare_types", "scfm_cross_domain"}.issubset(tasks)
