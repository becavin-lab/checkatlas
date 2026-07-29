"""Download TF paper evaluation datasets from CELLxGENE Census with TF embeddings.

Datasets:
  - tabula_sapiens_2: Tabula Sapiens 2.0 held-out donors (TSP17-TSP30)
  - covid19_lung: COVID-19 lung dataset
  - pbmc_10k: PBMC dataset (alternative source for scGPT benchmarking)

All saved with:
  - X_tf_sapiens (2048-dim, from Census)
  - X_tf_exemplar (2048-dim, from Census) -- where available
  - X_scvi_census (50-dim, from Census)
  - cell_type, donor_id, assay, dataset_id in obs
"""

from __future__ import annotations

import os
import time
import warnings

import cellxgene_census
import scanpy as sc

warnings.filterwarnings("ignore")

DATA_DIR = os.environ.get("SCFM_DATA_DIR", "/data/analysis/data_ganguly")
CENSUS_VERSION = "2025-11-08"


def download_tabula_sapiens_2() -> str:
    out_path = os.path.join(DATA_DIR, "tabula_sapiens_2_TF_scfm.h5ad")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 5e6:
        print(f"[skip] {out_path} already exists ({os.path.getsize(out_path)/1e6:.1f} MB)")
        return out_path
    t0 = time.time()
    print(f"[{time.time()-t0:.1f}s] Opening Census {CENSUS_VERSION}...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"[{time.time()-t0:.1f}s] Querying Tabula Sapiens 2.0 (donors TSP17-30)...")
        adata = cellxgene_census.get_anndata(
            census,
            organism="homo_sapiens",
            measurement_name="RNA",
            obs_value_filter=(
                "donor_id in ['TSP17', 'TSP18', 'TSP19', 'TSP20', 'TSP21', "
                "'TSP22', 'TSP23', 'TSP24', 'TSP25', 'TSP26', 'TSP27', "
                "'TSP28', 'TSP29', 'TSP30'] "
                "and is_primary_data == True "
                "and suspension_type == 'cell'"
            ),
            obs_column_names=[
                "cell_type", "tissue", "tissue_general", "donor_id",
                "assay", "dataset_id", "suspension_type",
            ],
            obs_embeddings=["tf-sapiens", "tf-exemplar-human", "scvi"],
        )
        print(f"[{time.time()-t0:.1f}s] Tabula Sapiens 2.0: {adata.shape}")
        print(f"  obsm: {[(k, adata.obsm[k].shape) for k in adata.obsm]}")
        print(f"  donors: {adata.obs['donor_id'].nunique()}")
        print(f"  cell types: {adata.obs['cell_type'].nunique()}")
        print(f"  tissues: {adata.obs['tissue'].nunique()}")
        # Rename for clarity
        adata.obsm["X_tf_sapiens"] = adata.obsm.pop("tf-sapiens")
        adata.obsm["X_tf_exemplar"] = adata.obsm.pop("tf-exemplar-human")
        adata.obsm["X_scvi_census"] = adata.obsm.pop("scvi")
        adata.uns["census_version"] = CENSUS_VERSION
        adata.uns["scfm_source"] = "CELLxGENE Census LTS 2025-11-08"
        adata.uns["scfm_atlas_name"] = "tabula_sapiens_2"
        adata.uns["scfm_tissue"] = "Tabula Sapiens 2.0 (multi-organ)"
        adata.uns["scfm_donors"] = "TSP17-TSP30 (held-out from TF pretraining)"
        print(f"[{time.time()-t0:.1f}s] Writing {out_path}...")
        adata.write_h5ad(out_path)
    print(f"  Done: {os.path.getsize(out_path)/1e6:.1f} MB")
    return out_path


def download_covid19() -> str:
    out_path = os.path.join(DATA_DIR, "covid19_TF_scfm.h5ad")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 5e6:
        print(f"[skip] {out_path} already exists ({os.path.getsize(out_path)/1e6:.1f} MB)")
        return out_path
    t0 = time.time()
    dataset_id = "01ad3cd7-3929-4654-84c0-6db05bd5fd59"  # 425k COVID-19 cells
    print(f"[{time.time()-t0:.1f}s] Opening Census {CENSUS_VERSION}...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"[{time.time()-t0:.1f}s] Querying COVID-19 dataset {dataset_id}...")
        adata = cellxgene_census.get_anndata(
            census,
            organism="homo_sapiens",
            measurement_name="RNA",
            obs_value_filter=(
                f"dataset_id == '{dataset_id}' "
                "and is_primary_data == True "
                "and suspension_type == 'cell'"
            ),
            obs_column_names=[
                "cell_type", "tissue", "tissue_general",
                "donor_id", "assay", "dataset_id",
                "disease", "suspension_type",
            ],
            obs_embeddings=["tf-sapiens", "tf-exemplar-human", "scvi"],
        )
        print(f"[{time.time()-t0:.1f}s] COVID-19: {adata.shape}")
        print(f"  obsm: {[(k, adata.obsm[k].shape) for k in adata.obsm]}")
        print(f"  cell types: {adata.obs['cell_type'].nunique()}")
        print(f"  donors: {adata.obs['donor_id'].nunique()}")
        adata.obsm["X_tf_sapiens"] = adata.obsm.pop("tf-sapiens")
        adata.obsm["X_tf_exemplar"] = adata.obsm.pop("tf-exemplar-human")
        adata.obsm["X_scvi_census"] = adata.obsm.pop("scvi")
        adata.uns["census_version"] = CENSUS_VERSION
        adata.uns["scfm_source"] = "CELLxGENE Census LTS 2025-11-08"
        adata.uns["scfm_atlas_name"] = "covid19"
        adata.uns["scfm_tissue"] = "lung (COVID-19 infected)"
        adata.uns["scfm_dataset_id"] = dataset_id
        print(f"[{time.time()-t0:.1f}s] Writing {out_path}...")
        adata.write_h5ad(out_path)
    print(f"  Done: {os.path.getsize(out_path)/1e6:.1f} MB")
    return out_path


def download_pbmc_10k() -> str:
    """Download a small PBMC dataset as a pbmc_10k analogue with TF embeddings."""
    out_path = os.path.join(DATA_DIR, "pbmc_10k_TF_scfm.h5ad")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 5e6:
        print(f"[skip] {out_path} already exists ({os.path.getsize(out_path)/1e6:.1f} MB)")
        return out_path
    t0 = time.time()
    print(f"[{time.time()-t0:.1f}s] Opening Census {CENSUS_VERSION}...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"[{time.time()-t0:.1f}s] Finding a PBMC 10k-sized dataset...")
        obs = cellxgene_census.get_obs(
            census,
            "homo_sapiens",
            value_filter=(
                "tissue == 'blood' "
                "and is_primary_data == True "
                "and suspension_type == 'cell' "
                "and assay == \"10x 3' v3\""
            ),
            column_names=["dataset_id", "cell_type"],
        )
        counts = obs.groupby("dataset_id").size().sort_values(ascending=False)
        # Find one close to 10k cells
        candidates = counts[(counts > 8000) & (counts < 15000)]
        if candidates.empty:
            candidates = counts[(counts > 5000) & (counts < 20000)]
        if candidates.empty:
            print("  [fail] no PBMC dataset in 5k-20k range")
            return None
        dataset_id = candidates.index[0]
        print(f"  Using dataset {dataset_id} ({counts[dataset_id]} cells)")
        print(f"[{time.time()-t0:.1f}s] Querying PBMC dataset...")
        adata = cellxgene_census.get_anndata(
            census,
            organism="homo_sapiens",
            measurement_name="RNA",
            obs_value_filter=(
                f"dataset_id == '{dataset_id}' "
                "and is_primary_data == True "
                "and suspension_type == 'cell'"
            ),
            obs_column_names=[
                "cell_type", "tissue", "tissue_general",
                "donor_id", "assay", "dataset_id",
                "suspension_type",
            ],
            obs_embeddings=["tf-sapiens", "tf-exemplar-human", "scvi"],
        )
        print(f"[{time.time()-t0:.1f}s] PBMC: {adata.shape}")
        print(f"  obsm: {[(k, adata.obsm[k].shape) for k in adata.obsm]}")
        adata.obsm["X_tf_sapiens"] = adata.obsm.pop("tf-sapiens")
        adata.obsm["X_tf_exemplar"] = adata.obsm.pop("tf-exemplar-human")
        adata.obsm["X_scvi_census"] = adata.obsm.pop("scvi")
        adata.uns["census_version"] = CENSUS_VERSION
        adata.uns["scfm_source"] = "CELLxGENE Census LTS 2025-11-08"
        adata.uns["scfm_atlas_name"] = "pbmc_10k"
        adata.uns["scfm_tissue"] = "PBMC (blood)"
        adata.uns["scfm_dataset_id"] = dataset_id
        print(f"[{time.time()-t0:.1f}s] Writing {out_path}...")
        adata.write_h5ad(out_path)
    print(f"  Done: {os.path.getsize(out_path)/1e6:.1f} MB")
    return out_path


if __name__ == "__main__":
    print("=" * 70)
    print("TF paper evaluation datasets from CELLxGENE Census")
    print("=" * 70)

    for fn, name in [
        (download_pbmc_10k, "pbmc_10k"),
        (download_covid19, "covid19"),
        (download_tabula_sapiens_2, "tabula_sapiens_2"),
    ]:
        print(f"\n{'='*70}")
        print(f"Dataset: {name}")
        print("=" * 70)
        try:
            fn()
        except Exception as exc:
            print(f"  [ERROR] {name}: {exc}")
            import traceback
            traceback.print_exc()
