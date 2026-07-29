"""Download single-cell atlases from CELLxGENE Census with TF-Sapiens embeddings.

For each tissue of interest, we:
  1. Find a representative Census dataset (CZI-curated, primary data, human)
  2. Pull the raw counts + tf-sapiens + scvi embeddings + cell type metadata
  3. Save to /data/analysis/data_ganguly/<atlas_name>_TF_scfm.h5ad
  4. Verify obsm keys are populated

Uses Census LTS 2025-11-08, which contains:
  - tf-sapiens        (CZI-CxG-czi-17, 159M cells, 2048-dim)
  - tf-exemplar-human (CZI-CxG-czi-18, 159M cells, 2048-dim)
  - tf-exemplar-mouse (CZI-CxG-czi-19,  44M cells, 2048-dim)
  - scvi              (CZI-CxG-czi-15, 159M cells,   50-dim, human)
  - scvi              (CZI-CxG-czi-16,  44M cells,   50-dim, mouse)
"""

from __future__ import annotations

import os
import time
import warnings

import cellxgene_census
from cellxgene_census.experimental import get_all_available_embeddings

warnings.filterwarnings("ignore")

DATA_DIR = "/data/analysis/data_ganguly"
TF_SFM_DIR = os.path.join(DATA_DIR, "TF_scfm")
os.makedirs(TF_SFM_DIR, exist_ok=True)

CENSUS_VERSION = "2025-11-08"
EMBEDDINGS = ["tf-sapiens", "scvi"]


def pick_dataset(census, tissue: str, max_cells: int = 80000) -> str | None:
    """Find a representative primary-data dataset for the given tissue, bounded by cell count."""
    obs = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter=(
            f"tissue == '{tissue}' "
            "and is_primary_data == True "
            "and suspension_type == 'cell'"
        ),
        column_names=["dataset_id", "cell_type", "assay"],
    )
    counts = obs.groupby("dataset_id").size().sort_values(ascending=False)
    counts = counts[(counts > 5000) & (counts < max_cells)]
    if counts.empty:
        return None
    return counts.index[0]


def download_one(tissue: str, atlas_label: str, max_cells: int = 80000) -> str | None:
    """Download one tissue atlas with TF embeddings; return output path."""
    out_path = os.path.join(DATA_DIR, f"{atlas_label}_TF_scfm.h5ad")
    if os.path.exists(out_path):
        size_mb = os.path.getsize(out_path) / 1e6
        if size_mb > 5:  # at least 5 MB
            print(f"  [skip] {out_path} already exists ({size_mb:.1f} MB)")
            return out_path

    t0 = time.time()
    print(f"[{time.time()-t0:.1f}s] Opening Census {CENSUS_VERSION}...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"[{time.time()-t0:.1f}s] Finding representative dataset for tissue='{tissue}'...")
        dataset_id = pick_dataset(census, tissue, max_cells=max_cells)
        if dataset_id is None:
            print(f"  [fail] no suitable dataset for tissue='{tissue}'")
            return None
        print(f"[{time.time()-t0:.1f}s] Using dataset {dataset_id}")
        print(f"[{time.time()-t0:.1f}s] Querying AnnData with embeddings={EMBEDDINGS}...")

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
                "suspension_type", "is_primary_data",
            ],
            obs_embeddings=EMBEDDINGS,
        )
        print(f"[{time.time()-t0:.1f}s] AnnData shape: {adata.shape}")
        print(f"  obsm keys: {list(adata.obsm.keys())}")
        for k in adata.obsm:
            print(f"    {k}: {adata.obsm[k].shape}, dtype={adata.obsm[k].dtype}")
        print(f"  obs columns: {list(adata.obs.columns)}")
        print(f"  cell_type n_unique: {adata.obs['cell_type'].nunique()}")
        print(f"  donor_id n_unique: {adata.obs['donor_id'].nunique()}")

        # Rename obsm keys to disambiguate from the user's existing baselines
        adata.obsm["X_tf_sapiens"] = adata.obsm.pop("tf-sapiens")
        adata.obsm["X_scvi_census"] = adata.obsm.pop("scvi")
        adata.uns["census_version"] = CENSUS_VERSION
        adata.uns["scfm_source"] = "CELLxGENE Census LTS 2025-11-08"
        adata.uns["scfm_atlas_name"] = atlas_label
        adata.uns["scfm_tissue"] = tissue
        adata.uns["scfm_dataset_id"] = dataset_id

        # If the user already has a same-tissue file, preserve its obsm baselines
        existing_path = os.path.join(DATA_DIR, f"{atlas_label}.h5ad")
        if os.path.exists(existing_path) and atlas_label in {
            "blood", "bone_marrow", "liver", "lung", "pbmc_10k",
        }:
            print(f"[{time.time()-t0:.1f}s] Merging with existing {atlas_label}.h5ad baselines...")
            import anndata as ad
            base = sc.read_h5ad(existing_path)
            print(f"  base shape: {base.shape}, obsm: {list(base.obsm.keys())}")
            # Try to align by donor_id + cell_type, but only if shapes match
            # Simpler: keep TF embeddings only for the Census slice, write as a separate file
            print(f"  (baselines from existing {atlas_label}.h5ad will not be merged - different datasets)")

        print(f"[{time.time()-t0:.1f}s] Writing {out_path}...")
        adata.write_h5ad(out_path)
        size_mb = os.path.getsize(out_path) / 1e6
        print(f"[{time.time()-t0:.1f}s] Wrote {out_path} ({size_mb:.1f} MB)")
        return out_path


if __name__ == "__main__":
    import scanpy as sc  # noqa: E402  (imported here to keep top-level light)

    tissues = [
        ("blood", "blood", 80000),
        ("bone marrow", "bone_marrow", 80000),
        ("liver", "liver", 80000),
        ("lung", "lung", 80000),
        # Skip hpancreas (pancreas), pbmc_10k, myeloid - these are external datasets
        # not directly available in Census
    ]

    for tissue, label, max_cells in tissues:
        print(f"\n{'='*70}")
        print(f"Tissue: {tissue} -> {label}_TF_scfm.h5ad (max {max_cells} cells)")
        print(f"{'='*70}")
        try:
            download_one(tissue, label, max_cells=max_cells)
        except Exception as exc:
            print(f"  [ERROR] {label}: {exc}")
            import traceback
            traceback.print_exc()
