"""Download pancreas + myeloid atlases from CELLxGENE Census with TF embeddings."""

from __future__ import annotations

import os
import time
import warnings

import cellxgene_census
import scanpy as sc

warnings.filterwarnings("ignore")

DATA_DIR = "/data/analysis/data_ganguly"
CENSUS_VERSION = "2025-11-08"


def download_one(dataset_id: str, atlas_label: str, tissue: str, max_cells: int = 80000):
    out_path = os.path.join(DATA_DIR, f"{atlas_label}_TF_scfm.h5ad")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 5e6:
        print(f"[skip] {out_path} already exists ({os.path.getsize(out_path)/1e6:.1f} MB)")
        return out_path
    t0 = time.time()
    print(f"[{time.time()-t0:.1f}s] Opening Census {CENSUS_VERSION}...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"[{time.time()-t0:.1f}s] Querying {atlas_label} (dataset {dataset_id}, max {max_cells} cells)...")
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
        # Subsample to max_cells if too large
        if adata.n_obs > max_cells:
            print(f"[{time.time()-t0:.1f}s] Subsampling {adata.n_obs} -> {max_cells}")
            sc.pp.subsample(adata, n_obs=max_cells, random_state=42)
        print(f"[{time.time()-t0:.1f}s] {atlas_label}: {adata.shape}")
        print(f"  obsm: {[(k, adata.obsm[k].shape) for k in adata.obsm]}")
        adata.obsm["X_tf_sapiens"] = adata.obsm.pop("tf-sapiens")
        adata.obsm["X_tf_exemplar"] = adata.obsm.pop("tf-exemplar-human")
        adata.obsm["X_scvi_census"] = adata.obsm.pop("scvi")
        adata.uns["census_version"] = CENSUS_VERSION
        adata.uns["scfm_source"] = "CELLxGENE Census LTS 2025-11-08"
        adata.uns["scfm_atlas_name"] = atlas_label
        adata.uns["scfm_tissue"] = tissue
        adata.uns["scfm_dataset_id"] = dataset_id
        print(f"[{time.time()-t0:.1f}s] Writing {out_path}...")
        adata.write_h5ad(out_path)
    print(f"  Done: {os.path.getsize(out_path)/1e6:.1f} MB")
    return out_path


if __name__ == "__main__":
    print("=" * 70)
    print("Pancreas + Myeloid atlases with TF embeddings")
    print("=" * 70)

    # Pancreas - hPancreas analogue (smaller dataset ~25k cells)
    print("\n" + "=" * 70)
    print("Pancreas (hPancreas analogue)")
    print("=" * 70)
    try:
        download_one(
            dataset_id="2adb1f8a-a6b1-4909-8ee8-484814e2d4bf",
            atlas_label="hpancreas",
            tissue="pancreas",
            max_cells=15000,
        )
    except Exception as exc:
        print(f"  [ERROR] hpancreas: {exc}")
        import traceback; traceback.print_exc()

    # Myeloid - smaller CD14+ monocyte from blood
    print("\n" + "=" * 70)
    print("Myeloid (CD14+ monocyte blood dataset)")
    print("=" * 70)
    try:
        download_one(
            dataset_id="456e8b9b-f872-488b-871d-94534090a865",
            atlas_label="myeloid",
            tissue="blood (CD14+ monocytes)",
            max_cells=10000,
        )
    except Exception as exc:
        print(f"  [ERROR] myeloid: {exc}")
        import traceback; traceback.print_exc()

    # Myeloid - CD14+ monocyte from blood, large
    print("\n" + "=" * 70)
    print("Myeloid (CD14+ monocyte blood dataset)")
    print("=" * 70)
    try:
        download_one(
            dataset_id="c838aec3-03ef-4398-b882-0e3912abfff0",
            atlas_label="myeloid",
            tissue="blood (CD14+ monocytes)",
            max_cells=20000,
        )
    except Exception as exc:
        print(f"  [ERROR] myeloid: {exc}")
        import traceback; traceback.print_exc()
