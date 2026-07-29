"""Final verification: load each *_TF_scfm.h5ad and confirm TF-Sapiens embeddings are present."""

from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")

DATA_DIR = os.environ.get("SCFM_DATA_DIR", "/data/analysis/data_ganguly")

files = [
    "blood_TF_scfm.h5ad",
    "bone_marrow_TF_scfm.h5ad",
    "covid19_TF_scfm.h5ad",
    "hpancreas_TF_scfm.h5ad",
    "liver_TF_scfm.h5ad",
    "lung_TF_scfm.h5ad",
    "myeloid_TF_scfm.h5ad",
    "pbmc_10k_TF_scfm.h5ad",
    "tabula_sapiens_2_TF_scfm.h5ad",
]


def verify(path: str) -> bool:
    name = os.path.basename(path)
    print(f"\n{'='*70}\n{name}\n{'='*70}")
    if not os.path.exists(path):
        print(f"  [FAIL] file does not exist")
        return False
    size_mb = os.path.getsize(path) / 1e6
    print(f"  size: {size_mb:.1f} MB")
    try:
        a = sc.read_h5ad(path)
    except Exception as exc:
        print(f"  [FAIL] cannot read: {exc}")
        return False
    print(f"  shape: {a.shape}")
    print(f"  obsm: {list(a.obsm.keys())}")

    has_tf_sapiens = "X_tf_sapiens" in a.obsm
    print(f"  has X_tf_sapiens: {has_tf_sapiens}")
    if has_tf_sapiens:
        emb = a.obsm["X_tf_sapiens"]
        print(f"    shape: {emb.shape}, dtype: {emb.dtype}")
        print(f"    mean: {emb.mean():.4f}, std: {emb.std():.4f}")
        print(f"    min: {emb.min():.4f}, max: {emb.max():.4f}")
        print(f"    nan_count: {np.isnan(emb).sum()}")
        if emb.shape[1] != 2048:
            print(f"    [WARN] expected 2048 dims, got {emb.shape[1]}")
            return False
    else:
        return False

    has_scvi = "X_scvi_census" in a.obsm
    print(f"  has X_scvi_census: {has_scvi}")
    if has_scvi:
        emb_scvi = a.obsm["X_scvi_census"]
        print(f"    shape: {emb_scvi.shape}, mean: {emb_scvi.mean():.4f}, std: {emb_scvi.std():.4f}")

    has_tf_exemplar = "X_tf_exemplar" in a.obsm
    print(f"  has X_tf_exemplar: {has_tf_exemplar}")
    if has_tf_exemplar:
        emb_te = a.obsm["X_tf_exemplar"]
        print(f"    shape: {emb_te.shape}, mean: {emb_te.mean():.4f}, std: {emb_te.std():.4f}")

    print(f"  uns: {list(a.uns.keys())}")
    if "scfm_source" in a.uns:
        print(f"    scfm_source: {a.uns['scfm_source']}")
    if "census_version" in a.uns:
        print(f"    census_version: {a.uns['census_version']}")
    if "scfm_atlas_name" in a.uns:
        print(f"    scfm_atlas_name: {a.uns['scfm_atlas_name']}")
    if "scfm_tissue" in a.uns:
        print(f"    scfm_tissue: {a.uns['scfm_tissue']}")
    if "scfm_dataset_id" in a.uns:
        print(f"    scfm_dataset_id: {a.uns['scfm_dataset_id']}")

    if "cell_type" in a.obs:
        print(f"  cell_type: {a.obs['cell_type'].nunique()} unique, top: {dict(a.obs['cell_type'].value_counts().head(3))}")
    if "donor_id" in a.obs:
        print(f"  donor_id: {a.obs['donor_id'].nunique()} unique")
    return True


def main():
    print("=" * 70)
    print("Verification: TF-Sapiens embeddings in CELLxGENE Census-derived h5ad")
    print("=" * 70)
    n_pass = 0
    n_fail = 0
    for f in files:
        path = os.path.join(DATA_DIR, f)
        ok = verify(path)
        if ok:
            n_pass += 1
        else:
            n_fail += 1
    print(f"\n{'='*70}\nSUMMARY: {n_pass} passed, {n_fail} failed\n{'='*70}")
    return n_fail


if __name__ == "__main__":
    sys.exit(main())
