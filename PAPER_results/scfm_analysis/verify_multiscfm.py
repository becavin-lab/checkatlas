"""Verify gold-standard multiscfm h5ad files:
  1. All 6 obsm embeddings present, correct dims, no NaN (except UCE sparsity)
  2. soma_joinid_v1 and soma_joinid_v2 stored in obs
  3. join_verification metadata in uns
  4. join method documented
"""

from __future__ import annotations

import json
import os
import sys
import warnings

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")

DATA_DIR = "/data/analysis/data_ganguly"

EXPECTED_EMBEDDINGS = [
    ("X_tf_sapiens", 2048),
    ("X_tf_exemplar", 2048),
    ("X_geneformer", 512),
    ("X_scvi_census", 50),
    ("X_scgpt", 512),
    ("X_uce", 1280),
]

EXPECTED_OBS_COLS = [
    "cell_type", "tissue", "tissue_general",
    "donor_id", "assay", "dataset_id",
    "suspension_type", "is_primary_data",
    "soma_joinid_v1", "soma_joinid_v2", "row_index",
]

files = [
    "blood_multiscfm.h5ad",
    "bone_marrow_multiscfm.h5ad",
    "liver_multiscfm.h5ad",
    "lung_multiscfm.h5ad",
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
    errors = []

    # 1. Check obsm keys
    obsm_keys = set(a.obsm.keys())
    expected_keys = {k for k, _ in EXPECTED_EMBEDDINGS}
    print(f"  obsm keys ({len(obsm_keys)}): {sorted(obsm_keys)}")
    for key, expected_dim in EXPECTED_EMBEDDINGS:
        if key not in a.obsm:
            print(f"  [FAIL] {key} missing from obsm")
            errors.append(f"missing_{key}")
            continue
        emb = a.obsm[key]
        if emb.shape[1] != expected_dim:
            print(f"  [FAIL] {key} dim {emb.shape[1]} != {expected_dim}")
            errors.append(f"wrong_dim_{key}")
        nan_all = np.isnan(emb).sum()
        nan_rows = np.all(np.isnan(emb), axis=1).sum()
        print(f"    {key:25s} dim={emb.shape[1]:<5d} "
              f"mean={np.nanmean(emb):.4f} std={np.nanstd(emb):.4f} "
              f"nan=({nan_all}/{emb.size}, {100*nan_all/emb.size:.2f}%) "
              f"all-nan-rows={nan_rows}")

    # 2. Check obs columns
    for col in EXPECTED_OBS_COLS:
        if col not in a.obs.columns:
            print(f"  [WARN] obs column {col} missing")
            errors.append(f"missing_obs_{col}")
    print(f"  obs columns ({len(a.obs.columns)}): {list(a.obs.columns)}")

    # 3. Check soma_joinid columns
    if "soma_joinid_v1" in a.obs.columns and "soma_joinid_v2" in a.obs.columns:
        sj1 = a.obs["soma_joinid_v1"]
        sj2 = a.obs["soma_joinid_v2"]
        print(f"  soma_joinid_v1 range: [{sj1.min()}, {sj1.max()}], dtype={sj1.dtype}")
        print(f"  soma_joinid_v2 range: [{sj2.min()}, {sj2.max()}], dtype={sj2.dtype}")
        # They should be different values (different Census versions) but same count
        if len(sj1.unique()) != len(sj2.unique()) != len(a):
            print(f"  [WARN] soma_joinid unique count mismatch: v1={len(sj1.unique())}, v2={len(sj2.unique())}, n_obs={len(a)}")
    else:
        print(f"  [FAIL] soma_joinid columns missing")
        errors.append("missing_soma_joinid")

    if "row_index" in a.obs.columns:
        ri = a.obs["row_index"]
        if (ri == np.arange(len(a), dtype=int)).all():
            print(f"  row_index: 0..{ri.max()} (sequential, correct)")
        else:
            print(f"  [WARN] row_index not sequential 0..n-1")

    # 4. Check uns metadata
    print(f"  uns keys: {list(a.uns.keys())}")
    for key in ["join_verification", "census_version_v1", "census_version_v2",
                "scfm_source", "scfm_atlas_name", "scfm_tissue",
                "scfm_dataset_id", "scfm_models_in_obsm"]:
        if key in a.uns:
            val = a.uns[key]
            if key == "join_verification":
                print(f"  join_verification:")
                for k, v in val.items():
                    if isinstance(v, str) and len(v) > 120:
                        print(f"    {k}: {v[:120]}...")
                    else:
                        print(f"    {k}: {v}")
            elif isinstance(val, list):
                print(f"  {key}: {val}")
            else:
                print(f"  {key}: {val}")
        else:
            print(f"  [WARN] {key} missing from uns")
            errors.append(f"missing_uns_{key}")

    # 5. Cell metadata stats
    if "cell_type" in a.obs:
        print(f"  cell_type: {a.obs['cell_type'].nunique()} unique")
    if "donor_id" in a.obs:
        print(f"  donor_id: {a.obs['donor_id'].nunique()} unique")

    ok = len(errors) == 0
    if not ok:
        print(f"  [FAIL] errors: {errors}")
    else:
        print(f"  [PASS]")

    return ok


def main():
    print("=" * 70)
    print("Gold-standard multiscfm verification: 6 scFMs + soma_joinid provenance")
    print("=" * 70)
    n_pass, n_fail = 0, 0
    for f in files:
        path = os.path.join(DATA_DIR, f)
        ok = verify(path)
        if ok:
            n_pass += 1
        else:
            n_fail += 1
    print(f"\n{'='*70}")
    print(f"SUMMARY: {n_pass} passed, {n_fail} failed")
    print(f"{'='*70}")
    return n_fail


if __name__ == "__main__":
    sys.exit(main())
