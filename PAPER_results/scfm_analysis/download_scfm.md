# Reproduce scFM-processed atlases

This document explains how to reproduce all scFM-processed single-cell atlases
used in the checkatlas-scfm evaluation. Every file was obtained by querying
the **CELLxGENE Census API** — no model weights were downloaded and no GPUs
are required. All embeddings are the original authors' published outputs hosted
on Census S3.

## Prerequisites

```bash
pip install cellxgene-census tiledbsoma scanpy
# Minimum version: cellxgene-census >= 1.9.1, Python >= 3.9
```

Census data is hosted on AWS S3 (`us-west-2`) with anonymous read access.
Internet connectivity is required; no AWS credentials are needed.

## Quick start

To reproduce all 13 atlases, run these scripts in order:

```bash
# 1. Single-TF atlases (9 files, ~48 GB total)
python download_census_tf.py
python download_tf_eval_datasets.py
python download_pancreas_myeloid.py

# 2. Multi-scFM atlases (4 files, ~1.6 GB total)
python download_multiscfm.py

# 3. Verify
python verify_all.py
python verify_multiscfm.py
```

All scripts are in this directory. Output files land in the directory
specified by `DATA_DIR` inside each script (default: `/data/analysis/data_ganguly`).
Change it if running on a different machine.

---

## Tier 1 — Single-TF atlases (`*_TF_scfm.h5ad`)

9 atlases, each with **TF-Sapiens (2048-dim) + scVI (50-dim)** embeddings.
5 of 9 also have **TF-Exemplar (2048-dim)**. All sourced from **Census LTS 2025-11-08**.

| Script | Creates | Embeddings |
|--------|---------|------------|
| `download_census_tf.py` | `blood_TF_scfm.h5ad`, `bone_marrow_TF_scfm.h5ad`, `liver_TF_scfm.h5ad`, `lung_TF_scfm.h5ad` | X_tf_sapiens + X_scvi_census |
| `download_tf_eval_datasets.py` | `pbmc_10k_TF_scfm.h5ad`, `covid19_TF_scfm.h5ad`, `tabula_sapiens_2_TF_scfm.h5ad` | X_tf_sapiens + X_tf_exemplar + X_scvi_census |
| `download_pancreas_myeloid.py` | `hpancreas_TF_scfm.h5ad`, `myeloid_TF_scfm.h5ad` | X_tf_sapiens + X_tf_exemplar + X_scvi_census |

### How they were obtained

Each script calls `cellxgene_census.get_anndata` against **Census LTS 2025-11-08**
with `obs_embeddings=['tf-sapiens', 'scvi']` (and `tf-exemplar-human` for the
newer downloads). The function returns an AnnData with the embedding arrays
already in `obsm`. No model inference is run.

The dataset selection function picks a representative primary-data dataset
per tissue (5,000-80,000 cells range, `suspension_type='cell'`). The Census
dataset UUIDs are:

| Atlas | Census dataset ID |
|-------|-------------------|
| `blood_TF_scfm.h5ad` | `19e46756-9100-4e01-8b0e-23b557558a4c` |
| `bone_marrow_TF_scfm.h5ad` | `e22482ee-19e8-40bc-9f6e-541dc3c82c20` |
| `hpancreas_TF_scfm.h5ad` | `f7c1c579-2dc0-47e2-ba19-8165c5a0e353` |
| `liver_TF_scfm.h5ad` | `e3ed2ba4-edf5-40ac-8750-8a417ad1eefe` |
| `lung_TF_scfm.h5ad` | `3dc61ca1-ce40-46b6-8337-f27260fd9a03` |
| `myeloid_TF_scfm.h5ad` | `456e8b9b-f872-488b-871d-94534090a865` |
| `pbmc_10k_TF_scfm.h5ad` | `db59611b-42de-4035-93aa-1ed39f38b467` |
| `covid19_TF_scfm.h5ad` | `01ad3cd7-3929-4654-84c0-6db05bd5fd59` |
| `tabula_sapiens_2_TF_scfm.h5ad` | Tabula Sapiens 2.0 donors TSP17-30 (multi-dataset) |

---

## Tier 2 — Multi-scFM atlases (`*_multiscfm.h5ad`)

4 atlases, each with **6 SOTA scFM embeddings** for the same cells.
Constructed by merging two Census versions because no single version
contains all 6 scFMs together.

| Script | Creates | Embeddings |
|--------|---------|------------|
| `download_multiscfm.py` | `blood_multiscfm.h5ad`, `bone_marrow_multiscfm.h5ad`, `liver_multiscfm.h5ad`, `lung_multiscfm.h5ad` | X_tf_sapiens (2048d), X_tf_exemplar (2048d), X_geneformer (512d), X_scvi_census (50d), X_scgpt (512d), X_uce (1280d) |

### How they were obtained

**Two Census versions are needed** because:

| Embedding | Available in |
|-----------|-------------|
| TF-Sapiens, TF-Exemplar, Geneformer, scVI | Census **2025-01-30** |
| scGPT, UCE | Census **2023-12-15** |

The construction (in `download_multiscfm.py`) proceeds as follows:

1. **Find a dataset** that exists in BOTH Census versions with an identical
   cell count (e.g., `c2a461b1-...` for blood). This guarantees a one-to-one
   row correspondence.

2. **Query Census 2025-01-30** for the dataset with
   `obs_embeddings=['tf-sapiens', 'tf-exemplar-human', 'geneformer', 'scvi']`.

3. **Query Census 2023-12-15** for the *same* dataset_id with
   `obs_embeddings=['scgpt', 'uce']`.

4. **Verify row alignment** by three independent checks:
   - `obs_names` (Census internal soma_joinid indices) match 100 %
   - `donor_id` matches 100 %
   - `assay` matches 100 %
   - Both queries return the identical cell count.
   
   If any check fails, the script aborts.

5. **Merge** the obsm arrays from the v2 AnnData into the v1 AnnData by row
   position. This is valid because the verification confirms the two
   AnnData objects have cells in the identical row order.

6. **Subsample** to 15,000 cells and drop raw counts to keep file size at
   roughly 400 MB per atlas.

The merged dataset UUIDs and their original Census cell counts:

| Atlas | Census dataset ID | Source cells | Subsampled to |
|-------|-------------------|-------------|---------------|
| `blood_multiscfm.h5ad` | `c2a461b1-0c15-4047-9fcb-1f966fe55100` | 97,499 | 15,000 |
| `bone_marrow_multiscfm.h5ad` | `965386e9-1e4f-466d-bf59-ebdca4b66b9b` | 92,676 | 15,000 |
| `liver_multiscfm.h5ad` | `2adb1f8a-a6b1-4909-8ee8-484814e2d4bf` | 598,266 | 15,000 |
| `lung_multiscfm.h5ad` | `576f193c-75d0-4a11-bd25-8676587e6dc2` | 147,137 | 15,000 |

---

## Tier 3 — pancreas_scib (original benchmark)

The file `pancreas_scib.h5ad` (16,382 cells) is the original scIB pancreas
benchmark from Luecken et al. 2022, obtained via the Kedzierska et al. 2024
figshare deposit. It is **NOT** processed via the CELLxGENE Census.

**Source:** `https://doi.org/10.6084/m9.figshare.12420968.v8`  
**Download:** `https://figshare.com/ndownloader/files/24539828`  
**Citation:** Luecken, M. et al. Benchmarking atlas-level data integration in
single-cell genomics. *Nature Methods* 19, 41-50 (2022).  
**Additional processing:** The `X_geneformer`, `X_scgpt_human`, and
`X_scvi_dropin` embeddings in `obsm` are **synthetic** embeddings generated by
the checkatlas scFM pipeline for QC testing; they are not real scFM outputs.
Only `X_pca` reflects computations on the real data.

---

## Census version availability

Not every Census version has every embedding. If a version is retired from
the Census directory, the corresponding script will fail. The availability
matrix as of July 2026:

| Census version    | scvi | geneformer | scgpt | uce | tf-sapiens | tf-exemplar |
|-------------------|------|------------|-------|-----|------------|-------------|
| 2023-12-15        | yes  | yes        | yes   | yes | no         | no          |
| 2024-07-01        | yes  | yes        | yes   | no  | no         | no          |
| 2025-01-30        | yes  | yes        | no    | no  | yes        | yes         |
| 2025-11-08 (LTS)  | yes  | no         | no    | no  | yes        | yes         |

**If 2023-12-15 or 2025-01-30 are retired**, the multi-scFM construction will
lose scGPT + UCE (2023-12-15) or Geneformer (2025-01-30). The fallback is:
use the 4 scFMs available in the remaining version. To get all 6, the
Census LTS would need to be updated to host all embeddings in a single
release.

---

## Verification

After downloading, run the verification scripts to confirm each file has the
expected embeddings with the correct dimensions and no missing data:

```bash
python verify_all.py      # checks all 9 single-TF atlases
python verify_multiscfm.py # checks all 4 multi-scFM atlases
```

Both scripts load each h5ad, check `obsm` keys, verify dimensions match
the expected values (TF=2048, Geneformer/scGPT=512, UCE=1280, scVI=50),
report NaN rates, and confirm the join verification metadata is present.
A passing run shows `[PASS]` for each file.
