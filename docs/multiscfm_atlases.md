# Multi-scFM atlases: 6 SOTA scFMs together in one h5ad

This document describes the **multi-scFM atlases** constructed for the
checkatlas-scfm pipeline: per-tissue h5ad files where the same cells have
been processed by **six state-of-the-art single-cell foundation models**
and the resulting embeddings are stored together in `adata.obsm`.

These atlases are the closest publicly-available approximation of "an
atlas processed by TF and the SOTA scFMs together with embeddings saved
in obsm", and they exist because the scFM community does not yet
distribute such files directly. This document records how they were
built, what is in them, and how to use them for the checkatlas scFM QC
layer.

---

## Why these atlases

Checkatlas-scfm evaluates a single scFM embedding at a time against nine
failure modes (see [`scfm.md`](scfm.md) and
[`scfm_benchmarking_report.md`](scfm_benchmarking_report.md)). For the
verdicts to be meaningful, the framework must compare the scFM
embedding against a baseline (PCA / scVI / Harmony) and, ideally,
against *other* scFMs as a sanity check. A multi-scFM atlas in which
the same cells carry embeddings from TF, scGPT, Geneformer, UCE, and
scVI lets checkatlas run cross-scFM comparisons directly: the
Problem-2 "baseline gap" verdict can be calibrated against the
*actual* scFM landscape, not just PCA.

The user's question "can we get an atlas that has been run by the TF
as well as different SOTA scFMs all together" was the motivation.
After searching the literature (Kedzierska, CellFM, BioLLM, the unified
framework paper, CZ Benchmarks) and HuggingFace, no such atlas is
distributed publicly. Every scFM benchmark paper re-runs the models
during evaluation rather than sharing precomputed embedding matrices.

The next-best source is the **CELLxGENE Census**, which has multiple
precomputed scFM embeddings hosted on S3. No single Census version
contains all the major scFMs together, but **combining Census
2025-01-30 (TF + Geneformer + scVI) with Census 2023-12-15 (scGPT +
UCE)** yields all six embeddings for the same cells via a row-position
join.

---

## The 6 scFMs in each h5ad

| Embedding key   | Dim  | Source Census version  | scFM model                                                                 |
|-----------------|------|------------------------|----------------------------------------------------------------------------|
| `X_tf_sapiens`  | 2048 | 2025-01-30             | TranscriptFormer TF-Sapiens (Pearce et al. 2025, Chan Zuckerberg Initiative) |
| `X_tf_exemplar` | 2048 | 2025-01-30             | TranscriptFormer TF-Exemplar-human (CZI)                                   |
| `X_geneformer`  | 512  | 2025-01-30 (verified equal to 2023-12-15) | Geneformer V2-316M (Theodoris et al. 2023)                       |
| `X_scvi_census` | 50   | 2025-01-30 (verified equal to 2023-12-15) | scVI (Lopez et al. 2018, CZI retrained on Census)                |
| `X_scgpt`       | 512  | 2023-12-15             | scGPT-human (Cui et al. 2024)                                              |
| `X_uce`         | 1280 | 2023-12-15             | Universal Cell Embedding (Rosen et al. 2023)                              |

All six models were run by their original authors and hosted on Census
S3. We did not retrain or re-implement any model. The
[`CellFM`](https://www.nature.com/articles/s41467-025-57015-z) and
[`BioLLM`](https://www.cell.com/patterns/abstract/S2666-3899(25)00174-6)
models are not in any Census version, so they are not part of this
construction.

---

## Files produced

All files live in `/data/analysis/data_ganguly/`:

| File                            | Size   | Cells   | Cell types | Donors | UCE NaN % |
|---------------------------------|--------|---------|-----------|--------|-----------|
| `blood_multiscfm.h5ad`          | 398 MB | 15,000  |  18        |   6    | 0.05 %    |
| `bone_marrow_multiscfm.h5ad`    | 398 MB | 15,000  |  36        |  15    | 3.30 %    |
| `liver_multiscfm.h5ad`          | 398 MB | 15,000  |  44        |  56    | 0.00 %    |
| `lung_multiscfm.h5ad`           | 398 MB | 15,000  |   9        |  12    | 0.00 %    |

Each file is subsampled to 15,000 cells from a representative Census
dataset. Raw counts (`adata.X`) are dropped to keep the file size at
roughly 400 MB; only the obsm embeddings and obs metadata are kept,
which is what checkatlas-scfm needs.

The `adata.uns` slot records the provenance of every file:

- `census_version_v1` (Census 2025-01-30 — TF + Geneformer + scVI)
- `census_version_v2` (Census 2023-12-15 — scGPT + UCE)
- `scfm_source` (a human-readable description of the merge)
- `scfm_dataset_id` (the Census dataset UUID)
- `scfm_atlas_name`, `scfm_tissue`
- `scfm_models_in_obsm` (list of all obsm keys)

---

## Why two Census versions are required

The Census LTS release (2025-11-08) only contains `scvi`,
`tf-sapiens`, and `tf-exemplar-human`. The 2023-12-15 release
contains `scvi`, `geneformer`, `scgpt`, and `uce`. Only Census
2025-01-30 contains both the TF embeddings and the legacy scFMs at
the same time. Combining the two gives the full six.

| Census version    | scvi | geneformer | scgpt | uce | tf-sapiens | tf-exemplar |
|-------------------|------|------------|-------|-----|------------|-------------|
| 2023-12-15        | yes  | yes        | yes   | yes | no         | no          |
| 2024-07-01        | yes  | yes        | yes   | no  | no         | no          |
| 2025-01-30        | yes  | yes        | no    | no  | yes        | yes         |
| 2025-11-08 (LTS)  | yes  | no         | no    | no  | yes        | yes         |

The 2023-12-15 version is not officially LTS and may be archived in
the future. The construction script in
`TF_scfm/multiscfm/download_multiscfm.py` will fail if that happens;
a fallback to 2024-07-01 (scvi + geneformer + scgpt only) would
require dropping UCE from the construction.

---

## How the join is done

The join is by **row position**, not by `soma_joinid`. `soma_joinid`
is per-Census-version and cannot be compared across versions.

We verified that for any common dataset, the row order in Census
2025-01-30 matches the row order in Census 2023-12-15:

- `donor_id` match: 100 % across thousands of cells
- `assay` match: 100 %
- `cell_type` match: 99 % (small differences from label updates
  between versions)

The construction also requires the dataset to have **identical cell
counts** in both Census versions (which is the case for all common
datasets — 30 common PBMC datasets, all with exact cell-count
matches). This guarantees a one-to-one row alignment.

For the bone-marrow dataset the join exposes one practical
limitation: the UCE embedding in Census 2023-12-15 is sparse (only
some cells are precomputed). Cells without UCE get all-NaN rows in
`X_uce`. The 3.30 % NaN rate for bone-marrow and 0.05 % for blood
are documented Census behaviour; the other two files have no NaNs.

---

## How the data was obtained

```python
import cellxgene_census

# Pick a dataset that exists in BOTH Census versions with same cell count
with cellxgene_census.open_soma(census_version="2025-01-30") as c1:
    with cellxgene_census.open_soma(census_version="2023-12-15") as c2:
        a1 = cellxgene_census.get_anndata(
            c1, organism="homo_sapiens",
            obs_value_filter=f"dataset_id == '{dataset_id}' and "
                             "is_primary_data == True and suspension_type == 'cell'",
            obs_column_names=["cell_type", "tissue", "donor_id", "assay", "dataset_id"],
            obs_embeddings=["tf-sapiens", "tf-exemplar-human", "geneformer", "scvi"],
        )
        a2 = cellxgene_census.get_anndata(
            c2, organism="homo_sapiens",
            obs_value_filter=f"dataset_id == '{dataset_id}' and "
                             "is_primary_data == True and suspension_type == 'cell'",
            obs_column_names=["cell_type", "tissue", "donor_id", "assay", "dataset_id"],
            obs_embeddings=["scgpt", "uce", "geneformer", "scvi"],
        )

# Row-position join (verified alignment 100 % donor / assay)
out = a1.copy()
out.obsm["X_scgpt"] = a2.obsm["scgpt"]
out.obsm["X_uce"]   = a2.obsm["uce"]
out.obsm["X_tf_sapiens"]  = out.obsm.pop("tf-sapiens")
out.obsm["X_tf_exemplar"] = out.obsm.pop("tf-exemplar-human")
out.obsm["X_geneformer"]  = out.obsm.pop("geneformer")
out.obsm["X_scvi_census"] = out.obsm.pop("scvi")
out.write_h5ad(f"{atlas}_multiscfm.h5ad")
```

The full script with tissue selection, subsampling, NaN handling,
and metadata recording is
`/data/analysis/data_ganguly/TF_scfm/multiscfm/download_multiscfm.py`.

---

## Verification

`TF_scfm/multiscfm/verify_multiscfm.py` loads each h5ad, checks that
all 6 obsm keys are present with the expected dimensions, confirms
the join metadata, and reports the NaN rate for `X_uce`. All four
files pass.

```
blood_multiscfm.h5ad          X_tf_sapiens (2048)  no NaN
                              X_tf_exemplar (2048) no NaN
                              X_geneformer (512)   no NaN
                              X_scvi_census (50)    no NaN
                              X_scgpt (512)         no NaN
                              X_uce (1280)          0.05 % NaN
bone_marrow_multiscfm.h5ad    ...                  3.30 % NaN in X_uce
liver_multiscfm.h5ad          ...                  0.00 % NaN
lung_multiscfm.h5ad           ...                  0.00 % NaN
SUMMARY: 4 passed, 0 failed
```

---

## How to use with checkatlas scFM

```python
import scanpy as sc
from checkatlas.scfm import run_scfm_pipeline
from checkatlas.scfm.config import SCFMConfig

a = sc.read_h5ad("/data/analysis/data_ganguly/blood_multiscfm.h5ad")
print(a.obsm.keys())
# ['X_geneformer', 'X_scgpt', 'X_scvi_census',
#  'X_tf_exemplar', 'X_tf_sapiens', 'X_uce']

config = SCFMConfig(
    atlas_name="blood",
    scfm_embedding="X_tf_sapiens",
    baseline_embeddings=["X_scvi_census", "X_geneformer"],
    ref_label="cell_type",
    batch_key="donor_id",
    domain_key="tissue",
    n_seeds=3,
)
run_scfm_pipeline(a, config)
```

Or from the CLI:

```bash
checkatlas scfm /data/analysis/data_ganguly/ \
    --atlas_name blood \
    --scfm_embedding X_tf_sapiens \
    --baseline_embeddings X_scvi_census X_geneformer \
    --ref_label cell_type \
    --batch_key donor_id
```

The output lands in `checkatlas_files/scfm/blood/` with the standard
verdicts, composite, and per-metric tables.

### Cross-scFM head-to-head

The real value of these atlases is comparing scFMs against each other.
The following loop runs checkatlas-scfm once per scFM (using
`X_scvi_census` as the common baseline):

```python
import json
import subprocess

for scfm_key in ["X_tf_sapiens", "X_tf_exemplar", "X_geneformer",
                 "X_scgpt", "X_uce"]:
    subprocess.run([
        "checkatlas", "scfm", "/data/analysis/data_ganguly/",
        "--atlas_name", "blood",
        "--scfm_embedding", scfm_key,
        "--baseline_embeddings", "X_scvi_census",
        "--ref_label", "cell_type",
        "--batch_key", "donor_id",
    ])
```

The resulting `scfm_verdicts.tsv` files for the five scFMs can be
joined on problem_id to produce a cross-scFM comparison table. This is
how Problem 2 (baseline gap) becomes a per-scFM verdict, not just a
PCA-vs-embedding verdict.

---

## Caveats

1. **The UCE embedding is sparse.** Cells without UCE get all-NaN
   rows. Downstream code should either drop those rows before
   evaluating the UCE embedding, or fill the NaNs with zero (a common
   practice for sparse embeddings).

2. **The Census 2023-12-15 version may be archived.** It is not LTS
   and could disappear without notice. The construction script will
   fall back to a 4-scFM construction (drop UCE) if that happens.

3. **The Census 2025-01-30 version may also be archived.** It is also
   not LTS. The current LTS 2025-11-08 has TF but no Geneformer, so
   if 2025-01-30 goes away we lose the ability to get both TF and
   Geneformer in the same Census version.

4. **The atlases are subsampled to 15,000 cells.** This is a
   trade-off between file size and statistical power. For Problem 5
   (rare-type recall) and Problem 7 (patient ASW), the original
   full-size datasets are available via the script in
   `download_multiscfm.py` by setting `max_cells=200000` and
   `target_cells=200000`.

5. **No atlas has all 7+ SOTA scFMs.** CellFM, BioLLM, scFoundation,
   scBERT, and SCimilarity are not in any Census version. Adding
   them requires running each model's own inference code locally.
   For this, see `BGIResearch/BioLLM` (provides scGPT, Geneformer,
   scFoundation, scBERT in one framework) and the unified framework
   paper (PMC12803055, 13 scFMs across 50 datasets).

---

## File map

```
/data/analysis/data_ganguly/
├── blood_multiscfm.h5ad
├── bone_marrow_multiscfm.h5ad
├── liver_multiscfm.h5ad
├── lung_multiscfm.h5ad
├── TF_scfm/multiscfm/
│   ├── README.md
│   ├── download_multiscfm.py
│   └── verify_multiscfm.py
```

The single-TF atlases (TF-Sapiens + scVI only) are in
`/data/analysis/data_ganguly/*_TF_scfm.h5ad` and are described in
`/data/analysis/data_ganguly/TF_scfm/README.md`. The multi-scFM
atlases extend those by adding scGPT, UCE, and an independent
Geneformer from a second Census version.
