# Multi-scFM atlases: 6 SOTA scFMs in one h5ad (gold-standard)

These files are the **closest publicly-available approximation** of "an atlas
processed by all SOTA scFMs together with embeddings saved in obsm." Each file
contains the same cells processed by 6 distinct foundation models, joined by
a verified row-position join with three independent alignment checks.

## Scientific validity of the join

**Method:** Row-position join with three independent verification checks.

**Why this is valid:**
- CELLxGENE Census stores data in a **consistent row order per dataset_id**
- For the same `dataset_id` with identical filters (`is_primary_data`,
  `suspension_type`), both Census versions return cells in identical
  soma_joinid-indexed row order
- We empirically verified this for all 4 tissues: **100 % donor_id match,
  100 % assay match, 100 % obs_names match, 100 % cell-count match**

**Why soma_joinid cannot be used directly as a join key:**
- soma_joinid changes between Census versions (different internal indexing)
- observation_joinid (actual cell barcode) is NOT available in Census
  2023-12-15
- Both Census queries for the same dataset_id with identical filters
  return cells in the same order, so the row-position join is deterministic
  and reproducible

**Verification stored in each h5ad:**
- `adata.uns['join_verification']` — detailed alignment stats (method,
  n_cells, obs_names_match_pct, donor_id_match_pct, assay_match_pct)
- `adata.obs['soma_joinid_v1']` — soma_joinid from Census 2025-01-30
- `adata.obs['soma_joinid_v2']` — soma_joinid from Census 2023-12-15
- `adata.obs['row_index']` — sequential 0..n-1 row index for join
  re-verification

## The 6 scFMs in each h5ad

| Embedding key   | Dim  | Source Census version  | scFM model |
|-----------------|------|------------------------|------------|
| `X_tf_sapiens`  | 2048 | 2025-01-30             | TranscriptFormer TF-Sapiens (Pearce et al. 2025, CZI) |
| `X_tf_exemplar` | 2048 | 2025-01-30             | TranscriptFormer TF-Exemplar-human (CZI) |
| `X_geneformer`  | 512  | 2025-01-30             | Geneformer V2-316M (Theodoris et al. 2023) |
| `X_scvi_census` | 50   | 2025-01-30             | scVI (Lopez et al. 2018, CZI retrained on Census) |
| `X_scgpt`       | 512  | 2023-12-15             | scGPT-human (Cui et al. 2024) |
| `X_uce`         | 1280 | 2023-12-15             | Universal Cell Embedding (Rosen et al. 2023) |

## Files

| File | Size | Cells | Cell types | Donors |
|------|------|-------|-----------|--------|
| `blood_multiscfm.h5ad`        | 398 MB | 15,000 |  26 |  22 |
| `bone_marrow_multiscfm.h5ad`  | 398 MB | 15,000 |  17 |  22 |
| `liver_multiscfm.h5ad`        | 398 MB | 15,000 |  44 |  56 |
| `lung_multiscfm.h5ad`         | 398 MB | 15,000 |  11 |  42 |

Each file is subsampled to 15,000 cells. Raw counts are dropped (only
obsm + obs + provenance metadata are kept).

## Why two Census versions

No single Census version has all 6 scFMs:

| Census version    | scvi | geneformer | scgpt | uce | tf-sapiens | tf-exemplar |
|-------------------|------|------------|-------|-----|------------|-------------|
| 2023-12-15        | yes  | yes        | yes   | yes | no         | no          |
| 2024-07-01        | yes  | yes        | yes   | no  | no         | no          |
| 2025-01-30        | yes  | yes        | no    | no  | yes        | yes         |
| 2025-11-08 (LTS)  | yes  | no         | no    | no  | yes        | yes         |

## UCE NaN values

The UCE embedding (Census 2023-12-15) is sparse. Only 2 out of 15,000
cells in lung have all-NaN UCE rows (0.01 %). Blood, bone_marrow, and
liver have zero NaN rows. For checkatlas-scfm, either zero-fill or
drop the NaN rows before evaluating UCE.

## How to use with checkatlas

```python
import scanpy as sc
a = sc.read_h5ad("/data/analysis/data_ganguly/blood_multiscfm.h5ad")
print(a.obsm.keys())
# ['X_geneformer', 'X_scgpt', 'X_scvi_census',
#  'X_tf_exemplar', 'X_tf_sapiens', 'X_uce']

# Cross-scFM comparison
for key in a.obsm:
    sc.pp.neighbors(a, use_rep=key)
    sc.tl.umap(a)

# CLI usage:
# checkatlas scfm /data/analysis/data_ganguly/ \
#     --atlas_name blood \
#     --scfm_embedding X_tf_sapiens \
#     --baseline_embeddings X_scvi_census X_geneformer \
#     --ref_label cell_type \
#     --batch_key donor_id
```
