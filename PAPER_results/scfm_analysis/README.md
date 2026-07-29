# TF-Scfm atlases: scFM-processed single-cell data from CELLxGENE Census

All files in this directory and the parent `/data/analysis/data_ganguly/` are derived from the **CZ CELLxGENE Census LTS 2025-11-08** with precomputed **TranscriptFormer (TF-Sapiens / TF-Exemplar)** cell embeddings (per the TranscriptFormer paper, L83).

## Available embeddings in Census 2025-11-08

| Embedding | Organism | Cells | Dim | Source |
|-----------|----------|-------|-----|--------|
| `tf-sapiens`        | homo_sapiens | 158,982,719 | 2048 | CZI (CZI-CxG-czi-17) |
| `tf-exemplar-human` | homo_sapiens | 158,982,719 | 2048 | CZI (CZI-CxG-czi-18) |
| `tf-exemplar-mouse` | mus_musculus |  43,653,561 | 2048 | CZI (CZI-CxG-czi-19) |
| `scvi`              | homo_sapiens | 158,982,719 | 50   | CZI (CZI-CxG-czi-15) |
| `scvi`              | mus_musculus |  43,653,561 | 50   | CZI (CZI-CxG-czi-16) |

## Files in `/data/analysis/data_ganguly/*_TF_scfm.h5ad`

Each h5ad contains:
- `adata.X`: raw counts (Census raw)
- `adata.obs`: cell_type, tissue, donor_id, assay, dataset_id, suspension_type
- `adata.obsm['X_tf_sapiens']`: TF-Sapiens 2048-dim (always present)
- `adata.obsm['X_tf_exemplar']`: TF-Exemplar-human 2048-dim (present for newer downloads)
- `adata.obsm['X_scvi_census']`: scVI 50-dim
- `adata.uns['scfm_source']`: "CELLxGENE Census LTS 2025-11-08"
- `adata.uns['scfm_atlas_name']`: name used in filename

| File | Size | Cells | Cell types | Donors | TF-Exemplar |
|------|------|-------|-----------|--------|-------------|
| `blood_TF_scfm.h5ad`            | 1.2 GB |  66,985 |   9 | 24 | no |
| `bone_marrow_TF_scfm.h5ad`      | 0.8 GB |  54,753 |   8 | 66 | no |
| `hpancreas_TF_scfm.h5ad`        | 0.4 GB |  20,000 |  25 | 10 | yes |
| `liver_TF_scfm.h5ad`            | 1.7 GB |  75,104 |  14 |  6 | no |
| `lung_TF_scfm.h5ad`             | 2.2 GB |  71,752 |  77 | 12 | no |
| `myeloid_TF_scfm.h5ad`          | 0.3 GB |  10,000 |  19 | 13 | yes |
| `pbmc_10k_TF_scfm.h5ad`         | 0.4 GB |  11,574 |   9 |  1 | yes |
| `covid19_TF_scfm.h5ad`          | 16.1 GB | 600,929 |  11 | 80 | yes |
| `tabula_sapiens_2_TF_scfm.h5ad` | 22.8 GB | 549,027 | 162 |  9 | yes |

## How these were obtained

All atlases were pulled with the CELLxGENE Census Python API (no model weights downloaded, no GPU needed):

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="dataset_id == '<uuid>' and is_primary_data == True",
        obs_column_names=["cell_type", "tissue", "donor_id", "assay", "dataset_id"],
        obs_embeddings=["tf-sapiens", "tf-exemplar-human", "scvi"],
    )
    adata.write_h5ad("<atlas>_TF_scfm.h5ad")
```

The scripts that produced these files are:
- `TF_scfm/download_census_tf.py` — tissue-matched analogues (blood, bone_marrow, liver, lung)
- `TF_scfm/download_tf_eval_datasets.py` — TF paper evaluation datasets (pbmc_10k, covid19, tabula_sapiens_2)
- `TF_scfm/download_pancreas_myeloid.py` — hpancreas, myeloid

## Verification

`TF_scfm/verify_all.py` loads each h5ad and confirms:
- `X_tf_sapiens` shape is (n_obs, 2048), no NaNs
- `X_scvi_census` shape is (n_obs, 50)
- `X_tf_exemplar` shape is (n_obs, 2048) where present
- Embedding values are sensible (mean ≈ 0, std > 0, no NaN)
- `uns['scfm_source']` = "CELLxGENE Census LTS 2025-11-08"

All 9 files passed.

## Datasets in the TF paper that we could NOT get from Census

The TF paper used additional evaluation datasets that are NOT in the CELLxGENE Census:

| Dataset | Source | Status |
|---------|--------|--------|
| Tabula Sapiens 2.0 (TSP17-30) | Census | done (549k cells) |
| Spermatogenesis (9 vertebrate species) | GEO + Bgee (E-MTAB-11063..73) | not in Census |
| Tabula Microcebus (mouse lemur) | Tabula Microcebus Consortium | not in Census |
| Tropical clawed frog cell atlas | GSE113074 | not in Census |
| Zebrafish cell landscape | GEO | not in Census |
| Stony coral cell atlas | GSE166901 | not in Census |
| Sea lamprey brain | E-MTAB-11087 | not in Census |
| Tahoe-100M (drug perturbations) | HuggingFace | not in Census |
| COVID-19 human lung | Census | done (600k cells) |
| Innate immune LPS (mouse, rat, rabbit, pig) | E-MTAB-6773 | not in Census |

For these, the TF model weights can be downloaded from `pip install transcriptformer` + `transcriptformer download tf-sapiens` and the model can be run locally, but each run needs a GPU.
