# scFM test datasets: recommended sources

This document lists the most useful scRNA-seq datasets for evaluating
the CheckAtlas scFM QC layer. Each entry is a published dataset where
a foundation model has been benchmarked and the data is publicly
available for download.

**Datasets are referenced from `docs/scfm_benchmarking_report.md` and
verified against the original papers via paperclip.**

The test suite in `tests/test_scfm_metrics.py` builds a synthetic
AnnData that mimics the structure of these real datasets (celltype
column, batch column, multiple embeddings, etc.), so the scFM
pipeline can be tested without downloading the 891 MB Kedzierska
data.

---

## 1. Pancreas (scIB benchmark) — Kedzierska 2024 / 2025

**Source:** Kedzierska et al. 2024 / 2025 — *Foundation models in
single-cell biology: evaluating zero-shot capabilities* (PMC12007350).

**Why it matters:** This is the canonical dataset used to benchmark
Geneformer, scGPT, scVI, and Harmony on batch integration. It is the
single dataset most often cited in the scFM failure literature.

**Citation:** Kedzierska KZ, Crawford L, Amini A, Lu AX. (2025)
"Zero-shot evaluation reveals limitations of single-cell foundation
models." *Genome Biology* PMC12007350.

**Data location:**
- **Processed h5ad**: `tests/data/scfm/pancreas_scib.h5ad` (local
  copy, 316 MB)
- **Original**: Kedzierska Zenodo deposit
  `10.5281/zenodo.15123462` (code + raw loom files)
- **scIB pancreas**: Luecken et al. 2020 figshare
  `10.6084/m9.figshare.12420968.v8`

**Schema:**
```
adata.shape = (16382, 19093)
adata.obs:
  celltype: 14 cell types (categorical)
  tech: 9 sequencing technologies (categorical)
  batch: 6 batch groups (categorical)
  size_factors: float (per-cell)
adata.layers:
  counts: raw UMI counts
```

**Embeddings to add for scFM evaluation:**
- `X_pca` — baseline (40-50 PCs)
- `X_scvi` — baseline (10-20 dim)
- `X_harmony` — baseline (corrected PCA, 20-30 dim)
- `X_geneformer` — scFM, 256-512 dim, cell-mode embedding
- `X_scgpt_human` — scFM, 512 dim, CLS token output
- `X_scgpt_blood`, `X_scgpt_kidney` — scFM variants

**How to obtain the scFM embeddings:** Use the Kedzierska code at
`https://github.com/microsoft/zero-shot-scfoundation` (also on
Zenodo). Both `notebooks/Geneformer_zero_shot.ipynb` and
`notebooks/scGPT_zero_shot.ipynb` produce the embeddings and add
them to `adata.obsm`.

**Eval to run:**
```bash
checkatlas scfm tests/data/scfm/ \
    --atlas_name pancreas_scib \
    --scfm_embedding X_geneformer \
    --baseline_embeddings X_pca X_scvi \
    --ref_label celltype \
    --batch_key batch \
    --n_seeds 5
```

This is the **highest-priority test case** for the scFM layer
because it covers Problems 1 (zero-shot), 2 (baselines), and 4
(batch effects) directly.

---

## 2. Tabula Sapiens — multiple papers

**Source:** Pisco et al. 2021 — *Tabula Sapiens single-cell dataset*
(figshare `10.6084/m9.figshare.14267219.v5`).

**Why it matters:** Multi-organ human cell atlas used as
benchmarking data for Geneformer, scGPT, scFoundation, UCE, and
SCimilarity. The Kedzierska 2024 paper used Tabula Sapiens to
evaluate cross-tissue generalization (Problem 6) and the
overlap with scGPT pretraining data is a known confound that the
report explicitly addresses.

**Citation:** Tabula Sapiens Consortium. (2021) "The Tabula Sapiens:
A multiple-organ, single-cell transcriptomic atlas of humans."
*bioRxiv*.

**Data location:**
- **figshare**: `10.6084/m9.figshare.14267219.v5` (h5ad files per
  organ)
- **CELLxGENE**: searchable as collection "Tabula Sapiens"
- **Kedzierska deposit**: included in the figshare bundle
  `43480497`

**Why this is good for CheckAtlas:**
- Multiple tissues (24+) — tests cross-domain ARI (Problem 6)
- Multiple donors — tests patient ASW (Problem 7)
- Author-curated cell types — tests ref-vs-pred ARI/F1 (Problem 1)
- ~500k cells total — realistic scale

**Embeddings to add:**
- `X_pca`, `X_scvi` (baselines)
- `X_geneformer`, `X_scgpt_human` (scFMs)

**Caveat:** Tabula Sapiens is in the scGPT pretraining corpus, so
zero-shot evaluation has a known overlap (Kedzierska 2024 reports
this as a confound).

---

## 3. PBMC 10k / 12k — multiple papers

**Source:** Zheng et al. 2017 (10x Genomics), widely re-distributed
via CELLxGENE and Scanpy tutorials.

**Why it matters:** The single most commonly used dataset for
scFM benchmarking. scGPT's strongest results are on PBMC 12k
(see Kedzierska 2024 Fig 1B).

**Data location:**
- **10x Genomics**: `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/`
- **CELLxGENE**: searchable as "pbmc_10k" or "pbmc_12k"
- **scanpy.datasets.pbmc3k()**: bundled with the Scanpy package

**Why this is good for CheckAtlas:**
- Multiple donors (8+ in 10k v3) — tests batch effects (Problem 4)
- 10-15 well-defined cell types — tests rare-type recall (Problem 5)
- Already preprocessed and h5ad-ready
- Small enough to download and test quickly (~50 MB)

**Embeddings to add:**
- `X_pca`, `X_scvi` (baselines)
- `X_geneformer`, `X_scgpt_human` (scFMs)

**Eval to run:**
```bash
checkatlas scfm tests/data/scfm_pbmc/ \
    --atlas_name pbmc_10k \
    --scfm_embedding X_geneformer \
    --baseline_embeddings X_pca X_scvi \
    --ref_label celltype \
    --batch_key donor \
    --n_seeds 5
```

This is the **second-highest-priority test case** for the scFM
layer because it covers Problem 1 (zero-shot F1) and Problem 5
(rare-type recall) directly.

---

## 4. Immune cell atlas (330k) — Kedzierska 2024

**Source:** Kedzierska 2024 — *Assessing the limits of zero-shot
foundation models in single-cell biology* (preprint, 2024).

**Why it matters:** Used by Kedzierska 2024 to show that scFM
embeddings fail to remove batch effects in large multi-donor
datasets. This is the empirical evidence behind Problem 4 in
our report.

**Data location:**
- **figshare**: `10.6084/m9.figshare.24747228` (bundled with the
  Kedzierska deposit)
- **CELLxGENE**: searchable by donor / tissue

**Why this is good for CheckAtlas:**
- 330k cells across 16 datasets — large enough for stability
  testing (Problem 8)
- Strong batch effect — clear kBET / iLISI signal (Problem 4)
- Multiple cell types including rare populations (Problem 5)

---

## 5. BCC and LIHC datasets — CellFM 2025

**Source:** CellFM 2025 (PMC12092794) — Zeng et al. *CellFM: a
large-scale foundation model pre-trained on transcriptomics of
100 million human cells*.

**Why it matters:** Explicitly deposited in h5ad format with
detailed cell-type annotations including exhausted and activated
CD8+ T cells. The CellFM paper compares against scGPT, Geneformer,
scFoundation, UCE.

**Data location:**
- **GEO**: `GSE123813` (BCC, basal cell carcinoma)
- **GEO**: `GSE140228` (LIHC, liver hepatocellular carcinoma)
- Both available as h5ad via GEO accession

**Why this is good for CheckAtlas:**
- h5ad format with rich cell-type annotations
- Multiple rare cell-type comparisons
- Direct head-to-head benchmark of 4 scFMs

---

## 6. PBMC, hPancreas, hPBMC — BioLLM 2025

**Source:** BioLLM 2025 (PMC12365531) — Qiu et al. *BioLLM: A
standardized framework for integrating and benchmarking
single-cell foundation models*.

**Why it matters:** The BioLLM paper uses 13 datasets including
hPancreas, hPBMC, humanDC, Zheng68K, COVID-19, Lung-Kim, etc. It
benchmarks scGPT, Geneformer, scFoundation, scBERT on cell
annotation, embedding quality, and drug response prediction.

**Data location:**
- **GitHub**: `https://github.com/BGIResearch/BioLLM` (code and
  config)
- **Zenodo**: model checkpoints also on Zenodo
- The 13 datasets are listed in the paper's `Table S2` and the
  download links are referenced in the supplementary

**Why this is good for CheckAtlas:**
- Direct head-to-head comparison of 4 scFMs
- Multiple tissue / dataset coverage
- Code is actively maintained (BGI Research)

---

## How the test fixtures in this repo map to the real datasets

| Local test fixture | Real dataset equivalent | Use case |
|---|---|---|
| `tests/test_scfm_metrics.py::small_adata` (synthetic) | PBMC 3-class synthetic | Unit test: imbalanced 3-class with batch effect |
| `tests/data/scfm/pancreas_scib.h5ad` (316 MB, real) | Kedzierska Pancreas | Integration test: real scFM-eval data |
| Synthetic, engine-built from `pancreas_scib.h5ad` subset | Kedzierska Tabula Sapiens | Future: tests cross-domain (Problem 6) |
| Synthetic, engine-built from `pancreas_scib.h5ad` subset | Kedzierska PBMC 12k | Future: tests rare-type recall (Problem 5) |

The `small_adata` fixture in `tests/test_scfm_metrics.py` is the
**synthetic** version used for unit tests. It mimics the structure
of Pancreas and PBMC: imbalanced 3-class labels, two batch
factors, two embeddings (one noisy, one clean). It runs without
downloading any external data, so the unit tests are fast and
self-contained.

The 316 MB `pancreas_scib.h5ad` is a **real** processed scFM-eval
dataset, included for integration tests and for users who want to
reproduce the Kedzierska results. It can be obtained from the
Kedzierska Zenodo deposit (see §1 above).

---

## Quick-start: getting real scFM-processed h5ad files

The minimum to test the scFM layer end-to-end:

1. **Download the Kedzierska Pancreas h5ad** (already in this
   repo at `tests/data/scfm/pancreas_scib.h5ad`):
   ```bash
   curl -L -A "Mozilla/5.0" "https://ndownloader.figshare.com/files/43480497" \
       -o kedzierska_data.zip
   unzip kedzierska_data.zip -d kedzierska_data
   cp kedzierska_data/data/datasets/pancreas_scib.h5ad \
      tests/data/scfm/
   ```

2. **Run Geneformer and scGPT zero-shot** using the Kedzierska
   notebooks (this requires the model weights and GPU):
   ```bash
   git clone https://github.com/microsoft/zero-shot-scfoundation
   cd zero-shot-scfoundation
   jupyter notebook notebooks/Geneformer_zero_shot.ipynb
   jupyter notebook notebooks/scGPT_zero_shot.ipynb
   ```
   Each notebook adds `X_geneformer` and `X_scgpt` to
   `adata.obsm` and saves the result as an h5ad.

3. **Run scVI and PCA** as baselines:
   ```python
   import scanpy as sc
   import scvi
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)
   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
   sc.pp.scale(adata, max_value=10)
   sc.tl.pca(adata, n_comps=50)
   adata.obsm['X_pca'] = adata.obsm['X_pca']
   scvi.model.SCVI.setup_anndata(adata, batch_key='batch')
   model = scvi.model.SCVI(adata)
   model.train()
   adata.obsm['X_scvi'] = model.get_latent_representation()
   adata.write_h5ad('pancreas_with_scfm.h5ad')
   ```

4. **Run CheckAtlas scFM QC**:
   ```bash
   checkatlas scfm path/to/atlas/ \
       --atlas_name pancreas_with_scfm \
       --scfm_embedding X_geneformer \
       --baseline_embeddings X_pca X_scvi \
       --ref_label celltype \
       --batch_key batch \
       --n_seeds 5
   ```

The output will land in `checkatlas_files/scfm/pancreas_with_scfm/`
with the full per-problem verdict, FMF/BF/PR scores, and a
MultiQC-friendly table.

---

## What is NOT in this repo (and why)

- **scFM model weights** (Geneformer-V2-316M, scGPT-human) — these
  are very large (multiple GB each) and require Hugging Face /
  Google Drive downloads with `git lfs`. They are not part of
  CheckAtlas. See the Kedzierska README for the exact download
  commands.

- **Tabula Sapiens full data** (~500k cells, ~5 GB) — too large to
  commit. The BioLLM GitHub repo and CELLxGENE provide direct
  download links.

- **Tabula Sapiens h5ad for the Kedzierska deposit** — the
  Kedzierska figshare `43480497` archive contains only the
  Pancreas h5ad; Tabula Sapiens is referenced as a separate
  download (see Kedzierska Table S3).

- **Other scFM outputs (e.g. scFoundation, UCE, scCello)** — the
  unified framework paper (PMC12803055) provides code to run all
  13 models, but the outputs themselves are not deposited as
  h5ad. They would need to be regenerated.
