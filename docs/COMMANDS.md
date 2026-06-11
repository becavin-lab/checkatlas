# CheckAtlas Commands Reference

## Quick navigation

- [Single Atlas](#single-atlas) — all task/process combinations for one `.h5ad` / `.rds` file
- [Multi Atlas (batch)](#multi-atlas-batch) — run the full pipeline across a folder of atlases
- [Available Metrics](#available-metrics) — reference lists of every metric per category
- [Global Options](#global-options) — common flags shared across all commands

---

## Single Atlas

Every single-atlas command expects the atlas file to live at `<search-dir>/<atlas-name>.h5ad`
(or `.rds` / `.qs`).  The folder structure is created automatically on first run.

### 1. Single task — all metrics

Run **one** process category with every available metric from its default list.

#### Cluster metrics (all 5)

```bash
checkatlas metric_cluster <search-dir> --atlas_name <atlas-name>
```

Runs: `silhouette`, `davies_bouldin`, `calinski_harabasz`, `dbcvi`, `kolmogorov_smirnov`
across every detected embedding × cluster-label combination.

#### Annotation metrics (all 16)

```bash
checkatlas metric_annot <search-dir> --atlas_name <atlas-name>
```

Runs ref-vs-pred, embedding-based, batch/integration, and graph-connectivity metrics
across every detected reference, predicted, embedding, and batch column.

#### Dimred metrics (all 9)

```bash
checkatlas metric_dimred <search-dir> --atlas_name <atlas-name>
```

Compares every `.obsm` embedding key against `adata.X` using:
`kruskal_stress`, `spearman_rho`, `entourage`, `coknn`, `dCor`, `den_pre`, `lcmc`,
`avg_jaccard_dis`, `ged`.

#### QC plots and tables

```bash
checkatlas qc <search-dir> --atlas_name <atlas-name>
```

#### Summary tables + UMAP/t-SNE figures

```bash
checkatlas summary <search-dir> --atlas_name <atlas-name>
```

---

### 2. Single task — selective metrics

Pick specific metrics from a category so only those are computed.

```bash
# Cluster: only silhouette and dbcvi
checkatlas metric_cluster <search-dir> --atlas_name <atlas-name> \
    --metric_cluster silhouette dbcvi

# Annotation: only rand_index and lisi
checkatlas metric_annot <search-dir> --atlas_name <atlas-name> \
    --metric_annot rand_index lisi

# Dimred: only spearman_rho
checkatlas metric_dimred <search-dir> --atlas_name <atlas-name> \
    --metric_dimred spearman_rho
```

Skip a category entirely with `none`:

```bash
checkatlas metric_annot <search-dir> --atlas_name <atlas-name> \
    --metric_annot none
```

---

### 3. Multiple tasks — selective or all metrics

#### All three metric categories together (cluster + annot + dimred)

```bash
checkatlas metric <search-dir> --atlas_name <atlas-name>
```

This is the most common command for a metrics run.  All metrics from every
category are computed with their defaults.  Supply selective metric lists to
narrow each category:

```bash
checkatlas metric <search-dir> --atlas_name <atlas-name> \
    --metric_cluster silhouette davies_bouldin \
    --metric_annot rand_index normalized_mutual_info \
    --metric_dimred spearman_rho
```

#### Full analysis (summary + QC + all metrics)

```bash
checkatlas analyse <search-dir> --atlas_name <atlas-name>
```

Runs summary tables, UMAP/t-SNE figures, QC tables/plots, and all three metric
categories.  Equivalent to running `summary`, `qc`, and `metric` sequentially.

---

### 4. Preprocess only

Builds the precomputed context (kNN graphs, neighbour graphs, distance matrices)
once so subsequent metric runs can reuse it from cache.

```bash
checkatlas preprocess <search-dir> --atlas_name <atlas-name>
```

Precomputed data is saved under `<search-dir>/checkatlas_files/temp/<atlas-name>/`
and is automatically detected by `metric`, `metric_cluster`, `metric_annot`,
and `analyse` on subsequent runs.

---

## Multi Atlas (batch)

Process every atlas found in a folder via the Nextflow pipeline.  The folder is
scanned recursively for `.h5ad`, `.rds`, and `.qs` files.

### Run the full Nextflow pipeline on a folder

```bash
checkatlas run <search-dir>
```

The pipeline auto-discovers every atlas and runs QC, summary, and metrics
for each one in parallel.  All metric categories use their full default lists.

Add Nextflow execution reports:

```bash
checkatlas run <search-dir> \
    --with_report my_report.html \
    --with_dag dag.dot \
    --with_timeline timeline.html
```

### List atlases in a folder (Nextflow workflow)

Use `nfcheckatlas` to inventory atlases before running the pipeline:

```bash
# List all Scanpy .h5ad atlases
nfcheckatlas list_scanpy <search-dir>

# List all CellRanger outputs
nfcheckatlas list_cellranger <search-dir>

# List all Seurat .rds / .qs atlases
nfcheckatlas list_seurat <search-dir>
```

### Generate HTML QC report from a completed run

```bash
nfcheckatlas html_report <search-dir>
```

---

## Available Metrics

### Cluster metrics (`--metric_cluster`)

```
silhouette            davies_bouldin        calinski_harabasz
dbcvi                 kolmogorov_smirnov
```

*Default: all 5 are run when `--metric_cluster` is omitted.*
*Use `--metric_cluster none` to skip.*

### Annotation metrics (`--metric_annot`)

```
rand_index            adj_rand_index        mutual_info
normalized_mutual_info  adj_mutual_info     fowlkes_mallows
vmeasure              isolated_f1_score
average_silhouette_width  dunn_index
graph_connectivity
cell_cycle_conservation  highly_variable_genes
kbet                  lisi                  pcr
```

*Default: all 16 are run when `--metric_annot` is omitted.*
*Use `--metric_annot none` to skip.*

### Dimensionality-reduction metrics (`--metric_dimred`)

```
kruskal_stress        spearman_rho          entourage
coknn                 dCor                  den_pre
lcmc                  avg_jaccard_dis       ged
```

*Default: all 9 are run when `--metric_dimred` is omitted.*
*Use `--metric_dimred none` to skip.*

> **Note:** `trustworthiness` and `continuity` are temporarily disabled due
> to GPU memory constraints on atlases with > 100 000 cells.  They will be
> re-enabled once the chunked GPU path supports float16.

---

## Global Options

These flags apply to all `checkatlas` processes:

| Flag | Description | Default |
|------|-------------|---------|
| `-d`, `--debug` | Print debug-level log messages | off |
| `-v`, `--version` | Display the CheckAtlas version | — |
| `--n_jobs N` | Max CPU threads for parallel work (capped at 48) | 48 |

### Distinguish `checkatlas` vs `nfcheckatlas`

| Command | Purpose |
|---------|---------|
| `checkatlas` | Single-atlas CLI — all processes except `run` and the list/html_report workflows |
| `checkatlas run` | Launches the **Nextflow pipeline** on a folder of atlases |
| `nfcheckatlas` | Nextflow helper — list atlases in a folder or generate HTML reports |
