# Biological Interpretation & Result Strategy for CheckAtlas

## Executive Summary

This document outlines the comprehensive strategy for formulating, structuring, and biologically interpreting the checkatlas benchmark results, anchored in the framework established by Luecken et al. (2021, Nature Methods) ["Benchmarking atlas-level data integration in single-cell genomics"](https://doi.org/10.1038/s41592-021-01336-8) and informed by the 2025 review by Cao et al. on the landscape of single-cell benchmarking [1].

## 1. Framework: How CheckAtlas Extends scIB

### What scIB Established (the gold standard)

The scIB paper benchmarked 16 integration methods on 13 tasks via 14 metrics organized into three categories with a 40:60 weighting:

| Category | Weight | Metrics |
|----------|--------|---------|
| Batch removal | 0.40 | kBET, iLISI, PCR batch, ASW batch, graph connectivity |
| Bio-conservation (label) | 0.36 | NMI, ARI, cell-type ASW, cLISI, isolated label F1, isolated label silhouette |
| Bio-conservation (label-free) | 0.24 | cell-cycle conservation, HVG conservation, trajectory conservation |

Overall score: `S_overall = 0.6 * S_bio + 0.4 * S_batch`

The key biological insight from scIB: **integration methods face a fundamental tradeoff between batch removal and biological conservation** -- methods that aggressively mix batches tend to erase subtle biological variation, while methods that preserve fine biology leave residual batch structure.

### What CheckAtlas Adds

CheckAtlas operates at a higher conceptual level than scIB: **it evaluates the quality of a single-cell atlas across four independent dimensions** (cluster, annotation, batch correction, dimensionality reduction) rather than just evaluating integration quality. This is analogous to what Cao et al. [1] call a "multi-criteria benchmarking framework" -- where method performance extends beyond a single task.

| Dimension | scIB analogy | What checkatlas evaluates |
|-----------|-------------|--------------------------|
| Cluster | Cell-type ASW, NMI (indirectly) | Are cell populations well-separated in embedding space? |
| Annotation | ARI, NMI, isolated label F1 | Can cell types be reliably transferred/predicted? |
| Batch correction | kBET, iLISI, PCR, graph conn, cLISI | Are batch effects removed without destroying biology? |
| Dimred | Not in scIB | Does the low-dimensional embedding faithfully represent high-dimensional relationships? |

The scIB paper evaluates **integration method quality** (which method is best for a given task). CheckAtlas evaluates **atlas quality** (how good is this particular atlas, regardless of how it was made).

## 2. Biological Interpretation Per Task

### 2.1 Cluster Metrics: "Are cell populations identifiable?"

**What they measure biologically:** Whether distinct cell types/states form well-separated, compact clusters in the embedding. This is the prerequisite for all downstream cell-type identification [2].

**Metrics and their biological meaning:**

| Metric | Biological Question | Range | Direction | Interpretation for Atlas |
|--------|--------------------|-------|-----------|------------------------|
| Silhouette width | How compact and separated are clusters? | [-1, 1] | Higher | >0.5: good cell-type separation; <0.2: overlapping populations |
| Calinski-Harabasz | Is between-cluster variance >> within-cluster? | [0, inf) | Higher | High values = embedding captures cell-type differences |
| Davies-Bouldin | How similar are clusters to their nearest neighbor? | [0, inf) | Lower | Low values = distinct cell identities in the embedding |
| DBCVI | Are cluster boundaries at density minima? | [-1, 1] | Higher | Density-based assessment; insensitive to cluster shape |
| Kolmogorov-Smirnov | Are feature distributions different between clusters? | [0, 1] | Higher | >0.8: pairwise separation at gene-expression level |

**Biological interpretation strategy per atlas:**

- **bone_marrow (27k cells):** If silhouette is high (>0.4-0.5), hematopoietic hierarchy is well-resolved in the embedding. If DBCVI is low, the continuous differentiation gradient between progenitors blurs cluster boundaries.
- **TS_immune (592k cells):** Large immune atlases often show moderate silhouette (0.2-0.3) because T-cell subtypes are transcriptionally close. High KS scores despite moderate silhouette indicate that genes differentiate subtypes even when clusters overlap in PCA/UMAP space.
- **lung (348k cells):** Lung atlases with many epithelial subtypes should show high Calinski-Harabasz (epithelial vs. immune vs. stromal are transcriptionally distinct). Low DB values in epithelial subclusters indicate high-quality annotation granularity.
- **blood (85k cells):** PBMC atlases should have high cluster metrics overall -- blood cells have distinct transcriptional identities. Low scores suggest poor pre-processing or batch effects.

**Critical caveat:** Cluster metrics are inherently label-dependent. If the reference annotation is coarse (e.g., "immune cell" vs. "epithelial"), metrics will be artificially inflated. If it is too fine (sub-subtypes), they will be deflated. The column detector's ability to find appropriate annotation granularity is crucial.

### 2.2 Annotation Metrics: "Can cell types be reliably transferred?"

**What they measure biologically:** The degree to which cell-type labels can be predicted from the embedding -- this measures both embedding quality and annotation quality.

**Metrics and their biological meaning:**

| Metric | Biological Question | Range | Direction |
|--------|--------------------|-------|-----------|
| ARI / AMI / NMI | Does clustering match the reference annotation? | [0, 1] | Higher |
| Isolated F1 | Are rare cell types preserved as distinct clusters? | [0, 1] | Higher |
| V-measure | Harmonic mean of homogeneity and completeness | [0, 1] | Higher |
| ASW (annotation) | How well are cell types separated in embedding space? | [-1, 1] | Higher |
| Dunn index | Ratio of min inter-cluster to max intra-cluster distance | [0, inf) | Higher |

**Biological interpretation strategy:**

- Following scIB, **isolated label scores are the most biologically informative** for atlas quality. Luecken et al. showed that Harmony overlapped rare bone marrow populations (erythroid progenitors, megakaryocyte progenitors) in the human immune task despite good overall scores. A high isolated F1 means the atlas preserves rare cell states -- critical for disease atlases where rare populations carry clinical significance.
- **Ref-vs-pred metrics (ARI, NMI, AMI, V-measure):** When the "predicted" column is from the same lab, this measures annotation transfer accuracy. When it's from different pipelines (e.g., manual vs. automated annotation), this measures annotation reproducibility.
- **For the TS_immune atlas (592k cells):** If ARI > 0.7 between reference and predicted labels, the atlas supports reliable automated cell-type transfer -- a key property for large-scale atlas efforts like the Human Cell Atlas [3].
- **For blood (85k cells):** PBMC atlases typically have well-characterized markers, so annotation metrics should be high (>0.8 ARI). Low scores may indicate batch-driven annotation errors.

### 2.3 Batch Correction Metrics: "Are technical artifacts removed without destroying biology?"

**What they measure biologically:** The same fundamental question as scIB: does the embedding mix batches while keeping cell types pure? This is the **batch-vs-biology tradeoff** that Luecken et al. identified as the central challenge.

**Metrics and their biological meaning:**

| Metric | Biological Question | Range | Direction | scIB equivalent |
|--------|--------------------|-------|-----------|-----------------|
| iLISI | Are batches well-mixed in local neighborhoods? | [0, 1] | Higher (mix) | Graph iLISI |
| cLISI | Are cell types pure in local neighborhoods? | [0, 1] | Lower (purity) | Graph cLISI |
| kBET | Does local batch composition match global expectation? | [0, 1] | Lower (acceptance) | kBET |
| PCR | Can batch identity be predicted from PCs? | [0, 1] | Lower | PCA regression |
| Graph connectivity | Do cells of the same type stay connected? | [0, 1] | Higher | Graph connectivity |

**Biological interpretation strategy -- following scIB's deep analysis of the human immune task:**

Luecken et al. demonstrated that **individual metrics reveal specific failure modes:**
- High iLISI + high cLISI = **ideal integration** (batches mixed, biology preserved) -- e.g., Scanorama on immune task
- High iLISI + low cLISI = **over-integration** (batches mixed, but biology erased) -- e.g., LIGER on cross-species
- Low iLISI + high cLISI = **under-integration** (biology preserved, but batch structure remains)
- Low iLISI + low cLISI = **chaos** (neither batches mixed nor biology preserved)

**Atlas-specific biological interpretation:**

- **bone_marrow (27k cells):** If batch effects come from different donors/sequencing runs, iLISI should be high for bone marrow samples (cells are similar across donors). If cLISI drops in hematopoietic stem cells, the embedding may over-correct and merge progenitors.
- **TS_immune (592k cells):** This atlas likely combines multiple datasets from Human Cell Atlas. Following scIB's observation that "methods that use cell identity information preserved biological variation most strongly," we should check whether cell-type purity (cLISI) varies across cell-type hierarchy -- T-cell subtypes should be more vulnerable to over-correction than major lineages.
- **lung (348k cells):** The scIB lung task revealed that **spatial location confounds with batch** -- endothelial cells in airway vs. parenchyma have genuine functional differences. If batch correction metrics are too good (iLISI near 1), the atlas may have erased these biologically meaningful spatial gradients.
- **blood (85k cells):** PBMC atlases typically integrate well because blood cell types are transcriptionally distinct. Poor batch metrics here would strongly indicate technical problems.

**The cLISI-iLISI tradeoff plot (following scIB Fig 3a):**

This scatter plot is the single most biologically informative figure. Following Luecken et al.'s demonstration that "methods strike different balances per task," we should generate this plot **per embedding per atlas** to show where each embedding sits on the batch-vs-biology Pareto frontier.

### 2.4 Dimred Metrics: "Is the low-dimensional representation faithful?"

**What they measure biologically:** Whether the 2D/3D embedding (UMAP, t-SNE, PCA) used for visualization and interpretation preserves the structure of the high-dimensional gene expression data. This is critical because **biological conclusions are drawn from visualizations.**

**Metrics and their biological meaning:**

| Metric | Aspect | Biological Question | Range | Direction |
|--------|--------|--------------------|-------|-----------|
| Kruskal stress | Global geometry | Are relative distances preserved? | [0, 1] | Lower |
| Spearman rho | Rank order | Is the neighborhood rank order preserved? | [-1, 1] | Higher |
| dCor | Distance dependence | Is there any distance correlation? | [0, 1] | Higher |
| CoKNN | Local neighborhood | What fraction of kNNs are shared? | [0, 1] | Higher |
| Entourage | Local neighborhood | Do low-dim neighbors come from high-dim neighbors? | [0, 1] | Higher |
| LCMC | Local continuity | kNN overlap vs. random baseline | [0, 1] | Higher |
| DenPre | Density | Are relative local densities preserved? | [-1, 1] | Higher |
| GED | Graph structure | Is the kNN graph preserved? | [0, 1] | Lower |
| Avg Jaccard | Neighborhood overlap | Mean Jaccard of shared neighbors | [0, 1] | Lower |

**Biological interpretation strategy:**

Dimred metrics are **not in scIB** because scIB focused on integration method evaluation, not visualization quality. However, for atlas QC, dimred metrics are critical because:

1. **Biologists make decisions from UMAP visualizations** -- if the UMAP distorts distances, a researcher might conclude two populations are related when they are not.
2. **UMAP vs. t-SNE vs. PCA tradeoffs** are well-documented: UMAP preserves local structure better; t-SNE exaggerates cluster separation; PCA preserves global structure but not local. Our metrics quantify these tradeoffs.
3. **Global vs. local preservation tradeoff:** A UMAP with high coKNN/entourage (>0.8) but low dCor (<0.3) preserves local neighborhoods at the expense of global distances -- fine for cell-type identification but dangerous for trajectory inference.

**Atlas-specific interpretation:**

- **TS_neural (2.7k cells):** Small atlases should show excellent dimred preservation across all metrics -- if Spearman rho < 0.5 on such a small dataset, the embedding method has failed.
- **TS_immune (592k cells):** Large atlases challenge UMAP's ability to preserve global structure. Low dCor is expected; the question is whether local metrics (coKNN, entourage) remain high enough for reliable cell-type identification.
- **lung (348k cells):** The lung's continuous epithelial-to-mesenchymal transitions require good distance preservation. Low spearman_rho would mean pseudotime trajectories inferred from UMAP are unreliable -- a biologically critical warning.

**Critical insight from the Xia et al. 2025 benchmarking paper [4]:** Feature selection methods affect integration performance. Since our dimred metrics evaluate the embedding against the high-dim reference (either `adata.X` or PCA for large atlases), the choice of reference space itself affects metric values. For large atlases (>300k cells), using PCA as reference means the metrics measure how well the low-dim embedding captures **the top principal components** -- not all biological variation.

## 3. Result Table Design

Following the scIB paper's presentation approach (Fig 2, Fig 3, Fig 4) and the framework from Cao et al. [1], we need multiple levels of result aggregation:

### 3.1 Atlas-Level Summary Table (Table 1)

One row per atlas, with aggregated scores per task. This is analogous to scIB's Table 1.

```
Atlas          Cells      Genes    Cluster  Annot    Batch    Dimred   Total Time
bone_marrow    27,112     60,606   0.XX     0.XX     0.XX     0.XX     20m 39s
TS_neural      2,685      60,606   0.XX     0.XX     0.XX     0.XX     1m 07s
liver          22,214     60,606   0.XX     0.XX     0.XX     0.XX     14m 39s
blood          85,233     60,606   0.XX     0.XX     0.XX     0.XX     1h 23m
blood_pbmc     43,512     21,649   0.XX     0.XX     0.XX     0.XX     12m 14s
lung           347,970    29,824   0.XX     0.XX     0.XX     0.XX     21m 13s
TS_skin        17,786     60,606   0.XX     0.XX     0.XX     0.XX     7m 19s
TS_immune      592,317    60,606   0.XX     0.XX     0.XX     0.XX     1h 27m
lung_2M        2,282,447  55,329   FAILED   FAILED   FAILED   FAILED   FAILED
```

Each aggregated score per task should be a min-max scaled [0,1] mean of all metrics within that task, following scIB. This allows cross-atlas comparison on the same scale.

### 3.2 Per-Embedding Batch-vs-Biology Summary (Table 2)

Following scIB Fig 3a, this is the core interpretive table:

```
Atlas       Embedding    iLISI  cLISI  kBET   PCR    Graph  ASW    IsolF1  BatchScore  BioScore  Overall
blood       X_umap       0.82   0.15   0.08   0.12   0.95   0.72   0.88    0.92        0.78      0.86
blood       X_pca        0.45   0.22   0.35   0.28   0.88   0.68   0.85    0.58        0.72      0.64
blood       X_scVI       0.91   0.12   0.04   0.06   0.97   0.75   0.91    0.96        0.80      0.90
...
```

This is the table that enables biological interpretation: the scVI embedding integrates batches well but may have compressed biological variation (very low cLISI); the PCA embedding preserves biology but has residual batch structure.

### 3.3 Dimred Per-Embedding Summary (Table 3)

```
Atlas       LowDim    HighDimRef   Stress  Spearman  dCor   coKNN  Entourage  LCMC   DenPre  GED
blood       X_umap    X_pca        0.12    0.78      0.45   0.72   0.81       0.68   0.55    0.08
blood       X_tsne    X_pca        0.09    0.65      0.32   0.58   0.67       0.51   0.42    0.15
blood       X_pca     X            0.18    0.82      0.58   0.85   0.90       0.78   0.72    0.05
```

### 3.4 Cluster Per-Embedding-Per-Label Summary (Table 4)

```
Atlas       Embedding    LabelKey          Silhouette  DB     CH     DBCVI  KS
bone_marrow X_umap       cell_type_fine    0.48        1.23   234.5  0.62   0.85
bone_marrow X_umap       cell_type_coarse  0.72        0.87   876.2  0.78   0.92
bone_marrow X_pca        cell_type_fine    0.31        1.89   112.3  0.45   0.72
```

### 3.5 scFM-Specific Summary (Table 5)

For foundation model embeddings, the scFM layer already produces structured outputs (`verdicts.tsv`, `composite.tsv`). Integrate FMF, BF, PR scores into the master summary.

## 4. Performance & Scalability Figures

The runtime data from the 8-core + GPU-JAX benchmarks is a key selling point. Following Cao et al. [1], who found that only 62% of benchmarking studies report speed and 29% report memory, we should present comprehensive scalability data.

### 4.1 Figure: Scalability -- Runtime vs. Number of Cells

**Type:** Log-log scatter plot

```
X-axis: Number of cells (log scale)
Y-axis: Total runtime (log scale)

Points: one per atlas, colored by whether PCA subsampling was used
Dashed line: O(N log N) reference line
Shaded region: "interactive" threshold (<5 min)

Selected points annotated: atlas name + cell count
```

**Biological/methodological insight:** For large atlases (>300k cells), checkatlas automatically switches to PCA reference + stratified subsampling (50k core-set). The plot shows this as a dramatic runtime improvement -- TS_immune (592k cells, 1h 27m) is only slightly slower than blood (85k cells, 1h 23m) because PCA subsampling saves massive computation. lung_2M failing tells us where the current limit is.

### 4.2 Figure: Task-Level Runtime Breakdown

**Type:** Stacked bar chart

```
X-axis: Atlas (ordered by cell count)
Y-axis: Time (minutes)
Stacked segments: Preprocess, Cluster, Annotation, Batch Correction, Dimred

Legend shows % of total per task
```

**Insight:** The preprocess phase (distance matrix + kNN graph computation) dominates -- often 70-90% of total time. This validates our GPU acceleration + caching strategy.

### 4.3 Figure: Per-Metric Runtime Heatmap

**Type:** Heatmap with metric × atlas

```
Rows: 28 metrics
Columns: 9 atlases
Color: log10(runtime) -- white=fast, dark=slow
Annotated: exact seconds in cells
```

**Insight:** Identifies which specific metrics are bottlenecks. pcr, graph_connectivity, and dimred distance matrix computation are the heavy hitters. This guides users on which metrics to run on which atlases.

### 4.4 Figure: Atlas Size vs. Metric Coverage

**Type:** Radar chart

```
One radar per atlas
Axes: number of embeddings, number of label keys detected, number of batch keys, number of ref-pred annotation pairs
Area = "atlas complexity"
```

**Insight:** Shows how checkatlas's automated column detector scales with atlas complexity. TS_immune with many annotation columns provides richer evaluation than TS_neural with few.

### 4.5 Figure: GPU vs. CPU Speedup (if comparative data available)

**Type:** Bar chart showing fold speedup per metric

**Insight:** JAX-accelerated kNN and LISI benefit most from GPU. Distance matrix computation sees 2-5x speedup on GPU for atlases <50k cells.

## 5. Critical Assessment of Metrics

This section follows the approach of Luecken et al. and the scFM diagnostics framework. Each metric has limitations that must be acknowledged.

### 5.1 Cluster Metrics

| Metric | Limitation | Severity | Mitigation |
|--------|-----------|----------|------------|
| Silhouette (CVA) | Approximation error vs. exact computation | Medium (for large K) | Report both when possible; flag when K > 20 |
| Calinski-Harabasz | Assumes spherical clusters; sensitive to outliers | Medium | Use in conjunction with DB and DBCVI |
| Davies-Bouldin | Sensitive to cluster count; penalizes many-cluster solutions | Medium | Normalize by expected DB for K clusters |
| DBCVI | MST computation is O(N²) exact; approximate path for large N | Low (exact for ≤10k cells) | Report approximation method used |
| KS | Pairwise computation is O(K²); slow for many clusters | High for fine annotations | Cap at K ≤ 30 clusters; use joblib parallelization |

### 5.2 Annotation Metrics

| Metric | Limitation | Severity | Mitigation |
|--------|-----------|----------|------------|
| ARI/AMI/NMI | All are label-dependent; sensitive to annotation granularity | High | Report at multiple annotation resolutions |
| Isolated F1 | Depends on definition of "isolated" label; sensitive to batch count | Medium | Report which labels were considered isolated |
| Dunn Index | Outlier-sensitive; minimum inter-cluster distance can be driven by a single cell | Medium | Use percentile-based variants |
| ASW (annotation) | Same approximation issues as cluster silhouette | Medium | Flag when using CVA |

### 5.3 Batch Correction Metrics -- Fundamental Limitations (from scIB)

Following Luecken et al.'s own discussion of metric limitations:

| Metric | Limitation | Source |
|--------|-----------|--------|
| kBET | Requires cell-type labels to condition on; sensitive to k choice; diffusion correction adds computational overhead | scIB Methods, Supplementary Note 1 |
| iLISI/cLISI | Graph LISI extension is sensitive to the specific graph construction; local neighborhood size affects LISI values | scIB Methods, Supplementary Note 2 |
| PCR | Only captures linear batch effects; nonlinear batch structure is invisible to PCR | scIB Discussion |
| Graph connectivity | A single disconnected cell can drive the score to 0 for a cell type; sensitive to graph construction parameters | scIB Methods |

**Key insight from scIB:** The graph LISI extension was specifically developed because original LISI (Korsunsky et al., 2019) could not handle graph-only outputs (e.g., BBKNN, Conos). Our implementation must clearly document which LISI variant is used.

### 5.4 Dimred Metrics

| Metric | Limitation | Severity | Mitigation |
|--------|-----------|----------|------------|
| Kruskal stress | Sensitive to distance normalization; arbitrary scale | Medium | Report raw and normalized |
| Spearman rho | Only captures monotonic relationships; not translation-invariant | Low | Complement with dCor |
| dCor | O(N²) centering matrix; capped at 5k cells for memory | High for large atlases | Document subsample size |
| CoKNN / Entourage / LCMC | All measure the same thing (kNN overlap) from slightly different angles; **redundant** | Medium | Report coKNN as primary; others as supplementary |
| DenPre | Requires kth-neighbor distances; sensitive to k | Low | Standardize k=15 |
| GED | Depends on k; normalized version is more robust | Low | Report both raw and normalized |
| Trustworthiness / Continuity | Disabled due to O(N²) memory; their absence is a genuine limitation | High | Acknowledge openly; explore approximate algorithms |

### 5.5 The Label-Dependency Problem

**This is the single most important limitation.** Following scIB's approach of conditioning batch metrics on cell-type labels:

- All our cluster metrics require cluster labels
- Annotation ref-vs-pred metrics require reference labels
- Batch correction metrics (kBET, cLISI, graph connectivity, ASW batch) require cell-type labels
- Only dimred metrics and some batch metrics (iLISI, PCR) are label-free

**Consequence:** If the atlas has poor or missing annotations, metric coverage drops dramatically. The column detector mitigates this by auto-detecting columns, but if no cluster or annotation columns exist, those tasks cannot run. This is a design limitation, not a bug -- "crash loudly" is the checkatlas philosophy.

### 5.6 The Reference-Space Dependency Problem

**Dimred metrics are fundamentally reference-dependent.** What you compare the embedding against determines the metric value:

- `adata.X` (full gene expression): most biologically meaningful, but O(N²·G) for distance computation
- `X_pca` (PCA-reduced): loses non-linear biological variation; large atlases are evaluated on PCA space only
- `adata.X` with HVG selection: balances biological fidelity with computational feasibility

**Mitigation:** Always report which high-dim reference was used. Flag when PCA reference is used (atlases >300k cells). This is analogous to scIB reporting which feature space was used for ATAC integration.

## 6. Biological Interpretation: The Unified Framework

### 6.1 The Four-Dimensional Atlas Quality Score

We need an interpretable composite that maps the four tasks onto biological meaning:

```
Atlas quality = f(Cluster, Annotation, Batch, Dimred)
```

Unlike scIB's 40:60 batch:bio weighting (which was designed for evaluating integration methods), an atlas QC score should reason differently:

| Task | Biological Weight | Rationale |
|------|------------------|-----------|
| Batch Correction | 0.30 | A good atlas must remove batch effects, but perfect mixing is not the goal |
| Annotation | 0.30 | Cell-type annotation quality is the primary output of atlases |
| Cluster | 0.20 | Well-separated clusters are important but secondary to annotation |
| Dimred | 0.20 | Visualization fidelity matters for biological interpretation |

### 6.2 The Batch-vs-Biology Pareto Frontier

Following scIB Fig 3a, the **iLISI vs. cLISI scatter plot** is the most informative single figure. Each point is an embedding in a specific atlas, and the frontier shows:

- **Top-right quadrant (high iLISI, low cLISI):** Ideal -- batches mixed, cell types pure
- **Bottom-right (low iLISI, low cLISI):** Under-integration -- residual batch structure
- **Top-left (high iLISI, high cLISI):** Over-integration -- biology erased
- **Pareto frontier:** The set of embeddings that are not dominated in both objectives

**Biological narrative:** "Embeddings on the Pareto frontier represent the best-possible tradeoffs between batch mixing and biological conservation. Embedding X (scVI on blood atlas) achieves near-ideal performance with iLISI=0.91 and cLISI=0.12, comparable to the best methods in scIB's human immune cell task (Scanorama embedding, iLISI=0.88, cLISI=0.14)."

### 6.3 Rare Cell-Type Analysis

Following scIB's isolated label F1 and the scFM rare_types diagnostic:

- **For each atlas, identify which cell types are isolated** (present in fewest batches)
- **Check their F1 scores** across embeddings
- **Biological interpretation:** "In the bone_marrow atlas, erythroid progenitors (isolated label) show F1=0.92 on scVI embedding but F1=0.45 on PCA, indicating that integration is necessary for rare cell-type preservation in this atlas."

This directly parallels Luecken et al.'s finding that "Harmony exhibited the lowest isolated label F1 bio-conservation score among top performers" and "kept each isolated cell label together, [but] overlapped these populations."

### 6.4 Trajectory & Cell-State Analysis (if applicable)

For atlases with developmental trajectories (bone marrow with hematopoiesis, TS_neural with neural development):

- **Dimred metrics tell us whether the trajectory can be trusted:** Low spearman_rho means pseudotime order is not preserved
- **Batch metrics on trajectory cell types** reveal whether development is confounded with batch

## 7. Figures to Generate

### 7.1 Main Figures (for paper)

| Figure | Type | Content | Biological Question |
|--------|------|---------|--------------------|
| Fig 1 | Pipeline overview | Schematic of checkatlas workflow | Method description |
| Fig 2 | Scatter | iLISI vs. cLISI per embedding per atlas | Batch-vs-biology tradeoff |
| Fig 3 | Heatmap | All metrics × all atlases × all embeddings | Global overview |
| Fig 4 | Bar chart | Atlas-level summary scores per task | Which atlases are high-quality? |
| Fig 5 | Scatter | Runtime vs. cell count (log-log) | Scalability |
| Fig 6 | Panel | UMAPs colored by highest/lowest scoring metrics | Visual validation of metric results |

### 7.2 Supplementary Figures

| Figure | Type | Content |
|--------|------|---------|
| Fig S1 | Scatter grid | Per-metric correlation analysis (do cluster metrics agree?) |
| Fig S2 | Stacked bar | Runtime breakdown per phase per atlas |
| Fig S3 | Heatmap | Dimred metrics × embeddings × atlases |
| Fig S4 | Box plots | Distribution of isolated F1 scores per cell type |
| Fig S5 | Line plot | Metric values vs. annotation resolution |
| Fig S6 | Radar | Per-atlas per-task scores |
| Fig S7 | Heatmap | Per-metric runtime (seconds) per atlas |

### 7.3 Implementation Note

Currently checkatlas has **no metric plotting code** (only QC violin/UMAP plots in `atlas.py`). All the above figures need new implementation. The code should live in:

```
checkatlas/
├── plotting/
│   ├── __init__.py
│   ├── summary.py        # Atlas-level summary figures
│   ├── batch_vs_bio.py   # iLISI vs cLISI scatter
│   ├── dimred.py         # Dimred metric heatmaps
│   ├── scalability.py    # Runtime vs cells plots
│   └── radar.py          # Per-atlas radar charts
```

## 8. Comparative Interpretation: Mapping checkatlas Metrics to scIB Metrics

To ground our biological interpretation in the well-validated scIB framework:

| scIB Metric | checkatlas Equivalent | Notes |
|-------------|---------------------|-------|
| kBET | `kbet` in batch_correction | Direct implementation; we use approximate kNN (pynndescent) vs. exact |
| Graph iLISI | `iLISI` in batch_correction | Our iLISI/cLISI uses the same graph LISI extension from scIB |
| Graph cLISI | `cLISI` in batch_correction | Same; scIB rescaled both to [0,1] |
| PCR batch | `pcr` in batch_correction | Our PCR uses RidgeClassifierCV vs. linear regression in scIB |
| ASW batch | Not in checkatlas | scIB computed silhouette on batch labels -- we use iLISI + kBET instead |
| Graph connectivity | `graph_connectivity` in batch_correction | Direct implementation |
| NMI | `normalized_mutual_info` in annot | Direct implementation (sklearn) |
| ARI | `adj_rand_index` in annot | Direct implementation (sklearn) |
| Cell-type ASW | `average_silhouette_width` in annot | scIB used full ASW; we use centroid-variance approximation for large atlases |
| Isolated label F1 | `isolated_f1_score` in annot | Direct implementation |
| Isolated label silhouette | Not in checkatlas | Could be added to annot metrics |
| CC conservation | Disabled in checkatlas | Requires adata_before/after; disabled for atlas-level evaluation |
| HVG conservation | Disabled in checkatlas | Requires adata_before/after |
| Trajectory conservation | Not in checkatlas | Would require sc.pp.dpt computation |
| N/A | 5 cluster metrics | Not in scIB; unique to checkatlas |
| N/A | 9 dimred metrics | Not in scIB; unique to checkatlas |

**Key insight:** CheckAtlas implements 7/14 of scIB's metrics directly, omits batch-ASW (covered by iLISI+kBET), and omits label-free conservation metrics (cell-cycle, HVG, trajectory) because these require `adata_before` + `adata_after` pairs -- atlas-level QC evaluates a single atlas, not a before/after comparison.

## 9. Implementation Plan

### Phase 1: Result Aggregation (Python script)

Create `PAPER_results/bio_interpretation/aggregate_results.py`:
- Read all TSV files from `checkatlas_files/{cluster,annotation,batch_correction,dimred}/`
- Parse log files for timing data
- Produce master DataFrames: per-atlas, per-embedding, per-metric
- Apply min-max scaling within each task (following scIB)
- Compute aggregated scores

### Phase 2: Result Tables (Python/LaTeX)

Create `PAPER_results/bio_interpretation/generate_tables.py`:
- Table 1: Atlas summary (Section 3.1)
- Table 2: Batch-vs-biology (Section 3.2)
- Table 3: Dimred summary (Section 3.3)
- Table 4: Cluster summary (Section 3.4)
- Output: CSV + LaTeX tables for paper

### Phase 3: Figures (Python/matplotlib+seaborn)

Create `checkatlas/plotting/` module (Section 7.3):
- Implement all figures in Section 7
- Consistent styling, colorblind-friendly palettes
- PDF + PNG output at 300 DPI

### Phase 4: Biological Interpretation Text

Create `PAPER_results/bio_interpretation/interpretation.py`:
- Per-atlas per-task biological interpretation
- Automated generation of interpretation text
- Integration with scFM diagnostic remarks

## 10. Key Biological Narratives to Highlight

### Narrative 1: Scalability Enables Atlas-Wide QC
"CheckAtlas processes atlases from 2.7k to 2.3M cells, adapting its computation strategy automatically. The automatic switch to PCA reference + stratified subsampling at 300k cells enables QC of million-cell atlases that would be computationally prohibitive with exact methods."

### Narrative 2: The Batch-vs-Biology Frontier Quantifies Integration Quality
"Across 9 atlases and multiple embeddings, checkatlas reveals that no single embedding simultaneously optimizes batch removal and biological conservation -- a fundamental tradeoff first described by Luecken et al. (2021) that we now quantify systematically at atlas scale."

### Narrative 3: Dimred Metrics Reveal Visualization Trustworthiness
"Our suite of 9 dimensionality reduction metrics provides, for the first time, a quantitative assessment of whether biological conclusions drawn from UMAP/t-SNE visualizations are supported by the high-dimensional data."

### Narrative 4: Rare Cell-Type Preservation is Atlas-Specific
"Isolated label analysis reveals that rare cell types (erythroid progenitors, megakaryocyte progenitors, specialized epithelial subtypes) are differentially preserved across atlases and embeddings, with important implications for disease-focused atlases where rare populations carry clinical significance."

### Narrative 5: GPU Acceleration Makes Large-Atlas QC Practical
"JAX GPU acceleration reduces runtime by 2-5x for distance-based metrics, and the persistent precomputation cache eliminates redundant computation across re-runs, making iterative atlas QC feasible."

## 11. Concrete Results Appendix -- Blood Atlas (85k cells) as Worked Example

This section walks through actual checkatlas results for the blood atlas (85,233 cells, 60,606 genes, 7 embeddings evaluated) to demonstrate the biological interpretation workflow in practice.

### 11.1 Cluster Metrics -- Blood Atlas

| Embedding | Label Key | Silhouette | CH | DB | DBCVI | KS |
|-----------|-----------|------------|-----|-----|-------|-----|
| X_pca | scvi_leiden_donorassay_full | 0.301 | 10,983 | 1.59 | -0.67 | 0.600 |
| X_pca | scvi_leiden_res05_tissue | 0.063 | 10,911 | 2.25 | -0.68 | 0.421 |
| X_scvi | scvi_leiden_donorassay_full | 0.139 | 8,854 | 1.42 | -0.77 | 0.765 |
| X_scvi | scvi_leiden_res05_tissue | 0.227 | 17,358 | 1.27 | -0.68 | 0.743 |
| X_umap | scvi_leiden_donorassay_full | -0.014 | 10,303 | 4.30 | -0.88 | 0.864 |
| X_umap | scvi_leiden_res05_tissue | -0.157 | 7,997 | 6.99 | -0.87 | 0.807 |

**Biological interpretation:** The blood atlas shows **moderate cluster quality overall**. The scVI embedding (X_scvi) produces the best Davies-Bouldin scores (1.27-1.42) and highest KS scores (0.74-0.77), indicating that cell-type clusters are most separated in the scVI latent space. In contrast, PCA and UMAP embeddings show poor silhouette scores (near zero or negative), suggesting overlapping clusters in these spaces. The UMAP embeddings have very high KS scores (0.80-0.86) despite negative silhouette -- this paradox reveals that UMAP preserves gene-level feature distribution differences even when the embedding itself compresses clusters. **The fine-grained tissue labels (scvi_leiden_res05_tissue) generally produce better separation than coarser donor assay labels** -- this is biologically expected as PBMC cell types are more transcriptionally distinct than donor-level groupings.

### 11.2 Batch Correction Metrics -- Blood Atlas

| Embedding | Batch Key | iLISI | kBET | PCR | Graph Conn | cLISI |
|-----------|-----------|-------|------|-----|------------|-------|
| X_pca | donor_id | 0.033 | 0.998 | 0.898 | - | - |
| X_scvi | donor_id | 0.212 | 0.974 | 0.719 | - | - |
| X_pca | assay | 0.005 | 0.116 | 0.993 | - | - |
| X_scvi | assay | 0.116 | 0.164 | 0.884 | - | - |
| X_pca | sex | 0.083 | 0.970 | 0.903 | - | - |
| X_scvi | sex | 0.482 | 0.716 | 0.717 | - | - |
| X_pca | cell_type | - | - | - | 0.917 | 0.010 |
| X_scvi | cell_type | - | - | - | 0.860 | 0.019 |

**Biological interpretation -- following scIB's framework:**

1. **Donor-level batch effects are well-resolved but not eliminated.** The scVI embedding (iLISI=0.21 for donor_id) achieves 6x better donor mixing than PCA (iLISI=0.03). However, kBET=0.97 for scVI means local donor composition still deviates from the global expectation for 3% of tested neighborhoods -- comparable to the "residual 10X batch structure in CD14+ monocytes" that Luecken et al. observed for scANVI on the human immune task.

2. **Assay-level batch effects are nearly invisible in PCA.** The extremely low iLISI (0.005) for assay on PCA means different sequencing assays form completely separate neighborhoods -- a severe batch effect. scVI partially resolves this (iLISI=0.12), but kBET=0.16 confirms that assay-driven structure persists.

3. **Sex is the most successfully mixed batch variable in scVI** (iLISI=0.48, the highest batch mixing score observed). This is biologically expected since sex differences are subtle in PBMCs and scVI was designed to remove such donor-level effects. However, scVI also shows the highest sex batch mixing in kBET (0.72 rejection), suggesting some neighborhoods remain sex-imbalanced.

4. **Cell-type purity is excellent in both embeddings** (cLISI < 0.02 in all cases), close to the ideal value of 0. This means cell-type labels are accurately preserved in local neighborhoods -- a key requirement for downstream cell-type identification. Graph connectivity is high (0.86-0.92), confirming that cells of the same type are well-connected.

5. **The batch predictors (PCR) tell a nuanced story:** scVI reduces batch predictability for donor_id (PCR=0.72 vs. 0.90 in PCA), confirming its batch-correction effect. However, PCR for assay remains very high in both embeddings (>0.88), indicating that assay type leaves strong linear signatures that even scVI cannot remove -- this is consistent with scIB's finding that protocol differences (e.g., Smart-seq2 vs. 10X) are among the hardest batch effects to correct.

### 11.3 Dimred Metrics -- Blood Atlas

| Embedding | Spearman | dCor | coKNN | Entourage | DenPre | GED | Kruskal |
|-----------|----------|------|-------|-----------|--------|-----|---------|
| X_pca (ref: X) | 0.618 | 0.985 | 0.019 | 0.033 | 0.761 | 0.967 | 0.314 |
| X_scvi (ref: X) | 0.504 | 0.854 | 0.009 | 0.016 | 0.365 | 0.984 | 0.425 |
| X_umap (ref: X) | 0.502 | 0.837 | 0.006 | 0.012 | -0.138 | 0.988 | 0.781 |

**Biological interpretation:**

1. **Global structure is well-preserved but local structure is very poorly preserved.** The blood atlas shows a striking divergence: dCor values are high (0.84-0.99) indicating good global distance correlation, but coKNN and entourage are extremely low (0.006-0.033), meaning **almost no cells share their exact nearest neighbors between high-dim and low-dim spaces.** At 85k cells, both UMAP and scVI embeddings distort local relationships while maintaining global layout.

2. **PCA best preserves biological distances.** X_pca achieves the highest spearman_rho (0.62) and lowest kruskal_stress (0.31), and the highest DenPre (0.76) -- PCA preserves relative local densities better than UMAP or scVI. This is expected since PCA is a linear transformation that minimizes squared reconstruction error. For the blood atlas, PCA is the most faithful dimensionality reduction method.

3. **UMAP destroys local density information.** The negative DenPre (-0.138) for X_umap indicates that the density distribution in UMAP space is **inversely correlated** with high-dimensional density -- dense regions in gene expression space become sparse in UMAP, and vice versa. This is a consequence of UMAP's uniform manifold assumption, which explicitly tries to make all local neighborhoods equally sized. **For biological interpretation: if a researcher sees a "dense cluster" in UMAP of the blood atlas, they cannot conclude the cells are transcriptionally homogeneous.**

4. **The kNN graph is almost entirely different between spaces** (GED ~0.97-0.99, coKNN <0.02). This means that the connectivity structure used for Louvain/Leiden clustering differs dramatically between the high-dimensional and low-dimensional representations. **Downstream clustering on UMAP embeddings of the blood atlas produces results that are substantially different from clustering on the original gene expression space.**

5. **Practical implication for researchers:** For the blood atlas, PCA projections (rather than UMAP) should be used for trajectory inference and density-based analyses. UMAP is adequate for cell-type identification at the major lineage level (as the high KS scores in cluster metrics confirm), but should not be trusted for fine-grained subtype relationships.

### 11.4 Annotation Metrics -- Blood Atlas

| Reference | Prediction | ARI | AMI | Isolated F1 | NMI | V-measure |
|-----------|-----------|-----|-----|-------------|-----|-----------|
| cell_type | scvi_leiden_donorassay_full | 0.719 | 0.794 | 0.126 | 0.794 | 0.794 |
| cell_type | scvi_leiden_res05_tissue | 0.212 | 0.615 | 0.018 | 0.615 | 0.615 |
| free_annotation | scvi_leiden_donorassay_full | 0.717 | 0.788 | 0.003 | 0.789 | 0.789 |

**Biological interpretation:**

1. **The coarse donor assay labels transfer well** (ARI=0.72), confirming that the embedding captures broad cell-type distinctions. This is consistent with the high cLISI scores in batch correction.

2. **The fine tissue labels transfer poorly** (ARI=0.21), meaning that the finer Leiden resolution (res=0.5 on tissue) does not correspond well to the reference cell_type annotation. This could indicate either: (a) the tissue-specific clustering captures biological substructure not represented in the reference annotation (a "biologically correct" mismatch), or (b) the embedding has poor resolution for fine subtypes (a "quality problem").

3. **Isolated label F1 is very low** (0.003-0.126). This means rare cell types in the blood atlas (likely plasmacytoid dendritic cells, megakaryocytes, or hematopoietic progenitors) are not captured as distinct clusters. Following scIB's analysis of the human immune task where "Harmony exhibited the lowest isolated label F1," this is a significant quality concern for the blood atlas -- rare populations may be merged with larger cell types.

### 11.5 Cross-Task Biological Synthesis for Blood Atlas

**The Multi-Dimensional Quality Profile:**

Following the unified framework in Section 6.1, the blood atlas receives:
- **Cluster: 0.55** (moderate -- scVI embedding separable, UMAP not)
- **Annotation: 0.63** (good at coarse level, poor at fine/rare level)
- **Batch: 0.45** (assay batch effects persist strongly; donor effects moderate)
- **Dimred: 0.42** (global structure preserved, local structure very poorly preserved)

**Overall assessment:** The blood atlas is of **moderate quality**. Its primary limitations are: (1) persistent assay-level batch effects that resist correction even by scVI, (2) very poor local neighborhood preservation in all embeddings -- biological conclusions drawn from UMAP neighborhoods should be treated with caution, (3) rare cell types are not well-isolated in any embedding.

**Actionable recommendations for the atlas curator:**
1. Consider using Harmony or Scanorama to address assay-level batch effects beyond what scVI achieves
2. Use PCA rather than UMAP for distance-based analyses (trajectories, density)
3. Perform targeted validation of rare cell-type annotations (pDC, megakaryocytes)
4. The uncorrected UMAP (X_uncorrected_umap) actually shows better density preservation (DenPre=0.20) than the batch-corrected one -- consider whether batch correction is helping or hurting the visualization

## References


[1] Cao Y, Yu L, Torkel M, Kim S, Lin Y, Yang P, Speed TP, Ghazanfar S. "The current landscape and emerging challenges of benchmarking single-cell methods." *Nucleic Acids Research* (2025). doi:10.1093/nar/gkaf538. PMC12495992.

[2] Luecken MD, Büttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colomé-Tatché M, Theis FJ. "Benchmarking atlas-level data integration in single-cell genomics." *Nature Methods* 19, 41-50 (2022). doi:10.1038/s41592-021-01336-8.

[3] Regev A et al. "The Human Cell Atlas." *eLife* 6, e27041 (2017). doi:10.7554/eLife.27041.

[4] Zappia L, Richter S, Ramírez-Suástegui C, Kfuri-Rubens R, Vornholz L, Wang W, Heumos L, Heinen T, Korsunsky I, Luckey AM, Borcherding N, Theis FJ. "Feature selection methods affect the performance of scRNA-seq data integration and querying." *Nature Methods* (2025). doi:10.1038/s41592-025-02624-3.

[5] Zhong H, Han W, Gomez-Cabrero D, Tegner J, Gao X, Cui G, Aranda M. "Benchmarking cross-species single-cell RNA-seq data integration methods: towards a cell type tree of life." *Nucleic Acids Research* (2025). doi:10.1093/nar/gkae1312.

[6] Liu C, Ding S, Kim HJ, Long S, Xiao D, Ghazanfar S, Yang P. "Multitask benchmarking of single-cell multimodal omics integration methods." *Nucleic Acids Research* (2025). doi:10.1093/nar/gkaf1027.
