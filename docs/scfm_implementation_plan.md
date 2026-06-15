# CheckAtlas scFM QC Implementation Plan

> **Status:** Approved design
> **Scope:** Extend CheckAtlas from atlas-level quality control to single-cell foundation model (scFM) embedding quality control.
> **Companion document:** `docs/scfm_benchmarking_report.md` catalogues the nine scFM failure modes this plan targets.
> **Mode:** CheckAtlas is a *fair critic*. It does not re-train or re-implement any scFM. It evaluates the embeddings the user has already produced and reports whether those embeddings exhibit the problems documented in the literature.

---

## 1. Motivation and Framing

The systematic review in `scfm_benchmarking_report.md` identifies **nine recurring failure modes** of current single-cell foundation models (scFMs) such as Geneformer, scGPT, scFoundation, and others. The literature consensus is striking: simpler baselines (PCA, scVI, HVG+logistic regression) match or beat scFMs on most downstream tasks, batch effects persist in scFM embeddings, tokenization destroys information, and the scaling laws that justify foundation models in NLP do not transfer to single-cell biology.

What the literature does *not* provide is a **reusable, runnable tool** that a researcher can use to test their own scFM embedding against these nine problems. Each paper builds an ad-hoc benchmark. The 40+ papers produce results, not infrastructure.

This plan turns CheckAtlas into that infrastructure. Researchers will:

1. Run their scFM (Geneformer, scGPT, scCello, custom) on their data.
2. Store the embedding in `adata.obsm["X_<scfm>"]` and the predicted cell types in `adata.obs["<scfm>_pred"]`.
3. Run `checkatlas scfm …` and receive a verdict on each of the nine problems, a single 0–100 fitness score, and literature-anchored explanations of every problem that was detected.

The tool is a *fair critic* in the same way CheckAtlas is already a fair critic of atlases. It does not solve the scFM problems; it tells the researcher honestly whether the embedding exhibits them.

---

## 2. The Nine Problems (Recap from the Report)

The following table is the master mapping. Every problem is assigned a stable numeric ID, a literature anchor, and the existing or new CheckAtlas metrics that detect it.

| ID | Problem (report) | Detection metrics (existing / new) | Source paper |
|----|-------------------|-------------------------------------|--------------|
| 1  | Zero-shot performance failure | `isolated_f1_score`, `fowlkes_mallows`, `adj_rand_index`, `normalized_mutual_info` (existing); `zero_shot_top3_acc` (new) | Kedzierska 2024, 2025; Boiarsky 2023 |
| 2  | Baselines outperform scFMs | All annotation and clustering metrics run on both scFM embedding and PCA/scVI baselines; the delta is the new `baseline_dominance_score` (composite) | Liu 2024; Souza & Mehta 2026; DenAdel 2025 |
| 3  | Pretraining scaling laws do not apply | New `scaling.py` module: re-evaluate silhouette/ARI/iLISI at multiple subsample fractions; new `plateau_ratio` | DenAdel 2025; Wang 2025 |
| 4  | Batch effects persist | `iLISI`, `kbet`, `pcr`, `graph_connectivity` (existing); new composite `batch_leakage_score` | Wang, Zhang & Zhang 2025 |
| 5  | Tokenization and representation failures | `dCor`, `kruskal_stress`, `spearman_rho` (existing); new `gaussian_robustness_drop`, `rare_type_recall`, `information_loss_ratio` | Atti & Subramaniam 2025; DenAdel 2025; Haber 2025 |
| 6  | Generalization failures across domains | All annotation metrics evaluated per `(train_domain, test_domain)` pair; new `cross_domain_ari_matrix` | Souza & Mehta 2026; Boiarsky 2023; Wenteler 2024 |
| 7  | Real-world clinical/translational failures | `silhouette`, `kbet` on patient labels; new `patient_separation_score`, optional `outcome_AUC` | Elmarakeby 2025; Csendes 2025 |
| 8  | Model stability and interpretability | All metrics re-evaluated with `n_seeds` random subsamples; new `seed_stability_index`, `metric_CV` | Liu 2024; Xu 2026 |
| 9  | The "denoising" hypothesis | Comparison of scFM vs PCA across the full metric set; new composite `denoiser_score` | Souza & Mehta 2026 |

The crucial insight is that **problems 1, 2, 4, 6, 7 are largely covered by metrics CheckAtlas already has**. The new work is concentrated on:

- New metric modules for problems 3, 5, 7, 8
- A new diagnostic layer that interprets the raw numbers as one of the nine problems
- A new narrative layer with composite scores
- A new CLI process and report writers

---

## 3. Architecture Overview

```
┌──────────────────────────────────────────────────────────────────────┐
│  Layer 3 — BIOLOGIST NARRATIVE                                       │
│  scfm_report.py   scfm_grade()   MultiQC HTML page                   │
├──────────────────────────────────────────────────────────────────────┤
│  Layer 2 — INTERPRETATION & DIAGNOSIS                                │
│  checkatlas/scfm/   diagnostics.py   rules.py   composite.py         │
├──────────────────────────────────────────────────────────────────────┤
│  Layer 1 — METRIC COLLECTION (re-uses existing code)                 │
│  cal_cluster()   cal_annot()   cal_dimred()   + new scfm modules     │
├──────────────────────────────────────────────────────────────────────┤
│  EXISTING FOUNDATION (unchanged)                                      │
│  atlas.py   col_detector.py   PreprocessContext   kNN cache   JAX    │
└──────────────────────────────────────────────────────────────────────┘
```

**Layer 1** collects raw metric values. It re-uses ~80% of CheckAtlas code and adds a small set of new metric modules.

**Layer 2** is the *new* intellectual contribution. It takes the long-format metric table from Layer 1 and applies nine rule-based verdicts, each anchored in a specific paper from the report. Each verdict carries a 0–1 score, a confidence band, a list of evidence metrics, an explanation string, and a citation.

**Layer 3** is what the user sees. It collapses the nine verdicts into three in-house scores and an A–F grade per problem, and writes MultiQC-friendly tables.

---

## 4. CLI and Configuration

### 4.1 New process

Add `"scfm"` to `PROCESS_TYPE` in `checkatlas/check.py:35` so it is a first-class subcommand alongside `metric_cluster`, `metric_annot`, and `metric_dimred`. The dispatch lives in `checkatlas/__main__.py:142-178` and is straightforward to extend.

### 4.2 New arguments (in `checkatlas/utils/checkatlas_arguments.py`)

A new argparse group, all optional except the obvious ones, forwarded to Nextflow via the existing `_metric_nf_arg` pattern.

```
--scfm_embedding NAME              required, e.g. X_geneformer
--baseline_embeddings NAME [NAME …]  e.g. X_pca X_scvi X_harmony
--scfm_predicted_label NAME        the scFM's predicted cell types
--scfm_batch_key NAME              e.g. donor_id
--scfm_ref_label NAME              e.g. celltype_author
--scfm_scaling_fractions F F F     default "0.01 0.10 0.50 1.00"
--scfm_n_seeds N                   default 5
--scfm_noise_sigma F               default 0.10
--scfm_min_domain_size N           default 50
--scfm_outdir DIR                  default checkatlas_files/scfm/
--scfm_weights PATH                optional JSON file overriding composite weights
```

### 4.3 Output structure

New folder registered in `checkatlas/utils/folders.py:22`:

```
checkatlas_files/scfm/
├── scfm_verdicts.tsv        long format, one row per (atlas, problem)
├── scfm_composite.tsv       one row per atlas, columns for FMF, BF, PR
├── scfm_per_metric.tsv      the raw metric values, one row per (atlas, metric, embedding, label)
├── scfm_problem_heatmaps/   per-atlas bar plots of the 9 problem scores
└── scfm_report.html         standalone HTML summary page
```

### 4.4 Example invocation

```bash
checkatlas scfm path/to/pbmc/ \
    --atlas_name pbmc_5k \
    --scfm_embedding X_geneformer \
    --baseline_embeddings X_pca X_scvi \
    --ref_label celltype_author \
    --predicted_label celltype_geneformer \
    --batch_key donor_id \
    --n_seeds 5
```

---

## 5. Layer 1 — Metric Modules

### 5.1 Re-use of existing code

The new `cal_scfm(...)` orchestrator (in `checkatlas/metrics/scfm/run.py`) is a thin wrapper that calls the existing pipelines with augmented inputs:

- `cal_annot(adata, metric_list=extended_list, ...)` where `extended_list` includes LISI, kBET, PCR, graph_connectivity in addition to ref-vs-pred metrics.
- `cal_cluster(adata, metric_list=cluster.__all__, ...)` with all embeddings.
- `cal_dimred(adata, metric_list=dimred.__all__, ...)` with all embeddings.

The scFM embedding is added to the auto-detected `embedding_keys` list, so every metric that supports embeddings will be evaluated on it. The user's baseline embeddings (PCA, scVI) are added to the same list. This is the key design choice: **scFM and baselines are evaluated on the exact same metric engine, so any comparison is apples-to-apples**.

##### Review Becavin: Lets add in our embedding detection every embedding from the scFM found (scGPT, GeneFormer, UCE, etc;..)

### 5.2 New modules in `checkatlas/metrics/scfm/`

```
scfm/
├── __init__.py
├── run.py                 # orchestrator: cal_scfm(adata, config) -> pd.DataFrame
├── scaling.py             # problem 3
├── stability.py           # problem 8
├── information_loss.py    # problem 5
├── rare_types.py          # problem 5 (overprediction of common types)
└── cross_domain.py        # problem 6
```

##### Review Becavin: Please if you add this it has to be connected to metrics.py; It should be ran after cla_cluster and other metrics. No re-processing of adata should be performed

`__init__.py` exposes:

```python
__all__ = ["scaling", "stability", "information_loss", "rare_types", "cross_domain"]
```

Each module follows the same convention as existing metric modules: a `run(...)` function with consistent signature, returning a float (or a small DataFrame for `cross_domain.py` and `scaling.py`). Every result is appended to the same long-format DataFrame so Layer 2 sees a uniform input.

#### `scaling.py` — Problem 3

```python
def run_scaling(adata, scfm_emb, ref_label, batch_key,
                fractions=(0.01, 0.05, 0.10, 0.50, 1.00),
                seed=42,
                metric_subset=("silhouette", "ari", "ilisi")):
    """For each fraction f, subsample cells, recompute metric on scfm_emb.
    Returns long DataFrame and derived 'plateau_ratio' = metric[1.00] / metric[0.01].
    """
```

A `plateau_ratio < 1.05` reproduces the DenAdel et al. (2025) finding that performance saturates at ~1% of pretraining data. The function re-uses existing silhouette and ARI implementations; iLISI is computed via the existing `lisi.run_with_neighbors` path. The result is one row per `(fraction, metric, embedding)` triple.

#### `stability.py` — Problem 8

```python
def run_stability(metrics_to_test, adata, scfm_emb, ref_label, batch_key,
                  n_seeds=5, base_seed=0):
    """For each seed, subsample 90% of cells with random_state=seed,
    recompute the requested metrics, return mean, std, CV per metric.
    """
```

Output is one row per metric with columns `mean`, `std`, `CV`. A high `CV` (e.g. `> 0.10`) across multiple metrics matches the Liu 2024 / Xu 2026 stability finding.

#### `information_loss.py` — Problem 5

```python
def gaussian_robustness(adata, scfm_emb, ref_label, pred_label,
                        noise_sigma=0.1, n_repeats=3, seed=42):
    """For each repeat, add N(0, sigma * std) noise to raw X,
    re-embed (or reuse the existing embedding) and recompute ARI.
    Report (ARI_clean - ARI_noisy) / ARI_clean.
    High drop = embedding fragile to input noise (tokenization problem).
    """

def information_loss_ratio(adata, scfm_emb, ref_label, n_samples=5000):
    """Compression-fidelity: ratio of (dCor scFM -> ref X) to
    (dCor PCA -> ref X). Lower means scFM loses more information
    than the simpler PCA projection.
    """
```

#### `rare_types.py` — Problem 5

```python
def rare_type_recall(adata, scfm_emb, pred_label, ref_label,
                     rare_quantile=0.10):
    """Cells in the bottom rare_quantile of class frequency are 'rare'.
    Per-class F1 on rare cells. If rare F1 << common F1, the model
    overpredicts common cell types (scGPT finding on CD14).
    """
```

#### `cross_domain.py` — Problem 6

```python
def cross_domain_ari(adata, scfm_emb, ref_label, domain_key,
                     min_domain_size=50, n_perm=10):
    """For each (train_domain, test_domain) pair with size >= min,
    project test cells into scfm_emb, compute ARI of test clustering
    vs test labels. Return domain x domain matrix plus summary
    statistics (mean off-diagonal, fraction below threshold).
    """
```

This generalises naturally to species transfer (the Souza & Mehta 2026 F1 < 0.5 finding) by passing `domain_key="species"`. It also works on `tissue`, `assay`, or any other categorical `obs` column.

### 5.3 Patient metrics (problem 7) — added in `checkatlas/metrics/scfm/patient.py`

```python
def patient_separation_score(adata, scfm_emb, patient_key,
                             outcome_key=None, n_jobs=-1):
    """ASW on patient labels (lower = patients better mixed; desirable
    for batch correction; higher = patients well separated; desirable
    for outcome prediction). Returns (asw, n_patients, optional outcome_AUC).
    Interpretation direction is recorded in a 'direction' field so
    Layer 2 knows which way to read the value.
    """
```

This is a thin wrapper around the existing `silhouette.run` and a small logistic-regression with cross-validation when `outcome_key` is provided.

---

## 6. Layer 2 — Diagnostics Engine

New directory `checkatlas/scfm/`. This is the layer that translates numbers into problems.

```
scfm/
├── __init__.py
├── config.py             # SCFMConfig dataclass: embeddings, baselines, batch keys, seeds
├── run.py                # entry: orchestrates Layer-1 metric collection
├── diagnostics.py        # per-problem diagnostic verdicts
├── rules.py              # threshold table per problem (literature-grounded)
├── composite.py          # the in-house aggregated scores
├── report.py             # MultiQC-friendly output writers
└── grades.py             # letter grade A–F per problem
```

### 6.1 `config.py` — `SCFMConfig`

A frozen dataclass that captures everything Layer 1 and Layer 2 need to know about a run. Constructed in `__main__.py` from `argparse.Namespace`.

```python
@dataclass(frozen=True)
class SCFMConfig:
    atlas_name: str
    scfm_embedding: str
    baseline_embeddings: tuple[str, ...]
    predicted_label: str | None
    ref_label: str | None
    batch_key: str | None
    domain_key: str | None
    patient_key: str | None
    outcome_key: str | None
    scaling_fractions: tuple[float, ...]
    n_seeds: int
    noise_sigma: float
    min_domain_size: int
    weights_path: str | None
```

### 6.2 `rules.py` — literature threshold table

A single dict (or YAML, see §11) mapping each problem ID to:

```python
{
    4: {
        "name": "Batch effects persist",
        "metric_thresholds": [
            ("kbet", ">", 0.50, "high"),
            ("ilisi", "<", 0.30, "high"),
            ("pcr", ">", 0.85, "medium"),
            ("graph_connectivity", "<", 0.50, "low"),
        ],
        "primary_reference": "Wang, Zhang & Zhang 2025",
    },
    ...
}
```
##### Review Becavin: also the user should be able to put its own thresholds in a modified yaml file
Threshold source policy: **hard-coded with cited sources in `rules.py`**. Each entry records the paper and (where available) the threshold the paper used. Threshold updates are a one-line change.

### 6.3 `diagnostics.py` — `ProblemVerdict` and the nine rules

```python
@dataclass
class ProblemVerdict:
    problem_id: int          # 1..9
    problem_name: str
    score: float             # 0..1, higher = more problematic
    confidence: str          # "high" | "medium" | "low" | "n/a"
    evidence: list[str]      # metric names that drove the verdict
    explanation: str         # human-readable sentence
    reference: str           # paper citation
```

`diagnose(metrics_df: pd.DataFrame, config: SCFMConfig) -> list[ProblemVerdict]`

The nine rule implementations:

| ID | Rule (informal) | Triggering metrics | Threshold source |
|----|------------------|---------------------|-------------------|
| 1 | Zero-shot failure | `isolated_f1` < 0.5 on at least one (ref, pred) pair | Kedzierska 2024 |
| 2 | Baseline dominates | For every metric, scFM ≤ PCA+δ | Liu 2024, Souza 2026 |
| 3 | Scaling plateau | `plateau_ratio` < 1.10 across 1% → 100% | DenAdel 2025 |
| 4 | Batch leakage | `kbet_rejection > 0.5` OR `iLISI < 0.3` OR `pcr > 0.85` | Wang 2025 |
| 5 | Tokenization loss | `gaussian_robustness_drop > 0.10` OR `rare_type_recall < common_F1 - 0.15` OR `dCor < 0.5` | Atti 2025 |
| 6 | No cross-domain | For ≥50% of (train, test) domain pairs `ARI < 0.3` | Souza 2026 |
| 7 | Patient collapse | `asw_patient < 0.05` AND `outcome_AUC < 0.6` | Elmarakeby 2025 |
| 8 | Instability | `CV > 0.10` on ≥ 2 metrics across seeds | Liu 2024, Xu 2026 |
| 9 | Denoiser only | scFM ≤ PCA+δ on ≥ 80% of metrics (overlap with #2, kept for narrative) | Souza 2026 |

**Missing-data policy:** if a required input is absent (e.g. no `batch_key`, no rare types, single domain only), the verdict is reported as `confidence="n/a"` with `evidence=[]` and an explanation that the problem could not be evaluated on this atlas. This matches the existing CheckAtlas behaviour where individual metrics warn-and-skip. The user always sees which problems were evaluable and which were not.

Each rule has a **low/medium/high confidence** band based on the magnitude of deviation:

- `high`: at least one threshold breached by a factor of ≥ 2
- `medium`: breach by 1.0–2.0×
- `low`: marginal breach (within 20% of threshold)

### 6.4 `composite.py` — three in-house scores

Three scores collapse the 30+ raw metrics and 9 verdicts into a small number a biologist can act on.

#### `FoundationModelFitness` (FMF) — overall, 0..100

The headline number. Geometric mean (not arithmetic) so any one failure pulls the whole score down.

```
FMF = 100 * geometric_mean(
    1 - batch_leakage_normalized,
    rare_type_recall,
    cross_domain_generalization,
    1 - scaling_plateau_severity,
    1 - instability_severity,
    baseline_competitiveness,
    gaussian_robustness,
)
```

Each input is normalised to [0, 1] using literature-derived scales. Geometric mean is the right aggregator: a "great batch corrector" that "cannot generalize" is *not* fit for production. This is the exact intuition biologists already have when they read a paper.

#### `BiologicalFaithfulness` (BF) — does the embedding still hold biology?

```
BF = w1 * cLISI_norm + w2 * graph_connectivity_norm
   + w3 * cell_cycle_conservation + w4 * HVG_overlap
   - w5 * iLISI_norm   (penalty for batch)
```

`cLISI` low = biology preserved, `iLISI` low = batch leaking. Subtraction is intentional: the metric is a *net* preservation-of-biology score.

#### `ProductionReadiness` (PR) — should you ship this?

```
PR = sigmoid(
    alpha * (1 - instability_CV) +
    beta  * (1 - tokenization_loss) +
    gamma * runtime_factor +
    delta * (1 - failure_rate_across_datasets)
)
```

`runtime_factor` is `1 / (median_runtime_minutes / budget)`. `failure_rate_across_datasets` is the fraction of test datasets where any critical metric fell below threshold — directly the "no SOTA across all tasks" finding from Liu 2024 and Xu 2026.

**Default weights:** hard-coded in `composite.py` as a module-level dict. A user-supplied JSON file at `--scfm_weights` overrides individual values; missing keys fall back to defaults. The JSON is also written to the output folder so every report is reproducible.

```python
DEFAULT_WEIGHTS = {
    "fmf": {
        "batch": 1.0,
        "rare_types": 1.0,
        "cross_domain": 1.0,
        "scaling": 1.0,
        "stability": 1.0,
        "baseline_competitiveness": 1.0,
        "gaussian_robustness": 1.0,
    },
    "bf": {"cLISI": 1.0, "graph_connectivity": 1.0,
           "cell_cycle": 0.5, "hvg_overlap": 0.5, "iLISI_penalty": 1.0},
    "pr": {"stability": 1.0, "tokenization": 1.0,
           "runtime": 0.5, "failure_rate": 1.5},
}
```

### 6.5 `grades.py` — A–F per problem

```python
def letter(score: float) -> str:
    if score < 0.15: return "A"   # No evidence of this problem
    if score < 0.35: return "B"   # Mild, monitor
    if score < 0.60: return "C"   # Moderate, address before publication
    if score < 0.80: return "D"   # Severe, paper claim likely overstated
    return "F"                    # Embedding demonstrates the problem
```

Mapping is also written to a small `grade_legend.md` in the output folder for self-documentation.

---

## 7. Layer 3 — Reporting and Output

### 7.1 `report.py` — file writers

Three writers, all idempotent:

```python
def write_verdicts(verdicts: list[ProblemVerdict], outdir: str) -> Path: ...
def write_composite(atlas_name: str, scores: dict, outdir: str) -> Path: ...
def write_per_metric(metrics_df: pd.DataFrame, outdir: str) -> Path: ...
```

`scfm_verdicts.tsv` columns: `atlas_name, problem_id, problem_name, score, grade, confidence, evidence, explanation, reference`.

`scfm_composite.tsv` columns: `atlas_name, fmf, bf, pr, batch_score, scaling_score, stability_score, ...`.

`scfm_per_metric.tsv` is the raw layer-1 output, one row per `(atlas, metric, embedding, label)`.

### 7.2 MultiQC integration

Two routes, both supported:

1. **Custom-content table** in the existing MultiQC report. `assets/multiqc_config.yml` gets a new `scfm_verdicts` and `scfm_composite` table definition, mirroring the existing `cluster`, `annotation`, `dimred` tables.
2. **Standalone HTML report** at `checkatlas_files/scfm/scfm_report.html`, generated by a new `generate_scfm_html(path, atlas_name)` function that mirrors the existing `check.generate_fig_html`. The page shows:
   - Headline: `FMF = NN/100, grade X`
   - Nine-problem bar plot (matplotlib, embedded as base64 PNG)
   - The three composite scores with one-line explanations
   - Per-problem expandable sections with evidence, explanation, and reference
##### Review Becavin: Start with route 1, and then i will help implement the second route, this is easy later.


A small `scfm_problem_heatmaps/` subfolder contains the per-atlas bar plots, so users who prefer the static image can also browse them.

And then i will help implement the second route, this is easy

### 7.3 Nextflow integration

`checkatlas/nextflow/main.nf:26-32` is extended to include a new `CHECKATLAS_SCFM` workflow. It runs the SCFM process on every atlas in the samplesheet, after the standard metric processes. The output TSVs are collected and forwarded to `CREATE_REPORT` and `MULTIQC` exactly as the other tables are.
##### Review Becavin: Not a new workflow but a new process as below invoking checkatlas scfm <atlas_path> --atlas_name <name> --scfm_embedding <emb> ...`. after all metrics are done.

The minimal change to `checkatlas/nextflow/workflows/checkatlas_scanpy.nf` is a new process module that invokes `checkatlas scfm <atlas_path> --atlas_name <name> --scfm_embedding <emb> ...`. The Nextflow wrapper shells out to the CLI subprocess just as the existing metric processes do.

---

## 8. Worked Example

Researcher *Aisha* has fine-tuned Geneformer on her PBMC dataset. She stores:

- `adata.obsm["X_geneformer"]` — the embedding
- `adata.obs["celltype_geneformer"]` — predicted cell types
- `adata.obs["celltype_author"]` — author-curated labels (the ground truth)
- `adata.obs["donor_id"]` — batch labels
- `adata.obsm["X_pca"]` and `adata.obsm["X_scvi"]` — her baseline embeddings
- `adata.uns["scfm_config"]` — JSON dump of the arguments she will pass

She runs:

```bash
checkatlas scfm path/to/pbmc/ \
    --atlas_name pbmc_5k \
    --scfm_embedding X_geneformer \
    --baseline_embeddings X_pca X_scvi \
    --ref_label celltype_author \
    --predicted_label celltype_geneformer \
    --batch_key donor_id \
    --n_seeds 5
```

What happens, in order:

1. **Preprocessing** (`atlas.preprocess_atlas`). All three embeddings (Geneformer, PCA, scVI) are fingerprinted. The column detector identifies `celltype_author` as a reference and `celltype_geneformer` as predicted. kNN graphs and distance matrices are precomputed once for the largest of the embeddings and reused.

2. **Layer 1 metric collection**. `cal_annot`, `cal_cluster`, and `cal_dimred` are called with all three embeddings and all metrics. Then `metrics.scfm.run.cal_scfm(...)` runs the new modules:
   - `scaling.run_scaling` on `X_geneformer` and `X_pca` at fractions `(0.01, 0.05, 0.10, 0.50, 1.00)`.
   - `stability.run_stability` on `X_geneformer` with `n_seeds=5`.
   - `information_loss.gaussian_robustness` with `noise_sigma=0.10`, `n_repeats=3`.
   - `rare_types.rare_type_recall` with `rare_quantile=0.10`.
   - `cross_domain` is skipped (no `domain_key` provided) and flagged as `n/a`.

3. **Layer 2 diagnostics**. `diagnostics.diagnose(df, config)` returns nine `ProblemVerdict` objects. For Aisha's Geneformer-on-PBMC run, the verdicts are:

##### Review Becavin: Add a genomic diagnostic and sequencing = normal QC on the quality of the single cell. Is there enough gene per cells, percentage of mito... It is already implemented in the QC process.


   | # | Problem | Grade | Evidence |
   |---|---------|-------|----------|
   | 1 | Zero-shot failure | C | `isolated_f1 = 0.42` on `celltype_geneformer` vs `celltype_author` |
   | 2 | Baseline dominates | B | `silhouette_GF = 0.18`, `silhouette_PCA = 0.16`, `silhouette_scVI = 0.21` — Geneformer is between baselines |
   | 3 | Scaling plateau | A | `plateau_ratio = 1.04` |
   | 4 | Batch leakage | **D** | `kBET = 0.62`, `iLISI = 0.21`, `PCR = 0.88` |
   | 5 | Tokenization | C | `gaussian_drop = 0.13`, `rare_F1 = 0.31` vs `common_F1 = 0.58` |
   | 6 | Cross-domain | n/a | only one tissue |
   | 7 | Patient | n/a | no patient IDs |
   | 8 | Instability | A | `CV = 0.04` across 5 seeds |
   | 9 | Denoiser | B | scFM marginally above PCA on 60% of metrics |

4. **Layer 3 composite**. The three scores are computed from the verdicts:
   - `FMF = 100 * geomean(0.38, 0.31, 1.0, 0.92, 0.96, 0.65, 0.87) ≈ 67`
   - `BF = 0.50` (cLISI good, graph_connectivity weak, iLISI penalty strong)
   - `PR = 0.72` (stable, no runtime data, 0 failure rate on a single dataset)

5. **Reporting**. `scfm_verdicts.tsv`, `scfm_composite.tsv`, `scfm_per_metric.tsv` are written. A bar plot of the 9 problem scores and a standalone HTML report are produced. MultiQC picks up the new tables.

Aisha's screen shows:

- **Headline**: "FMF = 67/100, grade C. Your embedding has a severe batch leakage (problem 4, grade D) and moderate tokenization fragility (problem 5, grade C). Evidence trail: kBET 0.62, iLISI 0.21 (Wang, Zhang & Zhang 2025); gaussian_robustness_drop 0.13 (Atti & Subramaniam 2025). Consider Harmony post-correction or scVI to fix problem 4."
- **Nine-problem bar plot** with each problem's grade and the underlying metric values.
- **Actionable suggestions** per problem, e.g. for problem 4: "Your scFM embedding behaves like the 8 scFMs Wang et al. (2025) evaluated. Their recommendation: do not rely on the scFM embedding alone for batch integration; apply Harmony or scVI as a post-correction step."

Aisha does not have to be a benchmarking expert to know what to do.

---

## 9. Implementation Phases

This is roughly **8–10 weeks of focused work for a single engineer**, organised so the package is always runnable. Every phase ends with a green test suite and a CHANGELOG entry.

### Phase 0 — Hooks (1 week)

- Add `"scfm"` to `PROCESS_TYPE` in `check.py:35`.
- New argparse group in `checkatlas_arguments.py`.
- Stub `cal_scfm` in `metrics/scfm/__init__.py`.
- Register `checkatlas_files/scfm/` in `folders.py:22`.
- Update `__main__.py:142-178` to dispatch `scfm` to a new top-level handler.

**Deliverable:** `checkatlas scfm path/ --atlas_name foo` returns "scfm pipeline not yet implemented" without crashing.

### Phase 1 — Layer 1 new metrics (3 weeks)

##### Review Becavin: New metrics or just the old on reimplemented ??? Be careful ! 

- `metrics/scfm/scaling.py`
- `metrics/scfm/stability.py`
- `metrics/scfm/information_loss.py`
- `metrics/scfm/rare_types.py`
- `metrics/scfm/cross_domain.py`
- `metrics/scfm/patient.py`
- `metrics/scfm/run.py` orchestrator
- Unit tests mirroring the existing `tests/test_metrics.py` style (the file is currently empty — it gets filled).

**Deliverable:** `from checkatlas.metrics.scfm import run; df = run.cal_scfm(adata, config)` returns a long-format DataFrame. Tests on synthetic atlases engineered to trigger each rule.

### Phase 2 — Layer 2 diagnostics (2 weeks)

- `scfm/config.py` (`SCFMConfig` dataclass)
- `scfm/rules.py` (literature threshold table)
- `scfm/diagnostics.py` (nine `ProblemVerdict` rules)
- Tests with synthetic atlases engineered to trigger each rule.

**Deliverable:** `from checkatlas.scfm.diagnostics import diagnose; verdicts = diagnose(df, config)` returns nine `ProblemVerdict` objects. Each rule has a test fixture.

### Phase 3 — Layer 3 composite + grades (1 week)

- `scfm/composite.py` — FMF, BF, PR with default weights and JSON override.
- `scfm/grades.py` — A–F mapping.
- Sensitivity tests (change one input, confirm the score moves in the right direction).

**Deliverable:** `from checkatlas.scfm.composite import fmf, bf, pr` and `from checkatlas.scfm.grades import letter`.

### Phase 4 — Reporting (1 week)

- `scfm/report.py` — TSV writers.
- `scfm/html.py` — `generate_scfm_html(path, atlas_name)` mirroring `check.generate_fig_html`.
- Extend `assets/multiqc_config.yml` with `scfm_verdicts` and `scfm_composite` table definitions.
- Update `checkatlas/nextflow/main.nf` to include `CHECKATLAS_SCFM` workflow.
- Update `checkatlas/nextflow/workflows/checkatlas_scanpy.nf` with a new process module.
- Update MultiQC to pick up new tables.

**Deliverable:** Running the full pipeline produces all four output files and renders correctly in MultiQC.

### Phase 5 — Docs and example (1 week)

- New page `docs/scfm.md` with the full architecture, the worked example, and the rule table.
- Worked example on a real dataset from `tests/data/`.
- Update `README.md` "Summary" bullet list to mention the new process.
- Update `docs/usage.md` to show the new `checkatlas scfm` invocation.

**Deliverable:** `mkdocs build` succeeds, new page renders, README updated.

---

## 10. Testing Strategy

Three layers of tests, mirroring the architecture:

1. **Unit tests** for each metric module in `metrics/scfm/`. Each function gets at least three synthetic atlases: one designed to trigger the problem, one designed to pass, one with the required input missing.
2. **Rule tests** in `tests/test_scfm_diagnostics.py`. One test per problem with a synthetic atlas engineered to fail the rule.
3. **End-to-end test** in `tests/test_scfm_pipeline.py`. Builds a small AnnData with a Geneformer-like embedding, runs the full `cal_scfm → diagnose → composite → report` chain, and checks the output TSVs and HTML report exist with the expected rows.

The existing `tests/test_metrics.py` is empty. We will fill it as part of Phase 1 — not just for the new SCFM metrics but also for the existing ones that currently have no coverage.

---

## 11. Configuration and Threshold Policy

### Composite weights

- **Default location:** hard-coded dict `DEFAULT_WEIGHTS` in `checkatlas/scfm/composite.py`.
- **Override path:** `--scfm_weights /path/to/weights.json`. The JSON has the same structure as `DEFAULT_WEIGHTS`; missing keys fall back to defaults. Unknown keys warn and are ignored.
- **Persisted copy:** the resolved weights (after merge) are written to `checkatlas_files/scfm/resolved_weights.json` so every report is reproducible.

### Threshold table

- **Default location:** hard-coded dict in `checkatlas/scfm/rules.py`.
- **Source policy:** every threshold entry records the paper it came from and (where available) the threshold the paper used. Updates are a one-line change.
- **No runtime override:** thresholds are not user-tunable in v1. If a user disagrees with a threshold, they are encouraged to file an issue with the paper citation that motivates a change. This keeps the tool from drifting into "everyone uses their own thresholds" territory.

---

## 12. What Makes This Impactful for Biologists

Six concrete deliverables turn CheckAtlas from "another metric dump" into something biologists actually use:

1. **A single headline number per scFM run** (FMF, 0–100). Fits in a paper abstract or a Slack message.
2. **A nine-problem grade table** (A–F). The report's findings become a checklist the researcher literally checks off.
3. **Literature-anchored explanations**. Every problem cites the paper that identified it. The verdict is reproducible: another group running the same data gets the same grade. This is what was missing from the 40+ papers in the report — none of them provided reusable tooling.
4. **Actionable suggestions** written into the verdict, not buried in a paper. "Problem 4 detected → consider Harmony / scVI post-correction."
5. **Baseline-relative framing**. The user always sees their scFM *next to* PCA, scVI, Harmony on the same metrics. The Kedzierska / Souza / DenAdel finding becomes a visible chart, not a one-line conclusion.
6. **Trust through transparency**. The composite scores are sums/geometric-means of the underlying metrics, all visible in `scfm_per_metric.tsv`. A reviewer can disagree with the aggregation, audit it, and propose alternative weights via the JSON override.

The biologist does not need to read 40 papers. They run `checkatlas scfm …`, see FMF = 67 with a D-grade batch problem and a citation to Wang 2025, and they have an evidence trail that scales to a manuscript.

---

## 13. Open Questions and Future Work

Items deferred to a future iteration:

- **Auto-calibration of thresholds**: a small reference panel of known-good (PCA on clean data) and known-bad (Geneformer on cross-species data) embeddings would let us auto-derive thresholds rather than hand-coding them. Significant additional work; not in v1.
- **Wrapper APIs for common scFMs**: Geneformer, scGPT, scCello all have Python APIs. Wrapping them would let `checkatlas scfm` run the model itself rather than expecting the user to do it. Breaks if upstream code becomes unmaintained. Out of scope.
- **Multi-dataset aggregation**: the current `cal_scfm` runs on one atlas at a time. A `checkatlas scfm_aggregate /path/to/reports/` would average FMF across a benchmark suite (e.g. the 5 datasets from Kedzierska 2024). This is a natural follow-up once the single-atlas pipeline is stable.
- **Active learning loop**: identify the *one* experiment (which metric, which subsample size) that would most improve the FMF confidence, and suggest it to the user. Powerful but speculative; deferred to v2.

---

## 14. Summary of New Code

| Path | Purpose | Status |
|------|---------|--------|
##### Review Becavin: Not sure it is usefull
| `checkatlas/metrics/scfm/__init__.py` | Package init, exports | New |
| `checkatlas/metrics/scfm/run.py` | `cal_scfm` orchestrator | New |
| `checkatlas/metrics/scfm/scaling.py` | Problem 3 | New |
| `checkatlas/metrics/scfm/stability.py` | Problem 8 | New |
| `checkatlas/metrics/scfm/information_loss.py` | Problem 5 | New |
| `checkatlas/metrics/scfm/rare_types.py` | Problem 5 | New |
| `checkatlas/metrics/scfm/cross_domain.py` | Problem 6 | New |
| `checkatlas/metrics/scfm/patient.py` | Problem 7 | New |
#####

| `checkatlas/scfm/__init__.py` | Package init, exports | New |
| `checkatlas/scfm/config.py` | `SCFMConfig` dataclass | New |
| `checkatlas/scfm/rules.py` | Literature threshold table | New |
| `checkatlas/scfm/diagnostics.py` | Nine `ProblemVerdict` rules | New |
| `checkatlas/scfm/composite.py` | FMF, BF, PR | New |
| `checkatlas/scfm/grades.py` | A–F mapping | New |
| `checkatlas/scfm/report.py` | TSV writers | New |
| `checkatlas/scfm/html.py` | Standalone HTML report | New |
| `checkatlas/check.py` | Add `"scfm"` to `PROCESS_TYPE` | Edit |
| `checkatlas/utils/checkatlas_arguments.py` | New argparse group | Edit |
| `checkatlas/utils/folders.py` | Register `scfm/` folder | Edit |
| `checkatlas/__main__.py` | Dispatch `scfm` process | Edit |

##### Review Becavin: No new workflow but a process in checkatlas_scanpy.nf
| `checkatlas/nextflow/main.nf` | Include `CHECKATLAS_SCFM` workflow | Edit |
| `checkatlas/nextflow/workflows/checkatlas_scanpy.nf` | New process module | Edit |
#####

| `checkatlas/nextflow/assets/multiqc_config.yml` | New table definitions | Edit |
| `tests/test_metrics.py` | Unit tests for new metric modules | Fill in |
| `tests/test_scfm_diagnostics.py` | One test per problem rule | New |
| `tests/test_scfm_pipeline.py` | End-to-end test | New |
| `docs/scfm.md` | Documentation page | New |
| `README.md` | Update Summary section | Edit |
| `docs/usage.md` | Add new `checkatlas scfm` example | Edit |

Roughly 17 new files, 8 edited files, ~2,500–3,000 lines of Python (including tests and docs).

---

*End of plan.*
