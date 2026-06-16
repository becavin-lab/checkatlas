# scFM QC: foundation-model quality control

CheckAtlas now includes a quality-control layer for single-cell foundation
model (scFM) embeddings. The layer evaluates the embeddings produced by
Geneformer, scGPT, scFoundation, UCE, scCello, and any custom scFM
against the **nine failure modes** documented in
[`scfm_benchmarking_report.md`](scfm_benchmarking_report.md).

The tool is a **fair critic** of the embedding. It does not retrain or
re-implement any scFM. It runs the same CheckAtlas metrics on the
embedding and on user-supplied baselines (PCA, scVI, Harmony, ...), then
translates the raw numbers into a per-problem verdict, an A–F grade, and
a single 0–100 fitness score.

---

## Quick start

```bash
checkatlas scfm path/to/atlas/ \
    --atlas_name pbmc_5k \
    --scfm_embedding X_geneformer \
    --baseline_embeddings X_pca X_scvi \
    --ref_label celltype_author \
    --predicted_label celltype_geneformer \
    --batch_key donor_id \
    --n_seeds 5
```

The scFM embedding is also auto-detected: any `adata.obsm` key whose
name matches a known scFM pattern (`X_geneformer`, `X_scgpt`, `X_uce`,
`X_scfoundation`, `X_sccello`, `X_scbert`, ... or any `X_*scfm` /
`X_*foundation` key) is picked up automatically. The
`--scfm_embedding` flag is a hint, not a requirement.

---

## What you get

Every run produces a small set of files in
`checkatlas_files/scfm/<atlas_name>/`:

| File | What it contains |
|------|------------------|
| `scfm_verdicts.tsv` | one row per (atlas, problem) with score, grade, evidence, and citation |
| `scfm_composite.tsv` | one row per atlas with FMF, BF, PR and per-problem sub-scores |
| `scfm_per_metric.tsv` | the raw metric values, one row per (atlas, metric, embedding) |
| `scfm_inputs.tsv` | the inputs that were used (Becavin comment 6 — input data sanity) |
| `resolved_weights.json` | the composite-score weights actually applied (post-merge) |
| `resolved_thresholds.yaml` | the diagnostic thresholds actually applied (post-merge) |
| `grade_legend.md` | the A–F legend, written to the output folder for self-documentation |

The first two files are also picked up by MultiQC and rendered into the
existing `Checkatlas-MultiQC.html` report alongside the
cluster/annotation/dimred tables.

---

## The nine problems

Each verdict corresponds to one of the failure modes in
`scfm_benchmarking_report.md`:

| ID | Problem | Detection metrics | Source |
|----|---------|-------------------|--------|
| 1  | Zero-shot performance failure | ref-vs-pred F1, ARI | Kedzierska 2024, 2025 |
| 2  | Baselines outperform scFM | delta vs PCA / scVI on every metric | Liu 2024, Souza & Mehta 2026 |
| 3  | Pretraining scaling laws do not apply | `plateau_ratio` from subsampling | DenAdel 2025, Wang 2025 |
| 4  | Batch effects persist | iLISI, kBET, PCR, graph_connectivity | Wang, Zhang & Zhang 2025 |
| 5  | Tokenization and representation failures | rare-vs-common F1, dCor to raw X | Atti & Subramaniam 2025 |
| 6  | Generalization failures across domains | per-pair cross-domain ARI | Souza & Mehta 2026, Boiarsky 2023 |
| 7  | Real-world clinical failures | patient ASW, outcome AUC | Elmarakeby 2025, Csendes 2025 |
| 8  | Model stability issues | metric CV across `n_seeds` random subsamples | Liu 2024, Xu 2026 |
| 9  | The "denoising" hypothesis | scFM ≤ PCA+δ on ≥80% of metrics | Souza & Mehta 2026 |

If a problem cannot be evaluated (e.g. no `batch_key` for problem 4),
the verdict is reported as `n/a` and the explanation tells the user
which input is missing.

---

## Composite scores

Three in-house scores collapse the nine verdicts into a small number a
biologist can act on.

### FoundationModelFitness (FMF, 0–100)

The headline number. A *weighted* geometric mean of seven
"fitness components" derived from the nine verdicts:

```
FMF = 100 * geomean(
    1 - batch_leakage,
    rare_type_recall,
    cross_domain_generalization,
    1 - scaling_plateau_severity,
    1 - instability_severity,
    baseline_competitiveness,
    1 - tokenization_loss,
)
```

Geometric mean (not arithmetic) so any one hard failure floors the
score at 0. A "great batch corrector" that "cannot generalize" is not
fit for production, and the score reflects that.

### BiologicalFaithfulness (BF, 0–1)

```
BF = w1 * cLISI + w2 * graph_connectivity
   + w3 * cell_cycle + w4 * hvg_overlap
   - w5 * iLISI     # penalty for batch
```

`cLISI` low = biology preserved, `iLISI` low = batch leaking. The
subtraction is intentional: the score is a *net* preservation-of-biology
value.

### ProductionReadiness (PR, 0–1)

```
PR = sigmoid(alpha * (1 - instability_CV)
           + beta  * (1 - tokenization_loss)
           + gamma * log(runtime_factor)
           - delta * failure_rate)
```

A sigmoid-based score that captures stability, tokenization fragility,
runtime cost, and the fraction of test datasets where any critical
metric fell below threshold (the "no SOTA across all tasks" finding
from Liu 2024 and Xu 2026).

---

## Grades

Each verdict is mapped to an A–F letter:

| Grade | Score band | Plain language |
|-------|------------|----------------|
| A | 0.00 – 0.15 | No evidence of this problem |
| B | 0.15 – 0.35 | Mild, monitor |
| C | 0.35 – 0.60 | Moderate, address before publication |
| D | 0.60 – 0.80 | Severe, paper claim likely overstated |
| F | 0.80 – 1.00 | Embedding *demonstrates* the problem |
| n/a | n/a | Problem not evaluable on this atlas |

The grade is on the same line as the verdict in `scfm_verdicts.tsv`.

---

## Customising thresholds and weights

Both the diagnostic thresholds (`rules.py` defaults) and the composite
weights (`composite.py` defaults) are user-overridable.

```bash
checkatlas scfm path/to/atlas/ --atlas_name foo \
    --scfm_thresholds /path/to/my_thresholds.yaml \
    --scfm_weights /path/to/my_weights.json
```

A thresholds YAML may contain a partial table (only the problems the
user wants to override). Missing problems fall back to the defaults.
Unknown keys warn and are ignored. The resolved thresholds are written
to `resolved_thresholds.yaml` for reproducibility.

A weights JSON has the same structure as the default
`DEFAULT_WEIGHTS` dict:

```json
{
  "fmf": {"batch": 2.0, "scaling": 0.5},
  "pr":  {"failure_rate": 2.0}
}
```

Missing keys fall back to the defaults. Unknown keys warn and are
ignored. The resolved weights are written to `resolved_weights.json`.

---

## Worked example

The full worked example (Aisha + Geneformer on PBMC) is in §8 of
[`scfm_implementation_plan.md`](scfm_implementation_plan.md). A
shorter version: she runs `checkatlas scfm` with a Geneformer
embedding, gets back `scfm_verdicts.tsv` with a D-grade batch
problem and a citation to Wang, Zhang & Zhang 2025, and a
MultiQC page that renders the verdict alongside her existing
atlas QC.

---

## Implementation notes

- **Becavin comment 1 (auto-detect)**: the column detector's embedding
  patterns include all known scFM names; the user does not need to
  spell out `--scfm_embedding`.

- **Becavin comment 2 (no re-processing)**: every new metric module
  accepts a `preprocess_context` and re-uses the precomputed kNN
  graphs and distance matrices. No re-loading, no re-detection, no
  re-computation.

- **Becavin comment 3 (thresholds in YAML)**: thresholds are
  user-overridable via `--scfm_thresholds`. Defaults are documented
  with their source paper.

- **Becavin comment 4 (MultiQC route 1 first)**: the new tables
  render in the existing MultiQC report via
  `assets/multiqc_config.yml`. The standalone HTML page is a
  post-v1 enhancement.

- **Becavin comment 6 (input data sanity)**: the new
  `scfm_inputs.tsv` records the inputs that were used. The
  existing `qc/` tables are unchanged — atlas-level QC and
  embedding-level QC are deliberately separated.

- **Becavin comment 7 (no thin wrappers)**: only four new metric
  modules were added (`scaling`, `stability`, `rare_types`,
  `cross_domain`); thin re-implementations of existing metrics
  were moved to `composite.py` as score-aggregation helpers.

- **Becavin comment 8 (process, not workflow)**: the scfm step
  is a new process inside `checkatlas_scanpy.nf` (per Becavin's
  instruction "for now don't touch the nextflow integration",
  the actual Nextflow wiring is deferred to a follow-up).

The full architecture, plan, and source code mapping are in
[`scfm_implementation_plan.md`](scfm_implementation_plan.md).
