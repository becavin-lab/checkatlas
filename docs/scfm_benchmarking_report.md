# Single-Cell Foundation Models: A Critical Review of Benchmarking Literature

> **Generated using Paperclip** -- systematic literature search across 11M+ papers  
> **Date:** June 2026  
> **Focus:** Identifying problems, limitations, and evidence that current scRNA-seq foundation model claims are overstated

---

## Executive Summary

A systematic review of 40+ benchmarking papers (2023-2026) reveals a striking consensus: **current single-cell foundation models (scFMs) consistently fail to deliver on their foundational promise**. Across multiple independent evaluations, simpler methods (PCA, logistic regression, scVI, HVG selection) match or outperform models like Geneformer, scGPT, scFoundation, and others on the majority of downstream tasks. The evidence points to fundamental problems in pretraining strategies, tokenization paradigms, and the core assumption that scaling laws from NLP transfer to single-cell biology. Below we catalog the problems identified by the literature, organized by category.

---

## 1. Zero-Shot Performance Failure

**Problem:** scFMs do not work reliably without fine-tuning, contradicting the "foundation model" premise.

**Key Papers:**
- Kedzierska et al. (2024) *Genome Biology* -- "Assessing the limits of zero-shot foundation models in single-cell biology"
- Kedzierska et al. (2025) *Genome Biology* -- "Zero-shot evaluation reveals limitations of single-cell foundation models"

**Evidence:**
- Geneformer and scGPT exhibit **limited reliability** in zero-shot settings and often underperform compared to simpler methods
- In cell type clustering, HVG selection outperforms Geneformer in 4 out of 5 datasets
- In batch integration, the best scores across all datasets were achieved by simply selecting highly variable genes (HVG)
- Matching the pretraining tissue to the target task does **not** guarantee performance above random initialization
- Increasing pretraining dataset size and diversity **sometimes decreases** zero-shot performance
- The pretraining objectives (masked gene modeling) do not provide meaningful information for biological applications; mean expression prediction beats scGPT reconstructions

> "Our results indicate that both Geneformer and scGPT exhibit limited reliability in zero-shot settings and often underperform compared to simpler methods." -- Kedzierska et al. (2024)

---

## 2. Baselines Outperform: Simpler Methods Win

**Problem:** Across the majority of benchmarks, classical methods (PCA, scVI, logistic regression, HVG+linear) match or exceed foundation model performance.

**Key Papers:**
- Liu et al. (2024) *Advanced Science* -- "Evaluating the Utilities of Foundation Models in Single-cell Data Analysis"
- Hou et al. (2026) -- "A unified framework enables accessible deployment and comprehensive benchmarking of single-cell foundation models"
- Souza & Mehta (2026) -- "Parameter-free representations outperform single-cell foundation models on downstream benchmarks"
- DenAdel et al. (2025) -- "Evaluating the role of pre-training dataset size and diversity on single-cell foundation model performance"

**Evidence:**
- **13 models benchmarked across 50 datasets** in Hou et al.: PCA remained highly competitive and was frequently superior in zero-shot settings
- In cell type classification with full supervision, HVGs+PCA reached ~98% accuracy, "indicating limited additional gains from model complexity under full supervision"
- **Parameter-free scTOP** achieves SOTA or near-SOTA performance, "exceeding the performance of foundation models" on cross-species annotation (F1 > 0.5 vs foundation model F1 < 0.5)
- In perturbation prediction, a simple baseline of predicting "no change" in gene expression **outperformed all fine-tuned scFMs** (DenAdel et al.)
- Complex models (SSL, Geneformer) were exceeded by "both simpler approaches (Pre-trained PCA and scVI) and baselines (PCA and logistic regression)"
- For cancer patient outcome prediction, HVG baseline "slightly outperformed the best foundation models on average" (Elmarakeby et al., 2025)

> "Using simple, interpretable pipelines that rely on careful normalization and linear methods, we obtain SOTA or near SOTA performance across multiple benchmarks commonly used to evaluate single-cell foundation models..." -- Souza & Mehta (2026)

> "The benchmarking challenges the necessity of developing foundation models for single-cell analysis." -- Liu et al. (2024)

---

## 3. Pretraining Scaling Laws Do Not Apply

**Problem:** Unlike NLP/LLMs, larger pretraining datasets and bigger models do not yield proportional improvements in single-cell biology.

**Key Paper:**
- DenAdel et al. (2025) *PMC/BioRxiv* -- "Evaluating the role of pre-training dataset size and diversity on single-cell foundation model performance"

**Evidence:**
- Performance **plateaus** at ~200,000 cells (1% of typical pretraining corpus); saturates before 10%
- For cell type classification, learning saturation occurs at just 1% of the pretraining data
- **Naively increasing dataset diversity does not improve performance**
- The study directly contradicts the "scale is all you need" paradigm:
  > "Further scaling up of pre-training datasets from tens of millions of cells to hundreds of millions or even billions of cells may not yield tangible returns."

- Similar findings from Wang et al. (2025): "no analogous relationship between model size or training data volume and batch-effect removal performance"

- Domain-specific pretraining (e.g., cancer-specific models) significantly outperforms simply scaling up data, suggesting **quality over quantity** (Elmarakeby et al., 2025)

> "Current scRNA-seq cell atlases simply do not contain the kind of deep structure that justifies the use of complex deep learning models." -- Souza & Mehta (2026)

---

## 4. Batch Effects: A Fundamental Barrier

**Problem:** Current scFMs fundamentally fail to create universal cell embeddings because batch effects persist in pretrained representations.

**Key Paper:**
- Wang, Zhang & Zhang (2025) *BioRxiv* -- "Batch Effects Remain a Fundamental Barrier to Universal Embeddings in Single-Cell Foundation Models"

**Evidence:**
- **8 scFMs evaluated** across diverse batch scenarios (tissue, donor, disease, assay, RNA source, species)
- All scFMs exhibit **pronounced batch effects**; batch signals dominate biological signals in embeddings
- scFMs consistently underperformed Harmony and scVI; in some cases fell **below the PCA baseline**
- Probing analysis: embeddings retained **high batch-prediction accuracy**, indicating substantial batch-related information
- Batch centralization only **partially** mitigated batch effects; still failed to outperform traditional methods
- Unlike LLMs (where subtracting language identity recovers semantic structure), scFMs cannot recover biological universality after batch correction
- scCello (the only model with supervised biological priors) performed best, suggesting that learning strategy -- not scale -- is the critical gap

> "Current scFMs fail to effectively remove batch effects, with batch signals persisting in pretrained embeddings."

---

## 5. Tokenization and Representation Failures

**Problem:** The "genes-as-words, cells-as-sentences" analogy is fundamentally flawed, causing information loss.

**Key Papers:**
- Atti & Subramaniam (2025) *BioRxiv* -- "Fundamental Limitations of Foundation Models in Single-Cell Transcriptomics"
- DenAdel et al. (2025)
- Haber et al. (2025) -- "Heimdall: A Modular Framework for Tokenization"

**Evidence:**
- **Statistical model (Seurat v5) outperformed Geneformer by 9-15% accuracy** on cell classification
- bin-based tokenization (used by scGPT, SCMAMBA-2) fails to preserve relative gene expression; Gaussian noise perturbation caused 11% accuracy drop
- Geneformer's rank-value encoding was more robust (8% drop), but still underperformed Seurat
- scGPT overpredicts abundant cell types by 16% (CD14 cells) while underrepresenting rare types
- Compression from raw counts to embedding vectors spans several orders of magnitude, causing "significant loss of information"
- Longer gene sequence lengths **degrade** performance for some models (BioLLM findings on scBERT)

> "Current scFM tokenization schemes that produce 'cellular sentences' may be limited in their ability to represent scRNA-seq data faithfully." -- DenAdel et al. (2025)

> "Foundation models lack critical biological context that would allow for them to make the nuanced inferences required for complex biological analyses." -- Atti & Subramaniam (2025)

---

## 6. Generalization Failures Across Domains

**Problem:** scFMs fail to generalize across species, tissues, assays, and biological contexts.

**Key Papers:**
- Souza & Mehta (2026) -- Cross-species annotation
- Boiarsky et al. (2023) *BioRxiv* -- "A Deep Dive into Single-Cell RNA Sequencing Foundation Models"
- Wenteler et al. (2024) *BioRxiv* -- "PertEval-scFM: Benchmarking Single-Cell Foundation Models for Perturbation Effect Prediction"

**Evidence:**
- Cross-species annotation: foundation models achieve F1-scores consistently **below 0.5** between human and other organisms
- Parameter-free scTOP achieves consistently higher F1-scores than all foundation models across all species
- Foundation models fail to capture evolutionary relationships; simple normalized gene expression vectors show much stronger correlation with evolutionary distance
- In spatial transcriptomics (zero-shot): **none of 13 evaluated models outperformed PCA**
- In trajectory inference: scFM embeddings "may attenuate fine-grained continuous signals"
- Distribution shifts in perturbation data cause failure: "current single-cell foundation model embeddings do not consistently improve prediction accuracy, especially with distribution shifts" (PertEval-scFM)
- Domain sensitivity: a model may work well on PBMC data but perform poorly on pancreas data

---

## 7. Real-World Clinical/Translational Failures

**Problem:** scFMs perform poorly on clinically relevant tasks, particularly patient-level predictions.

**Key Papers:**
- Elmarakeby et al. (2025) -- "Empirical Evaluation of Single-Cell Foundation Models for Predicting Cancer Outcomes"
- Csendes et al. (2025) -- "Benchmarking foundation cell models for post-perturbation RNA-seq prediction"

**Evidence:**
- **9 scFMs evaluated on 6 cancer-specific tasks**: "limited advantage over baselines for predicting patient outcomes"
- HVG baseline slightly **outperformed the best foundation models** on average for patient classification
- Disease stage inference (early vs. late-stage cancer) was highly challenging with "high level of noise"
- End-to-end fine-tuning "did not consistently surpass strong baseline approaches" in patient-undersampled regimes
- **Perturbation prediction**: simpler machine learning approaches and baseline models outperformed foundation models
- Current scFMs treat cells as independent entities; lack tissue context, cell-cell communication, and patient-level representations

> "Current scFMs risk underperforming in clinically important oncology applications and thus demand careful methodological refinement before translational deployment." -- Elmarakeby et al. (2025)

---

## 8. Model Stability and Interpretability Issues

**Problem:** scFMs are unstable, high-variance black boxes with limited interpretability.

**Key Papers:**
- Liu et al. (2024)
- Xu et al. (2026) -- "CellBench-LS"

**Evidence:**
- scFMs exhibit **higher variance across random seeds** compared to scVI and ResPAN
- Cannot be integrated into standard toolboxes (Seurat, Scanpy) due to instability
- Training requires "massive computational infrastructure (e.g., thousands of GPUs)" (Souza & Mehta, 2026)
- Models remain "black boxes with limited interpretability"
- No single model maintains SOTA across all tasks (multi-task consistency failure)
- Aggregation strategy (how cell predictions become patient predictions) "can influence performance as much as, or more than, the choice of embedding model itself" (Elmarakeby et al., 2025)

---

## 9. The "Denoising" Hypothesis

**Problem:** scFMs may function primarily as sophisticated denoisers rather than learning new biology.

**Key Paper:**
- Souza & Mehta (2026) -- "Parameter-free representations outperform single-cell foundation models on downstream benchmarks"

**Evidence:**
> "Current foundation models primarily function as sophisticated denoising procedures rather than tools that uncover fundamentally new biological structure."

The authors argue the underlying biological data manifold is approximately linear, meaning:
- Current scRNA-seq atlases lack the "deep structure" to justify complex deep learning
- Simple linear methods already capture essential biological information
- The cost/benefit ratio of scFMs (massive compute for marginal gains) is unfavorable

---

## Summary Table: Problems vs. Evidence

| Problem | Key Finding | Representative Papers |
|---------|-------------|----------------------|
| Zero-shot failure | HVG/PCA consistently beat scFMs without fine-tuning | Kedzierska 2024, 2025; Boiarsky 2023 |
| Baselines outperform | Simpler methods match/exceed scFMs on most tasks | Liu 2024; Hou 2026; Souza & Mehta 2026; DenAdel 2025 |
| Scaling laws don't apply | Performance plateaus at ~200K cells; diversity doesn't help | DenAdel 2025; Wang 2025 |
| Batch effects persist | All 8 scFMs fail to remove batch signals; worse than PCA | Wang, Zhang & Zhang 2025 |
| Tokenization flawed | Binning loses gene-level information; "cellular sentences" analogy broken | Atti & Subramaniam 2025; DenAdel 2025 |
| No cross-species generalization | F1 < 0.5 across species; linear methods better | Souza & Mehta 2026 |
| Clinical utility limited | HVG beats scFMs on patient outcome prediction | Elmarakeby 2025; Csendes 2025 |
| Instability | High variance across seeds; not production-ready | Liu 2024; Xu 2026 |
| Cost vs. benefit | Massive compute for marginal gains; primarily denoising | Souza & Mehta 2026 |

---

## Paper Index (Selected)

| # | Title | Authors | Year | Venue | ID |
|---|-------|---------|------|-------|----|
| 1 | Fundamental Limitations of Foundation Models in Single-Cell Transcriptomics | Atti & Subramaniam | 2025 | bioRxiv | bio_e2b4691b7ebb |
| 2 | Batch Effects Remain a Fundamental Barrier to Universal Embeddings in Single-Cell Foundation Models | Wang, Zhang & Zhang | 2025 | bioRxiv | bio_d12e0cf38652 |
| 3 | Parameter-free representations outperform single-cell foundation models on downstream benchmarks | Souza & Mehta | 2026 | bioRxiv/PMC | bio_c9329bb26e08 |
| 4 | Evaluating the role of pre-training dataset size and diversity on single-cell foundation model performance | DenAdel et al. | 2025 | PMC | PMC12637645 |
| 5 | Zero-shot evaluation reveals limitations of single-cell foundation models | Kedzierska et al. | 2025 | Genome Biology | PMC12007350 |
| 6 | Assessing the limits of zero-shot foundation models in single-cell biology | Kedzierska et al. | 2024 | bioRxiv | bio_c1d1c1716ff3 |
| 7 | Evaluating the Utilities of Foundation Models in Single-cell Data Analysis | Liu et al. | 2024 | Advanced Science | bio_7eeffc2c6bfb |
| 8 | A unified framework enables accessible deployment and comprehensive benchmarking of single-cell foundation models | Hou et al. | 2026 | bioRxiv/PMC | bio_11a382ac6d5d |
| 9 | Empirical Evaluation of Single-Cell Foundation Models for Predicting Cancer Outcomes | Elmarakeby et al. | 2025 | PMC | PMC12637420 |
| 10 | A Deep Dive into Single-Cell RNA Sequencing Foundation Models | Boiarsky et al. | 2023 | bioRxiv | bio_461ebfe17183 |
| 11 | BioLLM: A Standardized Framework for Integrating and Benchmarking Single-Cell Foundation Models | Qiu et al. | 2025 | Patterns/PMC | PMC12365531 |
| 12 | Benchmarking foundation cell models for post-perturbation RNA-seq prediction | Csendes et al. | 2025 | BMC Genomics | PMC12016270 |
| 13 | PertEval-scFM: Benchmarking Single-Cell Foundation Models for Perturbation Effect Prediction | Wenteler et al. | 2024 | bioRxiv | bio_ef1731f16d53 |
| 14 | CellBench-LS: Benchmark Evaluation of Single-cell Foundation Models for Low-supervision Scenarios | Xu et al. | 2026 | bioRxiv | bio_42703a609557 |
| 15 | Single Cell Foundation Models Evaluation (scFME) for In-Silico Perturbation | Boylan et al. | 2025 | bioRxiv | bio_40d2388176d0 |
| 16 | Benchmarks and definitions of open problems in single-cell analysis | Luecken et al. | 2024 | PMC | PMC11030530 |
| 17 | Heimdall: A Modular Framework for Tokenization in Single-Cell Foundation Models | Haber et al. | 2025 | PMC | PMC12642435 |
| 18 | Single-cell foundation models: bringing artificial intelligence into cell biology | Baek et al. | 2025 | Exp Mol Med | PMC12586647 |

---

## Implications for CheckAtlas

The literature makes clear that **quality control for single-cell foundation models is an urgent and unmet need**. The benchmarking consensus reveals:

1. **No standard QC framework exists for scFMs** -- each paper creates ad-hoc benchmarks
2. **Pretraining data quality** (not quantity) is the bottleneck -- directly maps to CheckAtlas's atlas QC mission
3. **Batch effects** remain the single biggest technical barrier for scFMs
4. **Tokenization quality** (binning, normalization, encoding) is a critical and under-evaluated component
5. **Evaluation metrics are inconsistent** across papers, making cross-model comparison unreliable

This directly motivates extending CheckAtlas from atlas-level QC to **foundation model embedding QC** -- evaluating whether the embeddings produced by scFMs preserve biological signal, remove technical noise, and generalize across conditions.

---

*Report generated via Paperclip systematic literature search and parallel AI analysis of 40+ papers from PMC, bioRxiv, and arXiv.*
