# Isolated Label F1 Score

## Description

The Isolated Label F1 score is a metric used to evaluate how well a data integration method preserves rare or isolated cell types after batch correction in single-cell RNA-seq data. 
These isolated labels typically represent biologically meaningful but underrepresented populations that are often lost or merged during integration. 
This metric focuses on classification performance (precision and recall) specifically for these rare cell types, ensuring that integration does not compromise their identity.

## Formulas

The Isolated Label F1 score is computed as the average F1-score over a subset of isolated cell type labels. 
These labels are defined as those with few neighbors in the integrated space or low frequency across batches.

For each isolated label $l$:

- Precision : $P_l=\frac{TP_l}{TP_l + FP_l}$

- Recall : $R_l=\frac{TP_l}{TP_l + FN_l}$

- F1-score : $F1_l=\frac{2 \cdot P_l \cdot R_l}{P_l + R_l}$

Then, the Isolated Label F1 score is the mean over all isolated labels :

$$
\text{Isolated Label F1 score} = \frac{1}{\left | L_{iso} \right |} \sum_{l \in L_{iso}} F1_l
$$

Where $L_{iso}$ is the set of isolated labels.

## Sources 

[OpenProblems.](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[Luecken, M.D. et al. *Benchmarking atlas-level data integration in single-cell genomics*. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-)

[Christensen E. et al. *Evaluation of single-cell RNAseq labelling algorithms using cancer datasets*. Brief Bioinform. 2023.](https://doi.org/10.1093/bib/bbac561)
