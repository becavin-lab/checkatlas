# High Variable Genes (HVG) overlap 

## Description 

The High Variable Gene (HVG) Overlap metric evaluates how well batch integration preserves biologically informative gene variability in single-cell RNA-seq data.
It compares the sets of highly variable genes (HVGs) identified before and after integration.
A high overlap score indicates that the integration method retains key biological signals across batches.

## Formulas 

Let : 

- $X$ be the set of HVGs before integration

- $Y$ be the set of HVGs after integration

- $\left | X \cap Y \right |$ the number of genes common to both sets

- $\min(X,Y)$ the size of the smaller set

The HVG Overlap is defined as :

$$
overlap(X,Y)=\frac{\left | X \cap Y \right |}{\min(\left |X \right |,\left |Y \right |)}
$$

This formulation emphasizes the preservation of informative genes, even when the integration process alters the number of HVGs. 
The overall HVG score is computed as the mean of per-batch HVG overlap scores.

## Sources 

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022)](https://doi.org/10.1038/s41592-021-01336-8)

[Open Problems in Single-Cell Analysis](https://openproblems.bio/results/batch_integration?version=v2.0.0)

## Code 

[SCIB](https://github.com/theislab/scib/blob/main/scib/metrics/highly_variable_genes.py)


