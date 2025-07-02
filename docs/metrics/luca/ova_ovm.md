# One‑vs‑All (OvA) & One‑vs‑Max (OvM) 

## Description 

One vs All and One vs Max specificity metrics are fundamental statistical tools in single-cell analysis for quantifying cell type-specific differential expression. 
The OvA metric compares gene expression in a given cell type to all other cell types combined, while OvM compares expression to that of the cell type where the gene is most highly expressed (excluding the studied type). 
hese metrics enable identification of specific marker genes and characterization of the unique transcriptomic signature of each cellular population.
By combining these two approaches, we obtain a robust measure of expression specificity that limits false positives and improves cell annotation accuracy.
These tools are particularly valuable for integrating GWAS data (genome-wide association study) with single-cell transcriptomic profiles, allowing association of genetic variants with specific cell types.
To sum up, 

- OvA metric favors genes with highly localized expression

- OvM metric identifies differentially expressed genes even if they are not exclusive to a single cell type

## Formulas 

### *One-vs-All (OvA)*

Mean expression calculation :

$$
x_{g,c}=\frac{\displaystyle\sum_{i\in c} r_{g,i}}{N_c}
$$

where, 

- $x_{g,c}$ is the average expression of gene $g$ in cell type $c$

- $r_{g,i}$ represents the expression of gene $g$ in cell $i$

- $N_c$ is the number of cells of type $c$

Specificity score :

$$
s_{g,c}=\frac{x_{g,c}}{\displaystyle\sum_{j=1}^{K} x_{g,j}}
$$

where, 

- $s_{g,c}$ represents the specificity score normalizing average expression by the sum of average expressions across all $K$ cell types

- Interval: $[0,1]$, where 1 indicates exclusive expression in cell type $c$

### *One-vs-Max (OvM)*

Ratio calculation :

$$
p_{g,c}=\frac{x_{g,c}}{x_{g,d}}
$$

where, 

- $p_{g,c}$ is the expression ratio between cell type $c$ and cell type $d$ where gene $g$ is most highly expressed (excluding $c$)

- $d$ is the cell type with highest expression of gene $g$ (excluding $c$)

- Interval: $[0,+\infty[$, where values $>1$ indicate that $c$ expresses the gene more strongly than the second highest expressing cell type


## Sources 

[Wikipedia](https://en.wikipedia.org/wiki/Multiclass_classification#One-vs.-rest)

[Skene NG et al. (2018) Genetic identification of brain cell types underlying schizophrenia. Nature Genetics, 50(6):825–833.](https://doi.org/10.1038/s41588-018-0129-5)

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

[Wikipedia](https://en.wikipedia.org/wiki/Genome-wide_association_study)

## Code 


