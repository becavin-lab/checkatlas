# Kendall's tau coefficient

## Description 

Kendall’s tau is a non‑parametric rank correlation coefficient introduced in 1938. 
It quantifies the ordinal association between two variables based on the number of concordant vs discordant pairs. 
The coefficient ranges from –1 (perfect inversion) to +1 (perfect agreement), with 0 indicating no association.
The version defined by Yanai et al. (2005) is not a rank‑correlation measure but a tissue‑specificity index for genes, bounded between 0 (housekeeping) and 1 (tissue‑specific)

## Formulas 

Let $x_{g,c}$ denote the expression of gene $g$ in tissue $c$, xith $k$ total tissues. 
We define the normalized expression : 

$$
\hat{x}_{g,c}=\frac{x_{g,c}}{\displaystyle\max_{1\leq i \leq k}x_{g,i}}
$$

For each gene $g$, compute : 

$$\tau_g=\frac{\displaystyle\sum_{i=1}^{k} (1 - \hat{x}_{g,i})}{k-1}$$

THis $\tau_g$ varies from 0 (ubiquitously expressed) to 1 (specific to a single tissue).

## Sources 

[Wikipedia](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient)

[Wikipedia](https://en.wikipedia.org/wiki/Kendall_tau_distance)

[Itai YANAI et al. « Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification ». In : Bioinformatics 21 (mar. 2005)](10.1093/bioinformatics/bti042)

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

## Code 


