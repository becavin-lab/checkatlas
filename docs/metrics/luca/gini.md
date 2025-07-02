# Gini coefficient

## Description 

The Gini coefficient is a statistical metric used to measure inequality within a distribution. 
Originally used in economics to measure income inequality, , it has been adapted in bioinformatics to assess how uniformly a gene is expressed across multiple samples.
A low Gini value indicates stable expression across all samples, while a high Gini indicates variability, making the gene less suitable for normalization.

## Formulas 

The Gini coefficient for a gene $g$, across $K$ samples, is defined as :

$$
G=\frac{K+1}{K} - \frac{2}{K \displaystyle\sum_{i=1}^{K} x_{g,i}} \cdot \displaystyle\sum_{i=1}^{K} (K+1-i)\cdot x_{g,i}
$$

Where : 
- $K$ is the number of samples
- $x_{g,i}$ : expression value of gene $g$ in sample $i$
- Values $x_{g,i}$ must be sorted in ascending order ($x_{g,1} \leq x_{g,2} \leq \dots \leq x_{g,K}$)

The formula captures cumulative expression weighted by rank. 

The first term $\frac{K+1}{K}$ represents the theoretical maximum (perfect equality), while the second term computes deviation from that ideal, normalized by the total expression $\sum x_{g,i}$.

## Sources 

[Wikipedia](https://en.wikipedia.org/wiki/Gini_coefficient)

[Marina WRIGHT MUELAS et al. « The role and robustness of the Gini coefficient as an unbiased tool for the selection of Gini genes for normalising expression profiling data ». In : Scientific Reports 9 (déc. 2019).](https://doi.org/10.1038/s41598-019-54288-7)

## Code
