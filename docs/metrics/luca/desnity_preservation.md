# Desnity Preservation 

## Description 


The Density Preservation metric evaluates how well a dimensionality reduction method maintains the local density of points between the original high-dimensional space and the reduced low-dimensional space. 
This is crucial for avoiding misleading visual artifacts where dense regions in the original data appear artificially sparse or overly compact after projection. 
It is especially relevant in single-cell RNA-seq data, where local cell density often reflects biologically meaningful states.

## Formulas 

*Local Density Estimation*

For each point $i$, the local density $\rho_i$ is estimated as the inverse of the average distance to its $k$-nearest neighbors:

$$\rho_i = \frac{1}{\frac{1}{k} \sum_{j \in \text{KNN}(i)} d(i, j)}$$

where $d(i,j)$ is the distance between points $i$ and $j$, and $\text{KNN}(i)$ is the set of the $k$-nearest neighbors of $i$.

*Spearman Correlation*

The final metric is the Spearman correlation between the local density vectors in the high-dimensional and low-dimensional spaces:

$$\text{Density Preservation} = \text{SpearmanCorr}(\rho^{\text{high-dim}}, \rho^{\text{low-dim}})$$

*Value Range*

The metric ranges from -1 (complete inversion of density) to 1 (perfect preservation). A value near 0 indicates no correlation between densities.

## Sources 

[Narayan A, Berger B, Cho H. *Assessing single-cell transcriptomic variability through density-preserving data visualization*. Nat Biotechnol. 2021](https://pubmed.ncbi.nlm.nih.gov/33462509/)

[OpenProblems](https://openproblems.bio/results/dimensionality_reduction?version=v1.0.0)

## Code 

[Density?](https://scikit-learn.org/stable/modules/density.html)
