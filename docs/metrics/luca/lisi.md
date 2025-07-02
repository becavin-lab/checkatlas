# Local Inverse Simpson's Index (iLISI & cLISI)

## Description 

The Local Inverse Simpson's Index (LISI) measures local mixing by estimating the effective number of classes in local neighborhoods of cells and constitutes a fundamental metric for evaluating single-cell data integration quality.
Two main variants exist :

- iLISI (integration LISI) which denotes the effective number of datasets in a neighborhood to assess batch mixing. Higher iLISI values indicate better batch mixing.

- cLISI (conservation LISI) which measures the preservation of cell types after integration. Higher cLISI values suggest better preservation of biological structure.

The LISI of a cell is defined as the effective number of batches, properly scaled, among its k nearest neighbors.
This metric allows objective quantification of integration algorithm performance by simultaneously measuring batch effect correction and biological information conservation.

## Formulas 

LISI is based on the inverse Simpson index applied locally. For a given cell :

### *Local Simpson Index Calculation* :

$$
Simspon_i=\displaystyle\sum_{j} p_{j}^{2}
$$

where $p_j$ is the proportion of cells from category $j$ (batch for iLISI, cell type for cLISI) in the k-nearest neighbor neighborhood of cell $i$.

*LISI (Local Inverse Simpson Index)* :

$$LISI_i=\frac{1}{Simpson_i}=\frac{1}{\displaystyle\sum_{j} p_{j}^{2}}$$

*Interpretation* : 

- iLISI: effective number of batches in the local neighborhood

  - Minimum value = 1 (homogeneous neighborhood, single batch)

  - Maximum value = total number of batches (perfect mixing)

- cLISI: effective number of cell types in the local neighborhood

  - Low value = good conservation (homogeneous neighborhood for one cell type)

  - High value = potential over-mixing (loss of biological structure)

*Normalized Version* :

By default, this function returns a value scaled between 0 and 1 instead of the original LISI range from 0 to the number of batches.

$$LISI_{normalized}=\frac{LISI-1}{N_{categories} - 1}$$

where $N_{categories}$ is the total number of categories (batches or cell types).

## Sources 

[Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019).](https://doi.org/10.1038/s41592-019-0619-0) &
[Github](https://github.com/immunogenomics/LISI)

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-8)

[SCIB](https://scib.readthedocs.io/en/latest/api/scib.metrics.ilisi_graph.html)

[Wikipedia](https://en.wikipedia.org/wiki/Diversity_index)

[Open Problems](https://openproblems.bio/results/batch_integration?version=v2.0.0)

## Code 


