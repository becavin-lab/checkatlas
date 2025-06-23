# Average Silhouette Width of Batch & Cell Identity Labels

## Description 

The Average Silhouette Width (ASW) measures the quality of cluster separation and homogeneity by calculating the average silhouette width. In the context of single-cell batch integration, two variants are used :
- Batch ASW: Evaluates whether cells from different batches mix properly after batch effect correction
- Cell Identity ASW: Evaluates whether cells of the same cell type remain grouped together after integration

For Batch ASW, a silhouette width close to 0 represents perfect overlap of batches, so the absolute value is used to measure mixing quality.
The silhouette width ranges between -1 and 1, with positive values indicating good cluster separation. 
The scaled version (default) provides scores between 0 and 1, where 1 indicates optimal label representation and 0 indicates suboptimal representation.


## Formulas 

> Silhouette Clustering Metric Analysis
> 
> We remind that the general silhouette width is define as follows :
> 
> $$S(i) = \frac{(b_i - a_i)}{max(a_i, b_i)}$$
>
> where,
> - $a_i$ is the average distance between cell $i$ and all other cells in the same cluster
> - $b_i$ is the minimum average distance between cell $i$ and all cells in other clusters
>
> [For further infomations](silhouette.md)

*Batch ASW (non-scaled version)* : 

For all cells $i$ of cell type $C_j$, 

$$ batchASW_j=\frac{1}{\left| C_j \right|} \cdot\displaystyle\sum_{i \in C_j} \left| silhouette(i) \right|$$

Final score averaged over all cell types M : 

$$batchASW=\frac{1}{\left| M \right|} \cdot\displaystyle\sum_{j \in M} \left| batchASW_j \right|$$


*Batch ASW (scaled version - default)* : 

$$batchASW_j=\frac{1}{\left| C_j \right|}\cdot\displaystyle\sum_{i \in C_j} \left| 1-silhouette(i) \right|$$

*NB : The non-scaled version uses the absolute silhouette values directly, where lower scores indicate better batch mixing. 
The scaled version inverts this interpretation, making higher scores indicate better integration, which is more intuitive and consistent with other evaluation metrics.*

*Cell identity batch* : 

Same formula as classical ASW but calculated using cell type labels as reference clusters instead of batches.

## Sources 

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-8)

[Open Problems](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[SCIB Documentation](https://scib.readthedocs.io/en/latest/api/scib.metrics.silhouette_batch.html)

## Code 

[Github](https://github.com/theislab/scib/blob/main/scib/metrics/silhouette.py)
