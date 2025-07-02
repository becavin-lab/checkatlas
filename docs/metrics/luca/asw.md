# Average Silhouette Width (ASW) Metrics 

## Description 

The Average Silhouette Width (ASW) metrics evaluate how well cells are grouped or mixed after batch integration in single-cell RNA-seq data. There are three main variants:

- Batch ASW: Measures how well cells from different batches are mixed within the same biological group (e.g., cell type).

- Label ASW: Measures how well cells of the same biological identity (e.g., cell type) cluster together.

- Isolated Label ASW: Focuses on rare or isolated cell types to assess whether they remain well-separated after integration.

Each metric uses the silhouette score, which quantifies how similar a cell is to its own cluster compared to other clusters.

## Formulas 

> Silhouette Clustering Metric Analysis
> 
> We remind that the general silhouette width is define, for a cell $i$, as follows :
> 
> $$S(i) = \frac{(b_i - a_i)}{max(a_i, b_i)}$$
>
> where,
>
> - $a_i$ is the average distance between cell $i$ and all other cells in the same cluster
>
> - $b_i$ is the minimum average distance between cell $i$ and all cells in other clusters
>
> [For further infomations](silhouette.md)

### *Batch ASW (non-scaled version)* : 

For all cells $i$ of cell type $C_j$, 

$$ batchASW_j=\frac{1}{\left| C_j \right|} \cdot\displaystyle\sum_{i \in C_j} \left| silhouette(i) \right|$$

Final score across all cell types $M$ :

$$batchASW=\frac{1}{\left| M \right|} \cdot\displaystyle\sum_{j \in M} \left| batchASW_j \right|$$


### *Batch ASW (scaled version - default)* : 

$$batchASW_j=\frac{1}{\left| C_j \right|}\cdot\displaystyle\sum_{i \in C_j} \left| 1-silhouette(i) \right|$$

*NB : The non-scaled version uses the absolute silhouette values directly, where lower scores indicate better batch mixing. 
The scaled version inverts this interpretation, making higher scores indicate better integration, which is more intuitive and consistent with other evaluation metrics.*

### *Cell identity batch* : 

Same formula as classical ASW but calculated using cell type labels as reference clusters instead of batches.

### *Isolated Label ASW*

Same formula as the general silhouette score, but clusters are defined by cell type labels instead of batches. High values indicate that cells of the same type are well-clustered.

## Sources 

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-8)

[Open Problems](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[SCIB Documentation](https://scib.readthedocs.io/en/latest/api/scib.metrics.silhouette_batch.html)

## Code 

[Github](https://github.com/theislab/scib/blob/main/scib/metrics/silhouette.py)
