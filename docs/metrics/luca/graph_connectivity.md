# Graph Connectivity

## Description 

The "Graph connectivity" metric evaluates data integration quality by measuring whether biologically similar cells remain connected in a k-nearest neighbors (kNN) graph after integration. 
This metric is particularly important for preserving local data structure during batch effect correction. 
It quantifies how well cells of the same biological type maintain their connectivity in the integrated space, regardless of their batch of origin. 
Good integration should maintain connections between biologically similar cells while removing batch-related artifacts. 
The metric typically uses a kNN graph constructed on the integrated space and compares observed connectivity with expected connectivity based on known biological labels.

## Formulas 

Graph connectivity is typically calculated as follows :

$$GC = \frac{1}{N} \cdot \sum_{i=1}^{N} \frac{|N_i \cap S_i|}{k}$$

Where:
- $N$ = total number of cells
- $N_i$ = set of k nearest neighbors of cell i in the integrated space
- $S_i$ = set of cells with the same biological type as cell i
- $k$ = number of neighbors considered
- $|N_i \cap S_i|$ = number of neighbors of i that belong to the same biological type

*Value Range and Interpretation*
- Range : 0 to 1
- Values close to 1 indicate excellent biological connectivity, where cells of the same biological type remain well-connected in the kNN graph after integration, preserving local biological structure. 
- Values near 0 suggest poor biological connectivity, with biologically similar cells being dispersed or disconnected in the integrated space. 
- Intermediate values represent moderate connectivity with variable performance across different cell types. This metric is considered a bio-conservation measure that evaluates whether natural biological relationships are maintained after integration.

  
## Sources 

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022)](https://doi.org/10.1038/s41592-021-01336-8)

[OpenProblems](https://openproblems.bio/results/batch_integration)

[scIB](https://github.com/theislab/scib/tree/main)

## Code 

[Github](https://github.com/theislab/scib/blob/main/scib/metrics/graph_connectivity.py)
