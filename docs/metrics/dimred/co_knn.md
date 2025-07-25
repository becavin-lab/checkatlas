# co-KNN AUC & co-KNN size

![Static Badge](https://img.shields.io/badge/Checkatlas-Disable-red)
![Static Badge](https://img.shields.io/badge/Scikit-Python-yellow)

## Description

The co-KNN AUC and co-KNN size metrics assess how well local neighborhood structures are preserved after dimensionality reduction.  
They are inspired by co-ranking theory and extend traditional nearest-neighbor preservation metrics by comparing local structures between the original and embedded spaces.

- co-KNN size measures the average number of shared neighbors between the high-dimensional and reduced representations. It captures local stability, reflecting to what extent individual neighborhoods are preserved.

- co-KNN AUC evaluates global fidelity of neighborhood preservation by computing the area under the curve (AUC) of a co-ranking–based function across different neighborhood sizes. This function, $QNN(k)$, quantifies how many of the top-$k$ neighbors in the original space are also in the top-$k$ in the reduced space, averaged over all data points.

> Note: While these simplified metrics follow the intuition behind co-ranking matrices (see Zhang et al., 2021), they may differ slightly in formulation from the original definitions (e.g., QNX or QNN(k) curves).

## Formulas

### ![Static Badge](https://img.shields.io/badge/Dimensionality_reduction-local-blue) co-KNN size

This metric estimates how many of the $k$ nearest neighbors are preserved between the high-dimensional space and the reduced space, averaged over all data points:

$$
\text{co-KNN size}(k) = \frac{1}{Nk} \sum_{i=1}^{N} \left| N_k^{\text{orig}}(i) \cap N_k^{\text{embed}}(i) \right|
$$

Where:

- $N$: total number of data points
- To normalize this value into the range $[0, 1]$, we divide by $k$
- $N_k^{\text{orig}}(i)$: the set of $k$ nearest neighbors of point $i$ in the original space  
- $N_k^{\text{embed}}(i)$: the set of $k$ nearest neighbors of point $i$ in the reduced space  

Range: $[0, 1]$ — higher values indicate better local preservation.

### ![Static Badge](https://img.shields.io/badge/Dimensionality_reduction-global-green) co-KNN AUC

This metric is based on the co-ranking matrix $Q$, which compares the ranks of neighbor relationships in the original and embedded spaces.

Let:

- $Q(i,j)$ be the number of points ranked $\leq i$ in the original space and $\leq j$ in the reduced space.  
- $QNN(k) = \frac{1}{N k} \sum_{i=1}^{k} \sum_{j=1}^{k} Q(i,j)$, normalized to fall in $[0,1]$

Then:

$$
\text{co-KNN AUC} = \int_{k=1}^{K} QNN(k) \, dk \approx \frac{1}{K} \sum_{k=1}^{K} QNN(k)
$$

- This corresponds to the area under the $QNN(k)$ curve, which summarizes neighbor preservation quality across different neighborhood sizes.
- A value close to 1 indicates excellent global fidelity of neighborhood preservation.


## Sources

[OpenProblems](https://openproblems.bio/results/dimensionality_reduction?version=v1.0.0)  

[Zhang, Y., Shang, Q., & Zhang, G. (2021). *pyDRMetrics - a python toolkit for dimensionality reduction quality assessment*. Heliyon, 7(2), e06199.](https://doi.org/10.1016/j.heliyon.2021.e06199)  

[Wikipedia - k-nearest neighbors](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm)

## Code

[Scikit - AUC](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score)  

[Scikit - Neighbors](https://scikit-learn.org/stable/modules/neighbors.html)
