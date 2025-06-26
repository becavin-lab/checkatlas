# co-KNN AUC & co-KNN size

## Description 

The co-KNN AUC and co-KNN size metrics assess how well local neighborhood structures are preserved after dimensionality reduction. 
They are based on the concept of common k-nearest neighbors (co-KNN), which are the neighbors shared between the original high-dimensional space and the reduced space.
- co-KNN AUC evaluates the global fidelity of local structure preservation using a ROC curve.
- co-KNN size measures the average number of shared neighbors, reflecting local stability.
These metrics are particularly useful for benchmarking dimensionality reduction methods in biological data analysis, such as single-cell RNA-seq.

## Formulas 

*co-KNN size* : 

$$\text{co-KNN size} = \frac{1}{N} \displaystyle\sum_{i=1}^{N} \left| N_k^{\text{orig}}(i) \cap N_k^{\text{embed}}(i) \right|$$

Where : 
- $N$: total number of data points.
- $N_k^{\text{orig}}(i)$: the set of $k$ nearest neighbors of point $i$ in the original space.
- $N_k^{\text{embed}}(i)$: the set of $k$ nearest neighbors of point $i$ in the reduced space.
- Range: $[0, k]$, where a higher values indicate better local structure preservation.

*co-KNN AUC* :

Construct a [ROC](https://en.wikipedia.org/wiki/Receiver_operating_characteristic) curve by treating each pair of points as a positive example if they are co-KNN, and negative otherwise.

Compute the Area Under the Curve (AUC):

$$
\text{co-KNN AUC} = \text{AUC}(\text{ROC}_{\text{co-KNN}})
$$

Range: $[0, 1]$ where a value of 1 indicates perfect neighborhood preservation.

## Sources 

[OpenProblems](https://openproblems.bio/results/dimensionality_reduction?version=v1.0.0)

[Zhang, Y., Shang, Q., & Zhang, G. (2021). *pyDRMetrics - a python toolkit for dimensionality reduction quality assessment*. Heliyon, 7(2), e06199.](https://doi.org/10.1016/j.heliyon.2021.e06199)

[Wikipedia](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm)

## Code 

[Scikit - AUC](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score)

[Scikit](https://scikit-learn.org/stable/modules/neighbors.html)
