# Fowlkes–Mallows index


## Description 

The Fowlkes-Mallows index is an external evaluation method used to determine the similarity between two classifications (clusters). 
This measure can compare either two hierarchical clusterings or a clustering with a reference classification. 
This index is particularly useful for evaluating the performance of clustering algorithms. 
A higher index value indicates greater similarity between the clusters and the reference classifications
The index ranges from 0 (worst possible classification) to 1 (perfect classification). 
It is the geometric mean of precision and recall, making it a robust metric for clustering evaluation. 

## Formulas 

The Fowlkes-Mallows index can be expressed in several ways depending on the context.

### *General Formulation* : 

$$
FM=TPTP+FP⋅TPTP+FNFM=\sqrt{\frac{TP}{TP+FP}\cdot\frac{TP}{TP+FN}}
$$

where :

- $TP$ is the number of True Positives

- $FP$ is the number of False Positives

- $FN$ is the number of False Negatives

  
## Sources 

[Wikipedia](https://en.wikipedia.org/wiki/Fowlkes%E2%80%93Mallows_index)

[E. B. FOWLKES et C. L. MALLOWS. « A method for comparing two hierarchical clusterings ». In : Journal of the American Statistical Association 78 (1983).](https://doi.org/10.2307/2288117)

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

## Code 

[Scikit](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.fowlkes_mallows_score.html#sklearn.metrics.fowlkes_mallows_score)

