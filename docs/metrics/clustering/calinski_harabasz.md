# Calinski–Harabasz Index

## Description

The Calinski–Harabasz index (CHI), also known as the Variance Ratio Criterion, evaluates clustering quality by comparing the dispersion between clusters to the dispersion within clusters. 
In single-cell analysis, it assesses how well transcriptionally distinct cell populations are separated after clustering. 
A higher CHI indicates dense, well-separated clusters—key for identifying cell types or states. 
It is sensitive to cluster compactness and number, often favoring globular, equally sized groups.

## Formulas : 

Let:
- $n$ be the total number of data points (e.g., cells),
- $k$ be the number of clusters,
- $X$ the full dataset,
- $c_i$ be the centroid of cluster $i$,
- $c$ be the global centroid,
- $n_i$ the number of samples in cluster $i$.

Define:

### **Between-cluster dispersion** :

  $$
  B_k = \sum_{i=1}^k n_i \|c_i - c\|^2
  $$
  
### **Within-cluster dispersion** :

  $$
  W_k = \sum_{i=1}^k \sum_{x \in C_i} \|x - c_i\|^2
  $$

Then the Calinski–Harabasz index is : 

  $$
  CHI = \frac{B_k / (k - 1)}{W_k / (n - k)}
  $$

### **Value range** : 

$[0, +\infty)$ — higher is better. No strict upper bound.

## Sources : 

[scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.calinski_harabasz_score.html)  

[Wikipedia](https://en.wikipedia.org/wiki/Calinski%E2%80%93Harabasz_index)

## Code : 

[scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.calinski_harabasz_score.html)  
