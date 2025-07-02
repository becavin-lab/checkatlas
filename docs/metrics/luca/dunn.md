# Dunn Index Metric Analysis

## Description 

The Dunn index is a metric for evaluating clustering algorithms. 
The objective is to identify sets of compact clusters (low variance between cluster members) and well-separated clusters (means of different clusters sufficiently far apart compared to within-cluster variance). 
For a given cluster assignment, a higher Dunn index indicates better clustering.
The Dunn index quantifies the correlation between cluster compactness and separation by calculating the shortest distance between two points from different clusters divided by the longest distance between points within a cluster.

*It is part of a group of validity indices including the Davies-Bouldin index or the Silhouette index, as an internal evaluation scheme based on the clustered data itself.*

## Formulas 

The Dunn index is defined as follows :

$$
DI=\frac{\displaystyle\min_{1\leq i \leq j \leq m} \delta (C_i,C_j)}{\displaystyle\max_{1 \leq k \leq m} \Delta_k}
$$

Where : 

- $m$ is the number of clusters in the set

- $\delta (C_i,C_j)$ is the inter-cluster distance between clusters $C_i$ and $C_j$

- $\Delta_k$ is the intra-cluster distance (diameter) of cluster k

## Sources

[J. C. DUNN. « A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters ». In : Journal of Cybernetics 3 (jan. 1973).](http://dx.doi.org/10.1080/01969727308546046)

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

[Wikipedia](https://en.wikipedia.org/wiki/Dunn_index)

## Code 

[Matlab](https://www.mathworks.com/matlabcentral/fileexchange/27859-dunn-s-index)

