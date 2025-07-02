# Davies-Bouldin Index 

## Description 

The Davies-Bouldin Index (DBI) is a metric introduced by David L. Davies and Donald W. Bouldin in 1979 to evaluate clustering algorithms.
It is an internal evaluation scheme that measures clustering quality using the characteristics of the dataset itself.
The index simultaneously evaluates two essential criteria : cluster compactness (low internal dispersion of points) and their separation (significant distance between cluster centers).
The lower the index value, the better the cluster separation and the greater the 'compactness' within clusters.

The score is defined as the average of the similarities of each cluster with its most similar cluster, where similarity is the ratio between intra-cluster and inter-cluster distances.

## Formulas

The Davies-Bouldin Index is calculated as follows :

$$DB=\frac{1}{N}\displaystyle\sum_{i=1}^{N} D_i$$

Where : 
- $N$ is the number of clusters,
- $D_i \equiv \displaystyle\max_{j\neq i} R_{i,j}$
- $R_{i,j}=\frac{S_i+S_j}{M_{i,j}}$ where
    -  $S_i$ is the intra-cluster dispersion (average distance of points to the center of cluster i), given by $S_i=(\frac{1}{T_i}\displaystyle\sum_{j=1}^{T_i}\lVert X_j - A_i \rVert_p^q)^{\frac{1}{q}}$
    -  $M_{i,j}$ is the inter-cluster separation (distance between the centers of clusters i and j)

## Sources 

Applying Deep Learning algorithm to perform lung cells annotation, A. COLLIN

[Wikipeda](https://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index)

[Donald W. DAVIES, DAVID L. ; BOULDIN. « A Cluster Separation Measure ». In : IEEE Transactions on Pattern Analysis and Machine Intelligence (1979).](https://doi.org/10.1109/TPAMI.1979.4766909)

## Code 

[Scikit](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.davies_bouldin_score.html#sklearn.metrics.davies_bouldin_score)
