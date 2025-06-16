## Silhouette Clustering Metric Analysis

### Description

The Silhouette is a method for interpretation and validation of consistency within data clusters. This technique provides a succinct graphical representation of how well each object has been classified.
The silhouette value measures how similar an object is to its own cluster (cohesion) compared to other clusters (separation). Values range from $-1$ to $+1$, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters.
The metric is specialized for measuring cluster quality when clusters are convex-shaped and may not perform well if data clusters have irregular shapes or varying sizes.

A clustering with an average silhouette width greater than $0.7$ is considered "strong", greater than $0.5$ "reasonable" and greater than $0.25$ "weak". 

### Formulas

For a data point $i$ in cluster $C_I$, we define:

*Intra-cluster distance (cohesion)*

$a(i)$ represents the average distance between $i$ and all other points in the same cluster.

$$ a(i)=\frac{1}{|Ci|-1}\displaystyle\sum{j \in C_i}^{j \neq i}d(i,j) $$

where : 
- $C_i$ is the number of points in the cluster
- $d(i,j)$ is the distance between points $i$ and $j$

*Inter-cluster distance (separation)*

$b(i)$ represents the minimum average distance from $i$ to all points in any other cluster.

$$b(i)=\displaystyle\min_{k\neq I}\{\frac{1}{|C_k|}\displaystyle\sum_{j \in C_k} d(i,j)\}$$

*Silhouette Value*

The silhouette value $s(i)$ is defined by :

$$s(i) =\frac{b(i) - a(i))}{\max(a(i),b(i))}$$


### Sources 

https://en.wikipedia.org/wiki/Silhouette_(clustering)

Peter J. ROUSSEEUW. « Silhouettes : A graphical aid to the interpretation and validation of cluster analysis ». In : Journal of Computational and Applied Mathematics 20 (nov. 1987).

Applying Deep Learning algorithm to perform lung cells annotation, A. Collin

### Code 

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html 

