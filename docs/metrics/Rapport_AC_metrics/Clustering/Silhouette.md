# Silhouette metrics
## V1 ##
## Summary

## Equation

$$s(i) =\frac{b(i) - a(i))}{\max(a(i),b(i))}$$

where $a(i)$ is the average dissimilarity of point $i$ to all other points in its cluster, and $b(i)$ is the minimum average dissimilarity of point $i$ to all points in any other cluster. The silhouette score for a clustering is the average of the silhouette scores for all data points.

## Description

Silhouette metrics are a family of measures used to evaluate the quality of clustering results. The silhouette score for a single data point measures how similar it is to its own cluster compared to other clusters. The silhouette score ranges from -1 to 1, where a score of 1 indicates that the data point is well-matched to its own cluster, while a score of -1 indicates that it is more similar to neighboring clusters.

A high silhouette score indicates that the clustering is appropriate, as data points are well-matched to their own cluster and dissimilar from other clusters. A low silhouette score suggests that the clustering is suboptimal, as data points are either poorly matched to their own cluster or too similar to neighboring clusters.

Different variations of the silhouette metric exist, such as the weighted silhouette and the silhouette coefficient. These variations adjust for issues such as unbalanced cluster sizes and differing densities between clusters. The silhouette metric is a useful tool for evaluating clustering results and comparing the performance of different clustering algorithms.

## Example
## Source

Peter J. ROUSSEEUW. « Silhouettes : A graphical aid to the interpretation and validation of cluster analysis ». In : Journal of Computational and Applied Mathematics
20 (nov. 1987).

Applying Deep Learning algorithm to perform lung cells annotation, A. Collin

## Implementation

## V2 ##

## Description ##

La Silhouette est une méthode d'interprétation et de validation de la cohérence au sein des clusters de données. Cette technique fournit une représentation graphique succincte de la qualité de classification de chaque objet. La valeur silhouette mesure la similarité d'un objet à son propre cluster (cohésion) comparée aux autres clusters (séparation). Les valeurs varient de -1 à +1, où une valeur élevée indique que l'objet est bien associé à son cluster et mal associé aux clusters voisins. Un clustering avec une largeur silhouette moyenne supérieure à 0.7 est considéré comme "fort", supérieure à 0.5 comme "raisonnable" et supérieure à 0.25 comme "faible". La métrique est spécialisée pour mesurer la qualité des clusters de forme convexe et peut ne pas bien fonctionner avec des clusters de formes irrégulières ou de tailles variables.

## Formule ##

Pour un point de données $i$ dans le cluster $C_I$, on définit :

*Distance intra-cluster (cohésion)*

$a(i)$ représente la distance moyenne entre $i$ et tous les autres points du même cluster.

$$ a(i)=\frac{1}{|C_i|-1}\displaystyle\sum_{j \in C_i}^{j \neq i}d(i,j) $$

où : 
- $C_i$ est le nombre de points dans le cluster
- $d(i,j)$ est la distance entre les points $i$ et $j$

*Distance inter-cluster (séparation)*

$b(i)$ représente la distance moyenne minimale de $i$ à tous les points de n'importe quel autre cluster.

$$b(i)=\displaystyle\min_{k\neq I}\{\frac{1}{|C_k|}\displaystyle\sum_{j \in C_k} d(i,j)\}$$

*Valeur Silhouette* 

La valeur silhouette $s(i)$ est définie par : 

$$s(i) =\frac{b(i) - a(i))}{\max(a(i),b(i))}$$

## Sources ##
https://en.wikipedia.org/wiki/Silhouette_(clustering) 

Peter J. ROUSSEEUW. « Silhouettes : A graphical aid to the interpretation and validation of cluster analysis ». In : Journal of Computational and Applied Mathematics 20 (nov. 1987).

Applying Deep Learning algorithm to perform lung cells annotation, A. Collin


