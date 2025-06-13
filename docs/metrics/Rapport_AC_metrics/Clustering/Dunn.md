## DUNN INDEX ##

## Description ##

Métrique interne, évaluant un compromis entre une distances caractéristique inter et intra cluster

L'indice de Dunn est une métrique pour évaluer les algorithmes de clustering. Il fait partie d'un groupe d'indices de validité incluant l'indice de Davies-Bouldin ou l'indice de Silhouette, en tant que schéma d'évaluation interne basé sur les données clusterisées elles-mêmes. 
L'objectif est d'identifier des ensembles de clusters compacts (faible variance entre les membres du cluster) et bien séparés (moyennes des différents clusters suffisamment éloignées par rapport à la variance intra-cluster). Pour un assignment donné de clusters, un indice de Dunn plus élevé indique un meilleur clustering. L'indice de Dunn quantifie la corrélation entre la compacité et la séparation des clusters en calculant la distance la plus courte entre deux points de clusters différents divisée par la distance la plus longue entre points à l'intérieur d'un cluster.

## Formule ##

L'indice de Dunn est défini comme suit :

$$ DI=\frac{\displaystyle\min_{1\leq i \leq j \leq m} \delta (C_i,C_j)}{\displaystyle\max_{1 \leq k \leq m} \Delta_k} $$

où : 
- $m$ est le nombre de clusters dans l'ensemble
- $\delta (C_i,C_j)$ est la distance inter-cluster entre les clusters $C_i$ et $C_j$
- $\Delta_k$ est la distance intra-cluster (diamètre) du cluster k

*NB* 
Ils existent plusieurs formulations possibles pour la distance inter- et la distance intra-cluster

## Sources ## 

J. C. DUNN. « A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters ». In : Journal of Cybernetics 3 (jan. 1973).

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

https://en.wikipedia.org/wiki/Dunn_index 
