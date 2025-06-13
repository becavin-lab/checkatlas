## Kruskal's stress ##

## Description ##

Métrique basée sur la différence entre les distances en haute et en basse dimension.

Le Stress de Kruskal est une métrique fondamentale pour évaluer la qualité de la représentation en Scaling Multidimensionnel (MDS). Il mesure l'écart entre les distances dans l'espace de représentation réduit et les dissimilarités originales des données. Cette métrique quantifie la perte d'information lors de la réduction dimensionnelle, permettant d'évaluer si la projection en dimensions réduites préserve fidèlement les relations de proximité entre les points. Le Stress varie entre 0 (représentation parfaite) et 1 (représentation très déformée).
Elle permet de valider la qualité des visualisations de données haute dimension (i.e. ensembles de données comportant un grand nombre de variables et posant des probmèmes).

## Formules ##

Stress : 

KS = sqrt( sum_{i<j} (d_ij - d̂_ij)^2 / sum_{i<j} d_ij^2 )

Où : 
- $d_{ij}$ est la distance euclidienne entre les points $i$ et $j$ dans l'espace de représentation réduit
- $\hat{d}_{ij}$ représente la dissimilarité transformée (ou proximité) entre les objets $i$ et $j$ dans l'espace original

** Critères d'interprétation **

- Stress < 0.05 : Excellente représentation
- 0.05 ≤ Stress < 0.1 : Bonne représentation
- 0.10 ≤ Stress < 0.15 : Représentation acceptable
- 0.15 ≤ Stress < 0.20 : Représentation médiocre
- Stress ≥ 0.20 : Représentation inacceptable

Sources : 

https://en.wikipedia.org/wiki/Multidimensional_scaling

https://www.normalesup.org/~carpenti/Notes/MDS/MDS-metrique.html

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

Antonio GRACIA et al. A methodology to compare Dimensionality Reduction algorithms in terms of loss of quality. Rapp. tech. 2014
  
