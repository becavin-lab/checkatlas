## Density-Based Clustering Validation ##

## Description ##

Métrique interne pour l’évalutation de méthodes de clustering à densité.

*Version plus "sérieuse"*

DBCV (Density-Based Clustering Validation) est une métrique conçue pour évaluer la qualité des solutions de clustering (notamment pour les algorithmes de clustering basés sur la densité). 
Cette métrique évalue la qualité du clustering basée sur la connexion de densité relative entre les paires d'objets.
DBCV est particulièrement adaptée pour identifier les clusters concaves et imbriqués. 
L'indice DBCV évalue la qualité du clustering en évaluant la séparation basée sur la densité entre les clusters et la cohésion au sein des clusters. 
Ses valeurs sont comprises entre -1 et +1, DBCV étant un indice de maximisation, des valeurs plus élevées correspondent à de meilleures partitions.

*Version plus explicative, compréhensible*

DBCV est comme un "juge" qui note la qualité de vos groupes (clusters) de données. Si on doit trier des données, DBCV nous dit si le tri est bon ou mauvais. 
Cette métrique est spécialement conçue pour les algorithmes comme DBSCAN qui trouvent des groupes de formes bizarres (pas forcément ronds). 
DBCV regarde deux choses : est-ce que les points dans un même groupe sont bien collés ensemble ? Et est-ce que les différents groupes sont bien séparés ? 
La note va de -1 (très mauvais) à +1 (très bon).

## Formules ##

*Étape 1*:

Pour chaque groupe, calculer sa "densité interne"

$\Rightarrow$ On regarde à quel point les points sont serrés dans le groupe

*Étape 2*:

Calculer la "séparation entre groupes"

$\Rightarrow$ On mesure l'espace vide entre les différents groupes

*Étape 3*:

Formule finale : DBCV = Moyenne de tous les scores des groupes

$\Rightarrow$ Chaque score de groupe = (Séparation du groupe) - (Densité interne du groupe)

$\Rightarrow$ Plus la séparation est grande ET plus la densité interne est forte $\rightarrow$ meilleur score

## Sources ##

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

https://en.wikipedia.org/wiki/Density-based_clustering_validation

Davoud MOULAVI et al. Density-Based Clustering Validation. Rapp. tech.
