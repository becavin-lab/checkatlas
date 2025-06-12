**DESCRIPTION** 

Métrique interne basée sur une notion de similarité

L'indice de Davies-Bouldin (DBI) est une métrique introduite par David L. Davies et Donald W. Bouldin en 1979 pour évaluer les algorithmes de clustering. Il s'agit d'un schéma d'évaluation interne qui mesure la qualité du clustering en utilisant les caractéristiques du jeu de données lui-même.
L'indice évalue simultanément deux critères essentiels : la compacité des clusters (faible dispersion interne des points) et leur séparation (distance importante entre les centres des clusters). Plus la valeur de l'indice est faible, meilleure est la séparation des clusters et plus la 'compacité' à l'intérieur des clusters est importante.
Le score est défini comme la moyenne des similarités de chaque cluster avec son cluster le plus similaire, où la similarité est le rapport entre les distances intra-cluster et inter-cluster.

**FORMULE** 

L'indice de Davies-Bouldin se calcule comme suit :

$$DB=\frac{1}{N}\displaystyle\sum_{i=1}^{N} D_i$$

Où : 

- $N$ est le nombre de clusters
- $D_i \equiv \displaystyle\max_{j\neq i} R_{i,j}$
- $R_{i,j}=\frac{S_i+S_j}{M_{i,j}}$ où
  - $S_i$ est la dispersion intra-cluster (distance moyenne des points au centre du cluster i) et est donné par
$S_i=(\frac{1}{T_i}\displaystyle\sum_{j=1}^{T_i}\lVert X_j - A_i \rVert_p^q)^{\frac{1}{q}}$
  - $M_{i,j}$ est la éparation inter-cluster (distance entre les centres des clusters i et j)

**SOURCE**

Applying Deep Learning algorithm to perform lung cells annotation, A. COLLIN
https://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index 
