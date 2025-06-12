## RAND INDEX ##

## Description ##

Métrique externe basée sur la consistance entre deux clustering

L'indice de Rand est une métrique statistique utilisée pour mesurer la similarité entre deux partitions d'un même ensemble de données, notamment pour évaluer la qualité des algorithmes de clustering.Il mesure la fréquence d'occurrence des accords sur l'ensemble des paires, ou la probabilité que deux partitions s'accordent sur une paire choisie aléatoirement.
Cette métrique peut être interprétée comme le pourcentage de décisions correctes prises par un algorithme de clustering. L'indice varie entre 0 et 1, où 0 indique que les deux clusterings ne s'accordent sur aucune paire de points, et 1 indique que les clusterings sont exactement identiques.
La métrique fonctionne en comptant quatre types de paires : les accords positifs ($a$), les accords négatifs ($b$), les désaccords de type 1 ($c$) et les désaccords de type 2 ($d$). L'indice de Rand représente intuitivement le rapport entre les accords totaux et le nombre total de paires possibles dans l'ensemble de données.

## Formule ##

$$R=\frac{(a+b)}{(a+b+c+d)}=\frac{(a+b)}{\binom{2}{n}}$$

Où : 
- $a$ est le nombre de paires d'éléments qui sont dans le même sous-ensemble dans $X$ et dans le même sous-ensemble dans $Y$
- $b$ est le nombre de paires d'éléments qui sont dans des sous-ensembles différents dans $X$ et dans des sous-ensembles différents dans $Y$
- $c$ est le nombre de paires d'éléments qui sont dans le même sous-ensemble dans $X$ et dans des sous-ensembles différents dans $Y$
- $d$ est le nombre de paires d'éléments qui sont dans des sous-ensembles différents dans $X$ et dans le même sous-ensemble dans $Y$
- $\binom{2}{n}$ est le nombre total de paires possibles, définit par : $\binom{2}{n}=\frac{n(n-1)}{2}$

## Sources ##
https://en.wikipedia.org/wiki/Rand_index

Lawrence HUBERT et Phipps ARABIE. « Comparing partitions ». In : Journal of
Classification 2 (déc. 1985).

Applying Deep Learning algorithm to perform lung cells annotation, A. Collin
