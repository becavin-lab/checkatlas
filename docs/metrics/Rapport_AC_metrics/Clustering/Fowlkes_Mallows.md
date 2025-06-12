## Fowlkes–Mallows index 

## Description ##

Métrique externe comparant les taux d’accord et de désaccord de deux clustering

L'indice de Fowlkes-Mallows est une méthode d'évaluation externe utilisée pour déterminer la similarité entre deux classifications (clusters). Cette mesure peut comparer soit deux classifications hiérarchiques, soit une classification avec une classification de référence. Cet indice est particulièrement utile pour évaluer la performance des algorithmes de clustering. Une valeur plus élevée de l'indice indique une plus grande similarité entre les clusters et les classifications de référence. L'indice varie entre 0 (pire classification possible) et 1 (classification parfaite). Il s'agit de la moyenne géométrique de la précision et du rappel, ce qui en fait une métrique robuste pour l'évaluation de clustering.

## Formule ##

L'indice de Fowlkes-Mallows peut être exprimé de plusieurs façons selon le contexte.

**Formulation Générale**

$$FM=\sqrt{\frac{TP}{TP+FP}\cdot\frac{TP}{TP+FN}}$$

où : 
- $TP$ est le nombre de vrais positifs (True Positives)
- $FP$ est le nombre de faux positifs (False Positives)
- $FN$ est le nombre de faux négatifs (False Negatives)

## Source ##

https://en.wikipedia.org/wiki/Fowlkes%E2%80%93Mallows_index

E. B. FOWLKES et C. L. MALLOWS. « A method for comparing two hierarchical clusterings ». In : Journal of the American Statistical Association 78 (1983).

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020
