## V-measure ##

## Description ##

La V-measure est une métrique externe fondée sur l’entropie conditionnelle qui évalue la qualité d’un clustering en comparant une partition prédite à une partition de référence.
Elle combine deux propriétés complémentaires :
L'homogénéité :  chaque cluster prédit contient uniquement des objets d’une même classe réelle.
La complétude : tous les objets d’une même classe réelle sont affectés au même cluster prédictif. 
Elle est utile notamment pour vérifier si un clustering de cellules restitue correctement les types cellulaires annotés.

## Formules ## 

Posons : 
- $C$ la vraie partition (classes)
- $K$ la partition prédite (clusters)
- $G(C)$ et $H(K)$ les entropies de chaque partition
- $\mathbb{H} (C \mid K)$ l'entropie conditionnelle de C sachant K
- $\mathbb{H} (K \mid C)$ l'entropie conditionnelle de K sachant C

L'homogénéité est définie par :

$$h = 1 - \frac{H(C \mid K)}{H(C)}$$

La complétude est définie par :

$$c = 1 - \frac{H(K \mid C)}{H(K)}$$

La V-measure est la moyenne harmonique de l'homogénéité et de la complétude :

$$V = \frac{2 \cdot h \cdot c}{h + c}$$

## NB ## 

Soit :

- $n_{i,j}​$ le nombre d'objets dans la classe $i$ de $C$ ET dans le cluster $j$ de $K$
- $n_i$​ le nombre total d'objets dans la classe $i$ de $C$
- $n_j$​ le nombre total d'objets dans le cluster $j$ de $K$
- $n$ le nombre total d'objets

Alors, 

$$H(C \mid K)=-\displaystyle\sum_{j} \frac{n_j}{n} \displaystyle\sum_{i}\frac{n_{i,j}}{n_j}\log{(\frac{n_{i,j}}{n_j})}$$

Plus elle est faible, plus les clusters sont homogènes. 

Et, 

$$H(K \mid C)=-\displaystyle\sum_{i} \frac{n_{i}}{n} \displaystyle\sum_{j}\frac{n_{i,j}}{n_i}\log{(\frac{n_{i,j}}{n_i})}$$

Plus elle est faible, plus le clustering est complet. 

On définit également $V_\beta$ comme étant la V-measure généralisée. Elle permet de pondérer différemment l'homogénéité et la complétude selon les besoins de l'application : 

$$V_\beta=\frac{(1+ \beta)\cdot h \cdot c}{\beta \cdot h + c}$$

où : 
- $h$ est l'homogénéité
- $c$ est la complétude
- $\beta > 0$ est le paramètre de pondération

On retrouve la formule de la V-measure 'standard' en prenant $\beta = 1$. 

## Sources ##

Andrew Rosenberg & Julia Hirschberg (2007). V-Measure: A Conditional Entropy-Based External Cluster Evaluation Measure. Technical Report. Columbia University.

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.v_measure_score.html

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020
