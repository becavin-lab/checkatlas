### V-measure


## Description 

The V-measure is an external metric based on conditional entropy that evaluates clustering quality by comparing a predicted partition to a reference partition.
It combines two complementary properties : 
- Homogeneity: each predicted cluster contains only objects from the same true class.
- Completeness: all objects from the same true class are assigned to the same predicted cluster.

It is particularly useful for verifying whether cell clustering correctly reconstructs annotated cell types.

## Formulas 

Let:
- $C$ be the true partition (classes)
- $K$ be the predicted partition (clusters)
- $H(C)$ and $H(K)$ be the entropies of each partition
- $H(C \mid K)$ be the conditional entropy of $C$ given $K$
- $H(K \mid C)$ be the conditional entropy of $K$ given $C$

Homogeneity is defined by :

$$h = 1 - \frac{H(C \mid K)}{H(C)}$$

Completeness is defined by : 

$$c = 1 - \frac{H(K \mid C)}{H(K)}$$

The V-measure is the harmonic mean of homogeneity and completeness :

$$V = \frac{2 \cdot h \cdot c}{h + c}$$

*NB* : 

Let:
- $n_{i,j}$ be the number of objects in class $i$ of $C$ AND in cluster $j$ of $K$
- $n_i$ be the total number of objects in class $i$ of $C$
- $n_j$ be the total number of objects in cluster $j$ of $K$
- $n$ be the total number of objects

Then, 

$$H(C \mid K)=-\displaystyle\sum{j} \frac{nj}{n} \displaystyle\sum{i}\frac{n_{i,j}}{nj}\log{(\frac{n{i,j}}{nj})}$$

The lower it is, the more complete the clustering is.

And,

$$H(K \mid C)=-\displaystyle\sum{i} \frac{n{i}}{n} \displaystyle\sum{j}\frac{n_{i,j}}{ni}\log{(\frac{n{i,j}}{ni})}$$


The lower it is, the more complete the clustering is.

We also define $V_\beta$ as the generalized V-measure. It allows different weighting of homogeneity and completeness according to application needs : 

$$V_\beta=\frac{(1+ \beta)\cdot h \cdot c}{\beta \cdot h + c}$$

where:
- $h$ is homogeneity
- $c$ is completeness
- $\beta > 0$ is the weighting parameter

We recover the 'standard' V-measure formula by taking $\beta = 1$.

## Sources 

Andrew Rosenberg & Julia Hirschberg (2007). V-Measure: A Conditional Entropy-Based External Cluster Evaluation Measure. Technical Report. Columbia University.

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.v_measure_score.html

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

## Code 

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_completeness_v_measure.html
