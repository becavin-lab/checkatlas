# Mutual Information

## Description 

Mutual Information (MI) is an external measure based on information theory that quantifies the dependence between two partitions or random variables.
It measures the amount of mutually shared information : in practice, how much knowledge of partition $U$ reduces uncertainty about partition $V$, and vice versa.
MI evaluates the similarity between a reference label set (for example, annotated cell types) and a partition computed by a clustering algorithm. This measure is symmetric, non-negative, and relies on the marginal and joint entropies of the partitions.
A high MI means that the two partitions share a lot of information (the clusters align well with the true cell types), while a null MI indicates that they are completely independent.

## Formulas 

For two partitions $U$ and $V$ of $n$ objects (with $R$ clusters in $U$ and $C$ in $V$), we define the marginal probabilities

$$
p(i)=\frac{|U_i|}{n} \text{ and } p(j)\frac{|V_j|}{n}
$$

,as well as the joint probability 

$$
p(i,j)=\frac{|U_i \cup V_j|}{n}
$$.

The entropy of $U$ is written as $H(U)=-\displaystyle\sum_{i=1}^{R} p(i) \log{p(i)}$ and, similarly, for V, $H(V)=-\displaystyle\sum_{j=1}^{C} p(j) \log{p(j)}$. 

Thus, the mutual information (MI) is expressed as :

$$
MI(U,V)=\displaystyle\sum_{i=1}^{R}\sum_{j=1}^{C} p(i,j) \log\frac{p(i,j)}{p(i)\cdot p(j)}
$$

## Sources 

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

[Nguyen Xuan VINH, Julien EPPS et James BAILEY. « Information theoretic measures for clusterings comparison : Is a correction for chance necessary ? » In : ACM International Conference Proceeding Series. T. 382. Association for Computing Machinery, 2009.](10.1145/1553374.1553511)

[Wikipedia](https://en.wikipedia.org/wiki/Mutual_information)

## Code 

[Scikit](https://scikit-learn.org/stable/modules/clustering.html#mutual-information-based-scores)
