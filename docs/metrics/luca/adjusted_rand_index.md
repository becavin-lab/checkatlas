# Adjusted Rand Index (ARI)

## Description 

> « Adjusted Rand Index compares clustering overlap, correcting for random labels and considering correct overlaps and disagreements. »
[Open Problems](https://openproblems.bio/results/batch_integration?version=v2.0.0#luecken2022benchmarking)

The Adjusted Rand Index (ARI) is a clustering evaluation metric that measures similarity between two data partitions while accounting for chance correction. 
Unlike the classical Rand Index, ARI normalizes the score to eliminate the effect of random groupings.
It allows comparison of clusters obtained by different integration methods with reference annotations (cell types, biological conditions).
ARI ranges between $-1$ and $1$, where 0 represents chance-level agreement, 1 indicates perfect correspondence, and negative values signal disagreement worse than chance.


## Formulas 

Let $U=(U_1,U_2,\dots,U_n)$ and $V=(V_1,V_2,\dots,V_n)$ be two paritions of a set of $n$ elements. 

We define, 
- $n_{i,j}=\left | U_i \cap V_j \right |$ the number of elements in the intersection
- $a_i=\left | A_i \right |$ the size of cluster $i$ in $U$
- $b_j=\left | V_j \right |$ the size of cluster $j$ in $V$

> *Rand Index*
> 
> We remind that the Rand index (RI) is define as follows : 
> $$RI=\frac{a+b}{\binom{2}{n}}$$
> where,
> - $a$ is the number of pairs of elements in the same cluster in $U$ and in $V$
> - $b$ is the number of pairs of elements in different clusters in $U$ and in $V$
> - $\binom{2}{n}=\frac{n(n-1)}{2}$ is the total number of pairs

> [For further information](rand.md)

Thus, the adjusted Rand index is given by :

$$ARI = \frac{\displaystyle\sum_{ij} \binom{n_{ij}}{2} - \frac{\displaystyle\sum_i \binom{a_i}{2} \displaystyle\sum_j \binom{b_j}{2}}{\binom{n}{2}}}{\frac{1}{2}\left[\displaystyle\sum_i \binom{a_i}{2} + \displaystyle\sum_j \binom{b_j}{2}\right] - \frac{\displaystyle\sum_i \binom{a_i}{2} \displaystyle\sum_j \binom{b_j}{2}}{\binom{n}{2}}}$$

## Sources 

[Open Problems](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[Benchmarking atlas-level data integration in single-cell genomics, Leucken et al., 2021 ](https://doi.org/10.1038/s41592-021-01336-8)

[Rand Index](https://en.wikipedia.org/wiki/Rand_index)

## Code 

[Scikit](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html) 

