# Normalized Mutual Information (NMI)

## Description 

Normalized Mutual Information (NMI) is a clustering evaluation metric used to compare predicted labels (i.e. after data integration) with ground-truth labels (i.e. known cell types). 
It quantifies the amount of shared information between two partitions, normalized to account for the complexity and size of the clusters. 
NMI is particularly useful in single-cell data integration to assess whether biological structure is preserved after batch correction.

## Formulas 

The NMI between two partitions $U$ and $V$ is defined as:

$$ NMI(U,V)=\frac{2 \cdot MI(U;V)}{H(U) + H(V)}$$

Where:
- $U$: clustering result (predicted labels)
- $V$: ground-truth labels
- $MI(U;V)$: mutual information between $U$ and $V$
- $H(U)$ and $H(V)$: entropies of $U$ and $V$

> Reminder
>
> The mutual information is computed as:
> $$ MI(U;V) = \displaystyle\sum_{u \in U} \displaystyle\sum_{v \in V} P(u, v) \log \left( \frac{P(u, v)}{P(u) P(v)} \right)$$
> 
> Where:
> - $P(u, v)$ is the joint probability of an element being in cluster $u$ and class $v$
> - $P(u)$ and $P(v)$ are the marginal probabilities


NMI ranges from 0 to 1:
- 0 : no mutual information (random clustering)
- 1 : perfect match between predicted and true labels

## Sources 

[OpenProblems.](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[Luecken, M.D. et al. *Benchmarking atlas-level data integration in single-cell genomics*. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-)

## Code 

[SCIB](https://github.com/theislab/scib/blob/main/scib/metrics/nmi.py)
