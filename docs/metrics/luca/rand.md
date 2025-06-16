### RAND INDEX

## Description 

The Rand index is a statistical metric used to measure the similarity between two partitions of the same dataset, particularly for evaluating the quality of clustering algorithms.
It measures the frequency of occurrence of agreements over all pairs, or the probability that two partitions will agree on a randomly chosen pair.
This metric can be interpreted as the percentage of correct decisions made by the clustering algorithm.
The index ranges between 0 and 1, where 0 indicates that the two clusterings do not agree on any pair of points, and 1 indicates that the clusterings are exactly identical.
The metric works by counting four types of pairs : positive agreements ($a$), negative agreements ($b$), type 1 disagreements ($c$), and type 2 disagreements ($d$).
The Rand index intuitively represents the ratio between total agreements and the total number of possible pairs in the dataset.

## Formula

$$R=\frac{(a+b)}{(a+b+c+d)}=\frac{(a+b)}{\binom{2}{n}}$$

Where : 
- $a$ is the number of pairs of elements that are in the same subset in $X$ and in the same subset in $Y$.
- $b$ is the number of pairs of elements that are in different subsets in $X$ and in different subsets in $Y$.
- $c$ is the number of pairs of elements that are in the same subset in $X$ and in different subsets in $Y$.
- $d$ is the number of pairs of elements that are in different subsets in $X$ and in the same subset in $Y$.
- $\binom{2}{n}$  is the total number of possible pairs, defined by : $\binom{2}{n}=\frac{n(n-1)}{2}$

 ## Sources 

https://en.wikipedia.org/wiki/Rand_index

Lawrence HUBERT et Phipps ARABIE. « Comparing partitions ». In : Journal of Classification 2 (déc. 1985).

Applying Deep Learning algorithm to perform lung cells annotation, A. Collin

 ## Code

 https://scikit-learn.org/stable/modules/clustering.html#rand-index

 
