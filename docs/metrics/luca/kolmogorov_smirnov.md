# Kolmogorov–Smirnov Statistic

## Description 

The Kolmogorov–Smirnov (KS) statistic measures the maximal difference between two empirical cumulative distribution functions (CDFs). 
In single-cell analysis, it is often used as a two-sample test to quantify separation between two distributions (e.g. inter- vs intra-type cell distance distributions). 
Intuitively, $D$ lies between 0 and 1 : a value near 0 implies the distributions largely overlap, while larger values indicate greater divergence. 
Pachter et al. use the two-sample KS statistic to quantify how well clusters or cell types separate in embedded vs ambient spaces (high $D$ means less overlap of the two distance distributions).


## Formulas 

### **Two-sample KS** : 

For two independent samples of sizes $n$ and $m$, with empirical CDFs $F_{1,n}(x)$ and $F_{2,m}(x)$, the KS statistic is  
 
$$
D_{n,m} = \sup_x \left|F_{1,n}(x) - F_{2,m}(x)\right|,
$$  
  
i.e., the maximum vertical distance between the two CDFs .  

### **One-sample KS** : 

Given a sample CDF $F_n(x)$ and a reference CDF $F(x)$, the one-sample KS statistic is  

$$
D_n = \sup_x \left|F_n(x) - F(x)\right|,
$$  

measuring the largest difference between the sample CDF and the reference distribution .  

### **Interpretation** : 

In both cases, $\sup_x$ denotes the supremum (maximum) over all real $x$.  
The KS statistic ranges from 0 (identical distributions) up to 1 (completely non-overlapping), with higher $D$ indicating stronger separation of the two distributions.

## Sources 

[Wikipedia](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#:~:text=The%20Kolmogorov%E2%80%93Smirnov%20statistic%20quantifies,distribution%20functions%20of%20two%20samples.)

[Chari & Pachter (2023), PLOS Comp. Biol. – Methods (KS statistic used to measure separation of distance distributions).](10.1371/journal.pcbi.1011288)


## Code 

[Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html)
