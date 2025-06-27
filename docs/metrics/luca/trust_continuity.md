# Trustworthiness and continuity

## Description 

These two metrics assess the quality of a dimensionality reduction by evaluating how well the data structure is preserved.
- Trustworthiness measures local fidelity : if two points are close in the low-dimensional space, they should also be close in the original space.
- Continuity does the reverse: it checks whether original neighbors remain close after projection. 

Both metrics range from 0 to 1. A value close to 1 indicates that the local structure is well preserved.
Trustworthiness penalizes intrusions (false neighbors added), while Continuity penalizes extrusions (true neighbors lost).
They are often used together to get a full picture of neighborhood preservation in embedding methods

## Formulas 

Let:
- $n$: the number of data points  
- $X \subset \mathbb{R}^p$: the high-dimensional input space  
- $Y \subset \mathbb{R}^q$: the low-dimensional embedded space  
- $r(i, j)$: the rank of point $j$ in terms of distance from point $i$ in the original space  
- $\hat{r}(i, j)$: the rank of point $j$ in terms of distance from point $i$ in the embedded space  
- $N_i^k(X)$: the set of the $k$ nearest neighbors of point $i$ in $X$  
- $N_i^k(Y)$: the set of the $k$ nearest neighbors of point $i$ in $Y$

### *Trustworthiness* : 

This metric penalizes neighbors in the low-dimensional space that were not true neighbors in the original space.

$$T(k) = 1 - \frac{2}{n k (2n - 3k - 1)} \sum_{i=1}^{n} \sum_{j \in N_i^k(Y)} \max(0,\; r(i, j) - k)$$

A high value means that most nearest neighbors in the embedded space were already close in the original space.

### *Continuity* : 

This metric penalizes neighbors in the original space that are lost in the low-dimensional space.

$$C(k) = 1 - \frac{2}{n k (2n - 3k - 1)} \sum_{i=1}^{n} \sum_{j \in N_i^k(X)} \max(0,\; \hat{r}(i, j) - k)$$

A high value indicates that original neighbors remain close after projection.

## Sources 

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

Lisha Chen & Andreas Buja. Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis. JASA, 2009.

Stasis et al. Semantically Controlled Adaptive Equalisation in Reduced Dimensionality Parameter Space. Applied Sciences, 2016.

https://scikit-learn.org/stable/modules/generated/sklearn.manifold.trustworthiness.html

## Code 

https://scikit-learn.org/stable/modules/generated/sklearn.manifold.trustworthiness.html
