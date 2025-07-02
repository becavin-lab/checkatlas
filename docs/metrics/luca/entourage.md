# Entourage

## Description 

The Entourage metric is a measure developed to evaluate the local quality of dimensionality reduction, unlike Kruskal stress which evaluates global deformation. 
It compares, for each instance, its k nearest neighbors in the reference space (obtained by SVD) with its k nearest neighbors in the reduced space after transformation. 
This metric quantifies the preservation of the local structure of the point cloud by calculating the percentage of common neighbors preserved after dimensional reduction. 
The higher the Entourage value, the better the original local structure is preserved. 
It allows identification of whether dimensionality reduction maintains proximity relationships between points, which is crucial for interpreting complex biological data.

## Formulas 

For each instance $\bar{X}_i$:

- $$N^{\text{ref}}_i$$ = $$k$$ nearest neighbors in reference space  

- $$N^{\text{new}}_i$$ = $$k$$ nearest neighbors in reduced space  

- $$G_i = \text{card}(N^{\text{ref}}_i \cap N^{\text{new}}_i)$$ = number of common neighbors

The trustworthiness-like metric is defined as:

$$
Ent_k=\frac{\displaystyle\sum_{i=1}^{n} G_i}{G}
$$

Where:

- $G = nk$ (normalization constant)

So:

- $$\text{Ent} \in (0, 1)$$

Interpretation:

- $0$: no common neighbors preserved  

- $1$: all local neighbors perfectly preserved

### *NB* : What is SVD ? 

Singular Value Decomposition (SVD) is a matrix factorization technique that decomposes any rectangular matrix into three matrices. 

