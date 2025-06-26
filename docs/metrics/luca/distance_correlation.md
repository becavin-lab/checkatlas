# Distance Correlation 

## Description 

Distance Correlation is a statistical metric that quantifies both linear and nonlinear dependence between two random vectors. 
Unlike Pearson correlation, which only captures linear relationships, distance correlation equals zero if and only if the variables are truly independent. 
This makes it especially useful in bioinformatics for evaluating how well dimensionality reduction methods preserve complex relationships between features or samples.

## Formulas 


Let $(X_k, Y_k)$, for $k = 1, 2, \dots, n$, be a sample of paired observations.

*Compute pairwise Euclidean distances* :

$$a_{j,k} = \|X_j - X_k\|,\quad b_{j,k} = \|Y_j - Y_k\|$$

*Double-center the distance matrices* :

$$A_{j,k} = a_{j,k} - \bar{a}_{j.} - \bar{a}_{.k} + \bar{a}_{..}$$

$$B_{j,k} = b_{j,k} - \bar{b}_{j.} - \bar{b}_{.k} + \bar{b}_{..}$$

*Compute distance covariance* :

$$\text{dCov}^2(X, Y) = \frac{1}{n^2} \sum_{j,k} A_{j,k} B_{j,k}$$

*Compute distance variances*:

$$\text{dVar}^2(X) = \frac{1}{n^2} \sum_{j,k} A_{j,k}^2,\quad \text{dVar}^2(Y) = \frac{1}{n^2} \sum_{j,k} B_{j,k}^2$$

*Distance correlation*:

$$\text{dCor}(X, Y) = \frac{\text{dCov}(X, Y)}{\sqrt{\text{dVar}(X) \cdot \text{dVar}(Y)}}$$

*Value Range*
- The metric ranges from 0 (independence) to 1 (perfect dependence).
- A value of 0 implies true independence, not just lack of linear correlation.


## Sources 

[Schober, P., Boer, C., & Schwarte, L. A. (2018). Correlation coefficients. Anesthesia & Analgesia, 126(5), 1763–1768.](https://doi.org/10.1213/ane.0000000000002864)

[Coifman, R. R., & Lafon, S. (2006). Diffusion maps. Applied and Computational Harmonic Analysis, 21(1), 5–30.](https://doi.org/10.1016/j.acha.2006.04.006)

[Wikipedia](https://en.wikipedia.org/wiki/Distance_correlation)

## Code 

[SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.correlation.html)
