# PCR (Principal Component Regression)

## Description 

The PCR metric evaluates how well continuous biological variation is preserved after data integration in single-cell genomics.
It is particularly useful for assessing whether continuous structures—such as cell cycle progression or developmental trajectories—remain intact after batch effect correction.
The method involves performing a linear regression between the principal components of the integrated data and a known continuous biological variable (i.e. cell cycle score).


## Formulas 

Principal Component Regression (PCR) combines Principal Component Analysis (PCA) with linear regression:
- PCA step : Apply PCA to the integrated data matrix $X$ (cells × genes) to obtain a set of orthogonal principal components (PCs), denoted $Z = XW$, where $W$ is the matrix of eigenvectors.
- Component selection: Select the top $k$ PCs that capture most of the variance.
- Regression step : Fit a linear regression model using the selected PCs $Z_k$ as predictors for a continuous biological variable $y$ (e.g., cell cycle score).
- Metric score : Compute the coefficient of determination $R^2$ from the regression:

$$R^2 = 1 - \frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{\sum_{i=1}^n (y_i - \bar{y})^2}$$

Where:
- $y_i$ is the true value of the biological variable for cell $i$,
- $\hat{y}_i$ is the predicted value from the regression,
- $\bar{y}$ is the mean of all $y_i$ values.

A score close to 1 indicates strong preservation of the biological signal in the integrated space.

## Sources 

[OpenProblems.](https://openproblems.bio/results/batch_integration?version=v2.0.0)

[Luecken, M.D. et al. *Benchmarking atlas-level data integration in single-cell genomics*. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-)

[Wikipedia](https://en.wikipedia.org/wiki/Principal_component_regression)
## Code 

[SCIB](https://github.com/theislab/scib/blob/main/scib/metrics/pcr.py)
