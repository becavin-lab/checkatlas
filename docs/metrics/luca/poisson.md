# Poisson Loss 

## Description 

The Poisson Loss is a cost function used to model count data, typically when observations follow a Poisson distribution. 
Unlike Mean Squared Error (MSE), it accounts for the intrinsic variance of count data. 
It is used in models like **Poisson regression** and in denoising benchmarks for scRNA-seq data. 
This loss evaluates the quality of reconstruction by comparing observed and predicted values while respecting the discrete and heteroscedastic nature of the data.

## Formulas 

The Poisson Loss between an observed value $y$ and a prediction $\hat{y}$ is given by :

$$
\mathcal{L}_{\text{Poisson}}(y, \hat{y}) = \hat{y} - y \log(\hat{y}) + \log(y!)
$$

Where :

- $y$: observed value (non-negative integer)

- $\hat{y}$: predicted value (strictly positive real number)

- $\log(y!)$ is constant with respect to $\hat{y}$ and is often ignored during optimization.

The loss is minimized when $\hat{y} \approx y$, and it penalizes errors more heavily for larger values of $y$, which aligns with the increasing variance of the Poisson distribution.

### *Value Range* :

- The Poisson Loss is non-negative : $\mathcal{L} \geq 0$

- It tends to infinity as $\hat{y} \to 0$ or when $\hat{y}$ is far from $y$

## Sources
[OpenProblems](https://openproblems.bio/results/denoising?version=v1.0.0)

[Metric Poisson](https://haibal.com/documentation/metric-poisson/)

[Wikipedia ](https://en.wikipedia.org/wiki/Poisson_regression)
  
[Batson, J., Royer, L., & Webber, J. (2019). *Molecular cross-validation for single-cell RNA-seq*. bioRxiv.](https://doi.org/10.1101/786269)

## Code 

[Scikit-learn](https://scikit-learn.org/stable/auto_examples/linear_model/plot_poisson_regression_non_normal_loss.html)

