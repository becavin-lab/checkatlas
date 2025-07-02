# Mean-squared Error 

## Description 

The Mean Squared Error (MSE) is a fundamental metric used in regression analysis and machine learning to evaluate the performance of predictive models. 
It measures the average of the squares of the differences between predicted values and actual values in a dataset. 
In statistics, MSE is defined as a risk function corresponding to the expected value of the squared error loss. 
By squaring the differences, MSE gives higher weight to larger errors, making it sensitive to outliers. 
MSE incorporates both the variance of the estimator (how widely spread the estimates are) and its bias (how far off the average estimated value is from the true value). 
A lower MSE value indicates that the model's predictions are closer to the actual values, reflecting better overall performance.


## Formulas 

### *Basic formula* : 

$$
MSE = \frac{1}{n} \displaystyle\sum_{i=1}^{n} (y_i - \hat{y}_i)^2
$$

Where:

- $n$ = number of data points

- $y_i$ = observed values

- $\hat{y}_i$ = predicted values

- $i$ = index from 1 to n

### *Matrix Notation* : 

$$MSE = \frac{1}{n} \left || y - \hat{y} \right || ^2$$


Where $y$ is the vector of observed values and $\hat{y}$ is the vector of predicted values.

### *MSE for an Estimator (Theoretical Definition)* : 

$$MSE(\hat{\theta}) = E[(\hat{\theta} - \theta)^2]$$

Where $\hat{\theta}$ is the estimator and $\theta$ is the true parameter.

### **Bias-Variance Decomposition:**

$$MSE = \text{Var}(\hat{\theta}) + [\text{Bias}(\hat{\theta})]^2 + \sigma^2$$

Where:

- $\text{Var}(\hat{\theta})$ = variance of the estimator

- $\text{Bias}(\hat{\theta})$ = bias of the estimator

- $\sigma^2$ = irreducible variance (noise)


## Sources 

[Open Problems](https://openproblems.bio/results/denoising?version=v1.0.0)

[Batson, J., Royer, L., & Webber, J. (2019). *Molecular cross-validation for single-cell RNA-seq*. bioRxiv.](https://doi.org/10.1101/786269)

[Wikipedia](https://en.wikipedia.org/wiki/Mean_squared_error)

[Encord Computer Vision Glossary](https://encord.com/glossary/mean-square-error-mse/#:~:text=In%20the%20fields%20of%20regression,target%20values%20within%20a%20dataset.)

## Code 

[Scikit](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mean_squared_error.html)
