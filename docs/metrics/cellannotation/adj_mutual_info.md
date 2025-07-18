# Adjusted Mutual Information

## Description

Adjusted Mutual Information (AMI) is a clustering evaluation metric that measures the agreement between two label assignments, corrected for chance. 
Unlike Normalized Mutual Information (NMI), AMI accounts for the fact that random clusterings can have non-zero mutual information. 
It is especially useful when comparing clustering algorithms or evaluating cluster quality against a ground truth. 
AMI is bounded : its maximum is 1 (perfect match), and its expected value is close to 0 when clusterings are randomly assigned, even if the number of clusters varies.

### *Reminder* : 
[Mutual Information](cellannotation/mutual_information.md)

## Formulas : 
The Adjusted Mutual Information is then defined as :

$$
AMI(U, V) = \frac{MI(U, V) - \mathbb{E}[MI(U, V)]}{\max(H(U), H(V)) - \mathbb{E}[MI(U, V)]}
$$

- $H(U)$ and $H(V)$ are the entropies of the labelings.
- $\mathbb{E}[MI(U, V)]$ is the expected mutual information of random clusterings with the same size distributions.

This normalization ensures that $AMI(U, V) = 1$ when the clusterings are identical, and close to 0 if the similarity is no better than random.

## Sources : 

[Wikipedia : Adjusted Mutual Information](https://en.wikipedia.org/wiki/Adjusted_mutual_information)  

[Scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html)

## Code : 

[Scikit-learn documentation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html)
