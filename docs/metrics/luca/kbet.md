# kBET (k-nearest neighbor Batch Effect Test)

## Description 

kBET is a statistical metric used to evaluate the effectiveness of batch effect correction in single-cell transcriptomics data. 
It is based on the idea that, in a well-integrated space, the $k$ nearest neighbors of a cell should be randomly distributed across batches. 
If batch labels remain distinguishable after integration, this indicates poor correction. 
kBET tests the null hypothesis that the batch distribution in a cell’s local neighborhood matches the global batch distribution. 
A high rejection rate of this hypothesis indicates poor integration.

## Formulas 

- Null hypothesis : The batch distribution among the $k$ nearest neighbors of a cell matches the global batch distribution.
- For each cell:
   - Identify its $k$ nearest neighbors.
   - Count the frequency of batch labels in this neighborhood.
   - Compare this local distribution to the global batch distribution using a chi-squared test.

The rejection rate is the proportion of cells for which the null hypothesis is rejected (typically at $p < 0.05$).


The final score is:

$$\text{kBET score} = 1 - \text{rejection rate}$$

A score close to 1 indicates good integration (few rejections), while a score close to 0 indicates poor integration.

## Sources 



## Code


