# Density-Based Clustering Validation

## Description 

DBCV (Density-Based Clustering Validation) is a metric designed to evaluate the quality of clustering solutions (particularly for density-based clustering algorithms).
This metric evaluates clustering quality based on the relative density connection between pairs of objects. DBCV is particularly suited for identifying concave and nested clusters.
The DBCV index evaluates clustering quality by assessing density-based separation between clusters and cohesion within clusters.
Its values range between -1 and +1, with DBCV being a maximization index, where higher values correspond to better partitions.


## Formulas 

### *Step 1* :

For each group, calculate its "internal density" 

$\Rightarrow$ We look at how tightly the points are packed within the group.

### *Step 2* :

Calculate the "separation between groups"

$\Rightarrow$ We measure the empty space between different groups

### *Step 3* :

Final formula: DBCV = Mean of all group scores

$\Rightarrow$ Each group score = (Group separation) - (Internal density of the group)

$\Rightarrow$ The greater the separation AND the stronger the internal density $\rightarrow$ better score.

## Sources 

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

[Wikipedia](https://en.wikipedia.org/wiki/Density-based_clustering_validation)

[Davoud MOU"Christopherjenness/DBCV"LAVI et al. Density-Based Clustering Validation. Rapp. tech.](https://doi.org/10.1137/1.9781611973440.96)

## Code

Wikipedia $\rightarrow$ Github : "Christopherjenness/DBCV" 
