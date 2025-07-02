# Average Jaccard Distance

## Description 

The Average Jaccard Distance (AJD) is a metric introduce to quantify the distortion introduced by dimensionality reduction in single-cell RNA-seq (scRNA-seq) data. 
It measures how well local neighborhoods of cells are preserved after projecting the data into a lower-dimensional space. 
For each cell, the set of its $k$ nearest neighbors is computed both in the original (high-dimensional) space and in the reduced space. 
The Jaccard distance is then calculated between these two sets, and the AJD is defined as the mean of these distances over all cells. 
A low AJD indicates that local neighborhoods are well preserved (low distortion), while a high AJD implies substantial structural changes caused by the dimensionality reduction.


## Formulas

### Jaccard Distance  
Given two sets $A$ and $B$, the **Jaccard distance** is defined as :

$$
d_J(A, B) = 1 - \frac{|A \cap B|}{|A \cup B|} = \frac{|A \cup B| - |A \cap B|}{|A \cup B|}
$$

where $|A \cap B|$ is the number of elements in the intersection and $|A \cup B|$ is the number in the union.  
- If $A = B$, then $d_J = 0$ (no distortion).
- If $A \cap B = \emptyset$, then $d_J = 1$ (complete distortion).

### Average Jaccard Distance (AJD)  
Let $N$ be the total number of cells in the dataset. For each cell $i$, let $A_i$ and $B_i$ be the sets of its $k$ nearest neighbors in the original and reduced space, respectively. The Average Jaccard Distance is given by :

$$
\mathrm{AJD} = \frac{1}{N} \sum_{i=1}^N d_J(A_i, B_i)
$$

### Value Range and Interpretation  
Since $0 \leq d_J(A, B) \leq 1$, it follows that $0 \leq \mathrm{AJD} \leq 1$.  
- AJD $= 0$ indicates that local neighborhoods are perfectly preserved (ideal embedding).  
- AJD $= 1$ indicates that all local neighborhoods have changed entirely (maximum distortion).

Intermediate values represent partial preservation of neighborhood structures.

## Sources 

[Cooley S. M. et al. (2021), *“A novel metric reveals previously unrecognized distortion in dimensionality reduction of scRNA-seq data”*.](https://doi.org/10.1101/689851)  

[Wikipedia](https://en.wikipedia.org/wiki/Jaccard_index)

## Code 

[Scikit](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.jaccard_score.html)
