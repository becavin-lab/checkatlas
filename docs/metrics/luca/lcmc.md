### Local Continuity Meta Criterion (LCMC)

## Description 

The Local Continuity Meta Criterion (LCMC) evaluates the quality of dimensionality reduction methods by measuring how well local neighborhood structures are preserved between the original high-dimensional space and the low-dimensional embedding.
For each data point, it compares the set of its K nearest neighbors before and after projection. The metric quantifies the overlap and adjusts for chance-level matches.
LCMC is particularly useful for identifying how well local geometry is maintained, complementing other metrics like Trustworthiness and Continuity.

## Formulas 

Let:
- $N$ = number of data points
- $K$ = neighborhood size
- $\mathcal{N}^K(i)$ = set of the $K$ nearest neighbors of point $i$ in the original space
- $\nu^K(i)$ = set of the $K$ nearest neighbors of point $i$ in the projected space

Then the LCMC is defined as :

$$
\mathrm{LCMC}(K) = \frac{1}{N K} \sum_{i=1}^{N} \left| \mathcal{N}^K(i) \cap \nu^K(i) \right| - \frac{K}{N - 1}
$$

This means:
- The first term is the average number of overlapping neighbors (shared between both spaces), normalized.
- The second term $\frac{K}{N-1}$ is the expected overlap by chance.


This adjustment means that:
- LCMC = 1 implies perfect local continuity (all true neighbors preserved),
- LCMC ≈ 0 implies performance similar to random neighbor selection,
- LCMC < 0 suggests worse-than-random behavior.

## Sources 

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

Jarkko VENNA et Samuel KASKI. « Local multidimensional scaling ». In : Neural Networks (2006).

John A. Lee. « Quality assessment of nonlinear dimensionality reduction based on K-ary neighborhoods ».

https://search.r-project.org/CRAN/refmans/coRanking/html/LCMC.html

## Code 

https://search.r-project.org/CRAN/refmans/coRanking/html/LCMC.html

https://github.com/hj-n/zadu
