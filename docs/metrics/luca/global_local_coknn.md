#

## Description 


The Global & Local co-KNN metric evaluates the quality of dimensionality reduction by comparing neighborhood relationships between the original high-dimensional space and the reduced low-dimensional space. 
It combines two perspectives :
- Global KNN : assesses how well the global structure (distant neighbors) is preserved.
- Local KNN : evaluates the preservation of local neighborhood relationships.

The co-KNN metric merges these two aspects to provide a balanced view of both local and global structure preservation after dimensionality reduction.

## Formulas 

Let:
- $N_k^{\text{high}}(i)$ be the set of the $k$ nearest neighbors of sample $i$ in the original space.
- $N_k^{\text{low}}(i)$ be the set of the $k$ nearest neighbors of sample $i$ in the reduced space.

Then:
- Local KNN Recall :

$$\text{Local}(i) = \frac{|N_k^{\text{high}}(i) \cap N_k^{\text{low}}(i)|}{k}$$
 
- Global KNN Recall :

Similar to the local recall but computed with a larger $k$ to capture global structure.


- co-KNN Score:

$$\text{co-KNN} = \alpha \cdot \text{Local} + (1 - \alpha) \cdot \text{Global}$$

Where $\alpha \in [0, 1]$ is a weighting parameter

## Sources 

[Zhang, Y., Shang, Q., & Zhang, G. (2021). *pyDRMetrics - a python toolkit for dimensionality reduction quality assessment*. Heliyon, 7(2), e06199.](https://doi.org/10.1016/j.heliyon.2021.e06199)

[Anava, O., & Levy, K. Y. (2016). *k-Nearest Neighbors: From Global to Local*. NeurIPS.](https://papers.nips.cc/paper/6373-k-nearest-neighbors-from-global-to-local)

[Laguna, V., & Lopes, A. A. de A. (2010). *Combining local and global KNN with cotraining*. ECAI 2010.](https://doi.org/10.3233/978-1-60750-606-5-815)

[OpenProblems](https://openproblems.bio/results/dimensionality_reduction?version=v1.0.0)

[Ultralytics](https://www.ultralytics.com/glossary/k-nearest-neighbors-knn#:~:text=K%2DNearest%20Neighbors%20(KNN)%20is%20a%20fundamental%20algorithm%20in,for%20understanding%20instance%2Dbased%20learning.)

## Code
