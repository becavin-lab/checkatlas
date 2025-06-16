### NOT READY YET !! TO REVIEW AND CHANGE SOME THINGS 


### One‑vs‑All (OvA) & One‑vs‑One (OvO) Classification

## Description 

In multiclass classification, these strategies decompose the task into multiple binary classifiers :
- One‑vs‑All (OvA) trains one binary classifier per class $k$, distinguishing class $k$ (positive) from all other classes (negative). At prediction, the class with the highest decision score is chosen.
- One‑vs‑One (OvO) trains a classifier for each pair of classes $(i,j)$, leading to $\frac{K(K-1)}{2}

Both methods allow applying binary classifiers to multiclass problems.
OvA is computationally lighter (linear in $K$) but may face class imbalance.
OvO provides better class separation but is more expensive (quadratic in $K$).

## Formulas 

*Notation*:
- $K$ is the number of classes
- Training data : $\{(x_i,y_i)\}, \text{ with } y_i \in \{1, 2, \dots, K\}$
- $f_k(x)$ : output of the classifier for class $k$
- $g_{i,j}(x)$ : output of the classifier distinguishing class $i$ vs class $j$.

*One-vs-All (OvA)*
For each class $k \in \{1,\dots,K\}$, train a binary classifier : 

$$f_k(x) =
\begin{cases}
+1 & \text{if } y = k \\
-1 & \text{if } y \ne k
\end{cases}$$

The predicted label is the one with the highest decision function :

$$\hat{y} = \arg\max_{k} f_k(x)$$

* One‑vs‑One (OvO)*
For each pair of classes $(i,j)$ such that $1 \leq i \leq j \leq K$, train :

$$g_{i,j}(x) =
\begin{cases}
+1 & \text{if } y = i \\
-1 & \text{if } y = j
\end{cases}$$ 

Each classifier votes for a class. The predicted label is the one with the most votes :

$$\hat{y} = \arg\max_{c} \sum_{i < j} \mathbf{1}(g_{ij}(x) = +1 \text{ and } c = i)$$


## Sources 

https://en.wikipedia.org/wiki/Multiclass_classification#One-vs.-rest

https://en.wikipedia.org/wiki/Multiclass_classification#One-vs.-one

Skene NG et al. (2018) Genetic identification of brain cell types underlying schizophrenia. Nature Genetics, 50(6):825–833.

“Applying Deep Learning algorithm to perform lung cells annotation”, A. Collin, 2020

## Code 


