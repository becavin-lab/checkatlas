### Specificity metrics 

| Metric                                                                 | Description                                                                 | Boundary             | Optimum              |
|------------------------------------------------------------------------|-----------------------------------------------------------------------------|----------------------|----------------------|
| [Gini coefficient](gini.md)                                           | Evaluates the inequality of the gene's distribution.                        | –                    | The lower the better |
| [Kendall's tau](tau.md)                                               | Similar to Gini and Shannon, but less restrictive.                          | [0;1] or [-1;1]      | –                    |
| [One-v-All](ova_ovm.md)                                              | Gene expression in one cell type versus the rest.                           | [0,1]                | 1                    |
| [Onevs-Max](ova_ovm.md)                                              | Gene expression in one cell type vs the second more highly expressed cell. | [0,+∞[               | –                    |
| [Shannon Entropy](shannon.md)                                         | Entropy of gene distribution across cell types.                             | [0;1]                | –                    |
