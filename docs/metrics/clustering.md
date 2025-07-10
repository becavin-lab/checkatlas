### Clustering metrics:


| Metric                                                   | Type                                                                 | Description                                                                                           | Boundary   | Optimum               |
|----------------------------------------------------------|----------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|------------|------------------------|
| [Davies-Bouldin](clustering/dbi.md)                                 | ![Static Badge](https://img.shields.io/badge/Internal-blue)          | Internal metric based on a similarity notion between clusters.                                         | –          | The lower the better   |
| [DBCV](clustering/dbcv.md)                                          | ![Static Badge](https://img.shields.io/badge/Internal-blue)          | Internal metric for evaluating density-based clustering methods.                                       | [-1, 1]    | 1  |
| [Dunn Index](clustering/dunn.md)                                    | ![Static Badge](https://img.shields.io/badge/Internal-blue)          | Trade-off between inter- and intra-cluster distances.                                                  | –          | The higher the better  |
| [Silhouette](clustering/silhouette.md)                              | ![Static Badge](https://img.shields.io/badge/Internal-blue)          | Trade-off between average inter- and intra-cluster distances.                                          | [-1, 1]    | 1  |
| [Kolmogorov–Smirnov](clustering/kolmogorov_smirnov.md)              | ![Static Badge](https://img.shields.io/badge/External-purple)        | Maximum difference between two empirical CDFs; close to 1 indicates highly dissimilar distributions.   | [0, 1]     | 1  |
