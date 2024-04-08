### Add your own metrics

You nedd to blablabla ...


### Add in checkatlas

Create a python file in checkatlas folder :

checkatlas/metrics

Example : checkatlas/metrics/annot/rand_index.py
or checkatlas/metrics/cluster/silhouette.py

Add the following code for a clustering metric:

```python
def run(count_repr, annotations):
    metric = silhouette_score(count_repr, annotations)
    return metric
```

Add the following code for a cell annotation metric:

```python
def run(annotation, ref_annotation):
    metric = adjusted_rand_score(annotation, ref_annotation)
    return metric
```

Add the following code for a dimensionality reduction metric:

```python
def run(high_dim_counts, low_dim_counts):
    low_dim_distances = euclidean_distances(low_dim_counts, low_dim_counts)
    high_dim_distances = euclidean_distances(high_dim_counts, high_dim_counts)
    metric = np.sqrt(
        np.square(high_dim_distances - low_dim_distances).sum()
        / np.square(low_dim_distances).sum()
    )
    return metric
```

Reference your metric in the __init__.py in the __all__ variable.
Example for `checkatlas.metrics.cluster` module:

__all__ = ["silhouette", "davies_bouldin", "calinski_harabasz"]

It will be added automatically to the arguments of checkatlas workflow.

Arguments:

. --metric_cluster
. --metric_annot
. --metric_dimred


### Documentation

Add a file my_metric.md in oone of the folder and the file :
- cellannotation - cell_annotation.md
- clustering
- dimred


