import anndata
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy.spatial import distance_matrix
from sklearn.cluster import dbscan, k_means
from sklearn.preprocessing import LabelEncoder

from .clust_plots import hist_internal

adata = anndata.read_h5ad(
    r"C:\Users\ipmc\Documents\Package_Metrics\scanpy\scanpy\datasets\
    HCA_Barbry_Hg19_seurat.h5ad"
)
adata.X = adata.X.todense()

sc.pp.neighbors(adata)

sc.tl.leiden(adata=adata)
sc.tl.louvain(adata=adata)

sc.pl.umap(adata, color=["CellType"])

hist_internal(
    adata=adata, partition_keys=["leiden", "louvain"], reference=["CellType"]
)

le = LabelEncoder()
le.fit(adata.obs["CellType"])
classes = le.classes_
n_c = len(classes)
# f, axes = plt.subplots(4,6)
i = 1
for type in classes:
    category = adata[adata.obs["CellType"] == type]
    sc.pl.umap(category, ax=plt.subplot(4, 6, i))
    plt.title(type)
    i += 1
sc.pl.umap(adata, color=["CellType"])


kmeans = k_means(adata.X, 28)
np.array([str(x) for x in kmeans[1]])
dbscan = dbscan(adata.X, 28)
adata.obs["dbscan"] = dbscan[1]


sc.pl.umap(adata, color=["kmeans"])
sc.pp.pca(adata, n_comps=20)

# phenograph.cluster(data = adata.obsm['X_pca'][:,:12])

phen_typing = sc.external.tl.phenograph(adata.obsm["X_pca"][:, :12])

x = adata.X.todense()

dist = distance_matrix(x[:10, :], x[:10, :])

# Create test set
brain = anndata.read_h5ad(
    r"C:\Users\ipmc\Documents\Package_Metrics\scanpy\scanpy\
    datasets\brain.h5ad"
)
sc.pl.highest_expr_genes(
    brain,
    n_top=20,
)

sc.pp.filter_cells(brain, min_genes=200)
sc.pp.filter_genes(brain, min_cells=3)

sc.pp.normalize_total(brain, target_sum=1e4)
sc.pp.log1p(brain)

sc.tl.pca(brain, svd_solver="arpack")

sc.pp.neighbors(brain, n_neighbors=10, n_pcs=40)

sc.tl.umap(brain)

sc.tl.louvain(brain)
sc.tl.leiden(brain)

sc.pl.umap(brain, color=["leiden", "louvain"])
sc.tl.tsne(brain)

brain.write_h5ad(
    r"C:\Users\ipmc\Documents\Package_Metrics\scanpy\scanpy\datasets\
    test_dataset.h5ad"
)

# Load test set

test_set = anndata.read_h5ad(
    r"C:\Users\ipmc\Documents\Package_Metrics\scanpy\scanpy\datasets\
    test_dataset.h5ad"
)
