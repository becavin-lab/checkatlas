import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector

"""
Module for management of Atlas n Seurat format

- Adata and sumarry tables are created from seuratà
- QC are created using Seurat
- UMAP, t-SNE are created using Seurat
- For creating metrics tables, obs tables are first converted from R to python and then scikitlearn is used.


"""


def install_library():
    # R package names
    packnames = ['Seurat']
    utils = rpackages.importr('utils')

    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1) # select the first mirror in the list

    # Selectively install what needs to be install.
    # We are fancy, just because we can.
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))