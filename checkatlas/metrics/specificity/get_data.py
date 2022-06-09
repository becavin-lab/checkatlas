"""
Created on Thu Jun  4 16:38:47 2020

@author: antoinecollin
"""

from typing import Literal

import anndata
import numpy as np
import pandas as pd

try:
    import compute  # Prevents circular import
except ImportError:
    from . import compute


def mean_celltype(
    adata: anndata, partition_key: str = "CellType", add_adata: bool = True
):
    """
    Computes the average gene expression by celltypes

    Parameters
    ----------
    adata
        Annotated data matrix.

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    gene_id_key
        The key in adata.obs corresponding to the gene ID column. Default
        will use adata.var.index.

    add_adata
        Indicate if the average matrix should be added to the varm field
        and its index to the uns field of adata.

    Returns
    -------
    average_by_celltype
        The mean expression by celltype matrix with format celltypes x genes.
    """
    celltypes = adata.obs[partition_key].cat.categories
    average_by_celltype = pd.DataFrame([], columns=list(adata.var.index))
    i = 0
    idx = []
    for cell in celltypes:
        if sum(adata.obs[partition_key] == cell) != 0:
            reduced_adata = adata[adata.obs[partition_key] == cell, :]
            mean_expr = np.asarray(reduced_adata.X.mean(axis=0))
            mean_expr = mean_expr.flatten()
            average_by_celltype.loc[i] = mean_expr
            idx.append(cell)
            i += 1
    average_by_celltype.index = idx
    if add_adata:
        adata.varm[f"ave_celltype_counts_{partition_key}"] = np.array(
            average_by_celltype.transpose()
        )
        adata.uns[
            f"ave_celltype_index_{partition_key}"
        ] = average_by_celltype.index
    return average_by_celltype


def get_average_celltype_counts(adata, partition_key: str = "CellType"):
    """
    Gets the mean expression by celltype matrix of adata. If it's already
    in the adata object, fetches it. If it's not, computes it and adds it
    to the adata object in varm with labels in uns

    Parameters
    ----------
    adata
        Annotated data matrix.

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    gene_id_key
        The key in adata.obs corresponding to the gene ID column. Default
        will use adata.var.index.

    Returns
    -------
    average_by_celltype
        The mean expression by celltype matrix.
    """
    try:
        adata.varm[f"ave_celltype_counts_{partition_key}"]
    except KeyError:
        average_by_celltype = mean_celltype(
            adata, partition_key=partition_key, add_adata=True
        )
    else:
        average_by_celltype = pd.DataFrame(
            adata.varm[f"ave_celltype_counts_{partition_key}"].transpose(),
            columns=adata.var.index,
        )
        average_by_celltype.index = adata.uns[
            f"ave_celltype_index_{partition_key}"
        ]
    return average_by_celltype


def get_anndata(adata_filename: str):
    """
    Fetches the anndata file

    Parameters
    ----------
    adata_filename

    Returns
    -------
    The Anndata object

    """
    adata = anndata.read_h5ad("")
    return adata


def get_markers(markers_filename: str):
    """
    Fetches the markers file

    Parameters
    markers_filename


    ----------
    markers_filename

    Returns
    -------
    A dict containing the markers list with their corresponding celltype

    """
    # markers = dict(csv.read(...))
    # return markers


def get_spe(
    adata, spe_metric: Literal["shannon", "tau", "gini"], partition_key
):
    specs = {
        "shannon": compute.shannon_average,
        "tau": compute.tau_average,
        "gini": compute.gini_average,
    }
    try:
        adata.var[f"{spe_metric}_{partition_key}"]
    except KeyError:
        spe_func = specs[spe_metric]
        spe_list = spe_func(adata, partition_key)
        adata.var[f"{spe_metric}_{partition_key}"] = spe_list
    else:
        spe_list = adata.var[f"{spe_metric}_{partition_key}"]
    return spe_list
