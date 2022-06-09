"""
Created on Thu Jun  4 16:38:47 2020

@author: antoinecollin
"""

from typing import Collection, Optional

import anndata
import numpy as np
import pandas as pd
from get_data import get_average_celltype_counts
from scipy.stats import entropy


def one_v_all_matrix(adata: anndata, partition_key: str = "CellType"):
    """
    Computes the one_v_all specificity of all genes for all celltype.
    It is defined as the average expression of a gene in a celltype
    divided by the sum of mean expressions in every celltypes.

    Parameters
    ----------
    adata : anndata
        Annotated data matrix.

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    spe_one_v_all : DataFrame
        DataFrame of shape n_genes x n_celltype containing the one_v_all
        specificities of all genes by celltypes

    """
    average_by_celltype = get_average_celltype_counts(
        adata=adata, partition_key=partition_key
    )
    spe_matrix = pd.DataFrame(
        columns=list(average_by_celltype.columns),
        index=list(average_by_celltype.index),
    )
    for celltype in average_by_celltype.index:
        spe_matrix.loc[celltype] = average_by_celltype.loc[
            celltype
        ] / average_by_celltype.sum(axis=0)
    return spe_matrix


def one_v_max_matrix(adata: anndata, partition_key: str = "CellType"):
    """
    Computes the one_v_max specifity of all genes for the celltypes which
    maximizes it. For a gene, one_v_max specificity is defined as the ratio
    between the average expression of this gene in the most expressed
    celltype divided by the second most expressed celltype.

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with the
        value of the specificity and the specific celltype name.

    """
    try:
        adata.uns[f"ovm_{partition_key}"]
    except KeyError:
        average_by_celltype = get_average_celltype_counts(
            adata=adata, partition_key=partition_key
        )
        spe_matrix = pd.DataFrame(
            columns=list(average_by_celltype.columns),
            index=list(average_by_celltype.index),
        )
        for celltype in average_by_celltype.index:
            spe_matrix.loc[celltype] = average_by_celltype.loc[
                celltype
            ] / average_by_celltype.loc[
                average_by_celltype.index != celltype
            ].max(
                axis=0
            )
        adata.uns[f"ovm_{partition_key}"] = spe_matrix
    else:
        spe_matrix = adata.uns[f"ovm_{partition_key}"]
    return spe_matrix


def shannon_average(adata: anndata, partition_key: str = "CellType"):
    """
    Computes the Shannon entropy of every gene based on the distribution of
    their average value across celltypes.

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with the
        value of the specificity and the specific celltype name.

    """
    prob_matrix = one_v_all_matrix(adata=adata, partition_key=partition_key)
    spe_shannon = pd.Series(
        index=prob_matrix.columns, name="Shannon_entropies"
    )
    for gene in prob_matrix.columns:
        spe_shannon[gene] = entropy(list(prob_matrix[gene]))
    return spe_shannon


# TODO
def shannon_sample(adata: anndata, genes: Optional[Collection[str]]):
    """
    Computes the Shannon entropy of the gene distribution across the
    celltypes based on a subsample of the expression matrix. The sampling
    has a normalizing effect, preventing the most represented celltypes
    to overshadow the rare ones.

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    genes : Collection[str]
        List of genes for which to compute the specificities. If not
        specified, it is computed for all the genes.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with
        the value of the specificity and the specific celltype name.


    """


def tau_average(adata: anndata, partition_key: str = "CellType"):
    """
    Computes tau coefficient of every gene based on the distribution of
    their average value across celltypes.

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with the
        value of the specificity and the specific celltype name.

    """
    average_by_celltype = get_average_celltype_counts(
        adata=adata, partition_key=partition_key
    )
    ratio_matrix = pd.DataFrame(
        columns=list(average_by_celltype.columns),
        index=list(average_by_celltype.index),
    )
    for celltype in average_by_celltype.index:
        ratio_matrix.loc[celltype] = average_by_celltype.loc[
            celltype
        ] / average_by_celltype.max(axis=0)
    inverted_ratio = 1 - ratio_matrix
    return inverted_ratio.sum(axis=0) / (len(average_by_celltype.index) - 1)


# TODO
def tau_sample(adata: anndata, genes: Optional[Collection[str]]):
    """

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    genes : Collection[str]
        List of genes for which to compute the specificities. If not
        specified, it is computed for all the genes.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with
        the value of the specificity and the specific celltype name.

    """


def gini_average(adata: anndata, partition_key: str = "CellType"):
    """
    Computes gini coefficient of every gene based on the distribution of
    their average value across celltypes.

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with
        the value of the specificity and the specific celltype name.

    """
    average_by_celltype = get_average_celltype_counts(
        adata=adata, partition_key=partition_key
    )
    spe_gini = pd.Series(
        index=average_by_celltype.columns, name="gini_Coefficient"
    )
    n = average_by_celltype.shape[0]
    index = np.arange(1, n + 1)
    for gene in average_by_celltype.columns:
        gini_list = average_by_celltype.loc[:, gene].sort_values()
        spe_gini[gene] = (np.sum((2 * index - n - 1) * gini_list)) / (
            n * np.sum(gini_list)
        )
    return spe_gini


# TODO
def gini_sample(adata: anndata, genes: Optional[Collection[str]]):
    """

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    genes : Collection[str]
        List of genes for which to compute the specificities. If not
        specified, it is computed for all the genes.

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with
        the value of the specificity and the specific celltype name.

    """


# TODO
def dist_gene_to_celltype(adata):
    """
    Just from the adata, identifies the marker genes and corresponding
    celltypes. Those are computed using the criteria of the fonction
    (threshold, metric etc...)

    Parameters
    ----------
    adata : anndata
        The corrected expression matrix

    Returns
    -------
    spe_one_v_max : dict
        Dictionary where keys are genes and attributes are lists with the
        value of the specificity and the specific celltype name.

    """


# To Move To Analysis


# if __name__ == "__main__" :
#     a = specificity_summary(adata=barbry,marker_genes =
#     mark_dict,partition_key='CellType')
