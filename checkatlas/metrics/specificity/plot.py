"""
Created on Thu Jun  4 16:38:47 2020

@author: antoinecollin
"""

from typing import Collection, Iterable, Optional

import anndata
import get_data
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from compute import one_v_max_matrix
from scanpy.plotting._tools.scatterplots import _panel_grid

ROOTDIR = r"C:\Users\ipmc\Documents\Package_Metrics\scanpy\scanpy"
FILENAME = "HCA_Barbry_Hg19_seurat.h5ad"
filename = "\\".join([ROOTDIR, "datasets", FILENAME])
annick = anndata.read_h5ad(filename)


def gene_expr_distribution(
    adata: anndata,
    gene: str,
    ax,
    partition_key: str = "CellType",
    celltype: Optional[str] = None,
):
    """
    Plot the distribution of the average gene expression across celltypes

    Parameters
    ----------
    adata
        The corrected expression data.

    gene
        Gene to plot

    ax
        axe on which to plot

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    celltype
        Celltype to highlight in the rugplot.

    Returns
    -------

    """
    average_by_celltype = get_data.get_average_celltype_counts(
        adata=adata, partition_key=partition_key
    )
    out = sns.distplot(
        average_by_celltype[gene],
        kde=False,
        rug=True,
        ax=ax,
        hist_kws={
            "histtype": "stepfilled",
            "facecolor": "royalblue",
            "edgecolor": "black",
            "alpha": 0.8,
            "linewidth": 1,
        },
    )
    indent = ""
    for spe in ["shannon", "tau", "gini"]:
        spe_to_plot = round(
            get_data.get_spe(adata, spe_metric=spe, partition_key="CellType")[
                gene
            ],
            2,
        )
        out = plt.annotate(
            spe + " = " + str(spe_to_plot) + indent,
            xy=(0.79, 0.70),
            xycoords="axes fraction",
        )
        indent += "\n"
    if isinstance(celltype, str):
        out = sns.rugplot(
            [average_by_celltype.loc[celltype, gene]], color="red"
        )
        out = plt.annotate(
            celltype,
            xy=[average_by_celltype.loc[celltype, gene], 0.06],
            xycoords=("data", "axes fraction"),
            annotation_clip=False,
        )
    ax.set(xlabel=None, ylabel="Number of celltypes")
    ax.set_title(gene, fontsize=10)
    return out


def marker_genes_distribution(
    adata: anndata,
    gene_list: list,
    celltype: Optional[str],
    partition_key: str = "CellType",
    style="seaborn-white",
):
    """
    Display the distribution of the expressions of the genes across all
    celltypes
    Parameters
    ----------
    adata
        The corrected expression data.

    gene_list
        List of genes to plot, usually marker genes

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    celltype
        Celltype to highlight in the rugplot

    Returns
    -------

    """
    plt.style.use(style)
    n_panels = len(gene_list)
    f, grid = _panel_grid(
        hspace=0.25, wspace=0.25, ncols=4, num_panels=n_panels
    )
    for i in range(n_panels):
        gene = gene_list[i]
        ax = f.add_subplot(grid[i])
        gene_expr_distribution(
            adata=adata,
            gene=gene,
            ax=ax,
            partition_key=partition_key,
            celltype=celltype,
        )


def distrib_one_v_max(
    adata: anndata,
    celltype: str,
    ax,
    gene_highlight: Iterable[str],
    partition_key: str = "CellType",
):
    """

    Parameters
    ----------
    adata
        The corrected expression data.

    celltype
        Celltype to be plotted

    gene_highlight
        List of genes to highlight on the plot, usually marker genes

    partition_key
        The key in adata.obs corresponding to the annotations to be used.


    Returns
    -------

    """
    one_v_max = one_v_max_matrix(adata, partition_key=partition_key)
    to_plot = one_v_max.loc[celltype]
    to_plot_trunc = to_plot[
        to_plot > 1
    ]  # The only relevant one_v_max are the ones >1
    gene_highlight = list(set(gene_highlight) & set(list(to_plot_trunc.index)))
    if gene_highlight:
        to_plot_highlight = to_plot_trunc[gene_highlight]
        to_plot_trunc = to_plot_trunc.drop(gene_highlight)
        out = sns.stripplot(
            y=to_plot_trunc, ax=ax, orient="v", s=3, linewidth=0.25
        )
        out = sns.stripplot(
            y=to_plot_highlight, ax=ax, orient="v", color="red", linewidth=0.5
        )
        out.axes.set_xticks([])
        if len(gene_highlight) == 1:
            spe_ovm = round(to_plot_highlight[0], 2)
            out.axes.set_title(
                label=f"{gene_highlight[0]} \n "
                f"{to_plot_trunc.name} onevmax={spe_ovm}",
                loc="center",
                fontsize="small",
            )
        else:
            out.axes.set_title(
                label=f"{to_plot_trunc.name}", loc="center", fontsize="small"
            )
        out.axes.set_ylabel(ylabel="")
    else:
        out = sns.stripplot(
            y=to_plot_trunc, ax=ax, orient="v", s=3, linewidth=0.25
        )
        out.axes.set_title(
            label=f"{to_plot_trunc.name}", loc="center", fontsize="small"
        )
        out.axes.set_ylabel(ylabel="")
    return out


# Not used
def one_v_max_celltypes(
    adata: anndata,
    celltype_list: Optional[Collection[str]] = None,
    partition_key: str = "CellType",
):
    """
    Plots the one_v_max stripplot across the specified celltypes
    in the same plot


    Parameters
    ----------
    adata
        The corrected expression data.

    celltype_list
        celltypes to plot. Default plots every celltypes

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------

    """
    one_v_max = one_v_max_matrix(adata, partition_key=partition_key)
    one_v_max["celltype"] = one_v_max.index
    mean_melt = pd.melt(
        one_v_max,
        id_vars="celltype",
        value_vars=list(one_v_max.columns[:-1]),
        var_name="genes",
        value_name="expression",
    )
    mean_melt = mean_melt[mean_melt["expression"] > 1]
    out = sns.stripplot(
        x=mean_melt["celltype"], y=mean_melt["expression"], s=2
    )
    out.set(yscale="log")
    return out


# Used
def one_v_max_genelist(
    adata: anndata, gene_list: list, partition_key: str = "CellType"
):
    """
    For each gene, plots the one_v_max stripplot in the celltype where it is
    maximum.

    Parameters
    ----------
    adata
    gene_list
    partition_key

    Returns
    -------

    """
    one_v_max = one_v_max_matrix(adata, partition_key=partition_key)
    n_panels = len(gene_list)
    f, grid = _panel_grid(
        hspace=0.25, wspace=0.25, ncols=4, num_panels=n_panels
    )
    for i in range(n_panels):
        gene = gene_list[i]
        ax = f.add_subplot(grid[i])
        ax.set(yscale="log")
        ax.set_yticks = [5 ** (n % 2) * 10 ** (n // 2) for n in range(10)]
        max_celltype = pd.to_numeric(one_v_max[gene]).idxmax()
        distrib_one_v_max(
            adata=adata,
            celltype=max_celltype,
            ax=ax,
            gene_highlight=[gene],
            partition_key=partition_key,
        )


# Not used
def one_v_max_celltypes_sep(
    adata: anndata,
    celltype_list: list,
    gene_highlight=None,
    partition_key: str = "CellType",
):
    """
    Plots the one_v_max distribution across the specified celltypes in
    separate plots

    Parameters
    ----------
    adata
        The corrected expression data.

    celltype_list
        celltypes to plot

    partition_key
        The key in adata.obs corresponding to the annotations to be used.


    Returns
    -------

    """
    n_panels = len(celltype_list)
    f, grid = _panel_grid(
        hspace=0.25, wspace=0.25, ncols=4, num_panels=n_panels
    )
    for i in range(n_panels):
        celltype = celltype_list[i]
        ax = f.add_subplot(grid[i])
        ax.set(yscale="log")
        distrib_one_v_max(
            adata=adata,
            celltype=celltype,
            ax=ax,
            gene_highlight=gene_highlight,
            partition_key=partition_key,
        )


def spec_distrib(
    adata: anndata,
    spe_name: str,
    ax,
    partition_key: str = "CellType",
):
    """
    Plots the distribution of the specified gene specificity

    Parameters
    ----------
    adata


    spe_name
        Specificity to plot

    ax

    partition_key


    Returns
    -------

    """
    dist_to_plot = get_data.get_spe(
        adata=adata, spe_metric=spe_name, partition_key=partition_key
    )
    out = sns.distplot(dist_to_plot, ax=ax)
    return out


def all_spec_distrib(adata: anndata, partition_key: str = "CellType"):
    f, (ax1, ax2, ax3) = plt.subplots(3, 1)
    for ax, spe in [(ax1, "shannon"), (ax2, "tau"), (ax3, "gini")]:
        spec_distrib(
            adata=adata, spe_name=spe, ax=ax, partition_key=partition_key
        )


# TODO
def marked_celltype_UMAP(
    adata: anndata, gene: str, marked_celltypes: Collection[str]
):
    """

    Parameters
    ----------
    adata
    gene
    marked_celltypes

    Returns
    -------

    """


# TODO
def gene_distrib_celltype(adata, gene, celltype):
    """

    Parameters
    ----------
    adata
    gene
    celltype

    Returns
    -------

    """


# if __name__ == "__main__" :
#     a = ['SERPINB4','SERPINB3','LY6D','CLCA2','MT1X','CCND1','AKR1C3',
#     'AQP3','DAPL1','NUPR1','S100A2','KRT5','CSRP2','IGFBP3','CALML3',
#     'TSPAN1','GPX2','FAM3B','EPAS1','TXN','SERPINB13','SUSD4','FAM84A']
#
#     one_v_max_celltypes(barbry,celltype_list=cl,partition_key='CellType')
#
#     ovm = one_v_max_matrix(barbry)
#
#     tp = ovm.loc['Basal']
#     tpt =tp[tp>1]
#
#     one_v_max_genelist(barbry, a, partition_key='CellType')
#
#     all_spec_distrib(barbry, partition_key='CellType')
