"""
Created on Thu Jun  4 16:38:47 2020

@author: antoinecollin
"""

import datetime
import os
from pathlib import Path

import anndata
import compute
import matplotlib.pyplot as plt
import pandas as pd
import plot
from get_data import get_average_celltype_counts, get_spe


def specificity_summary(adata, marker_genes, partition_key: str = "CellType"):
    """
    Returns a summary of the specificity metrics for the dataset.

    Parameters
    ----------
    adata
        The corrected expression matrix

    marker_genes
        A dictionary containing the celltypes and their marker genes to analyse
        in the following format {celltype1 : [gene1,gene2,...],celltype2 : ...}

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    Returns
    -------
    summary
        The summary as a pandas DataFrame with the following columns :
        gene : name of the gene
        expected celltype : supposedly marked celltype by this gene
        most expressed celltype : celltype in which the gene is most expressed
        one_v_max,shannon,tau,gini : specificity values for this gene

    """
    shannon = get_spe(adata, spe_metric="shannon", partition_key=partition_key)
    tau = get_spe(adata, spe_metric="tau", partition_key=partition_key)
    gini = get_spe(adata, spe_metric="gini", partition_key=partition_key)
    one_v_max = compute.one_v_max_matrix(adata, partition_key=partition_key)
    summary = pd.DataFrame(
        columns=[
            "gene",
            "expected celltype",
            "most expressed celltype",
            "one_v_max",
            "shannon",
            "tau",
            "gini",
        ]
    )

    for celltype, genes in marker_genes.items():
        for gene in genes:
            if (
                gene in shannon.keys()
            ):  # Check why some of them return an error
                summary.loc[gene, :] = [
                    gene,
                    celltype,
                    pd.to_numeric(one_v_max[gene]).idxmax(),
                    max(one_v_max[gene]),
                    shannon[gene],
                    tau[gene],
                    gini[gene],
                ]
    return summary


# TODO
def specificity_quality_control(
    adata: anndata,
    marker_genes: dict,
    partition_key: str,
    project_dir: str,
):
    """
    Performs an analysis of the identified marker genes with regard to
    specificity and save the files of the analysis.

    Parameters
    ----------
    adata
        Annotated expression matrix

    marker_genes
        Marker genes dict to analyze

    partition_key
        The key in adata.obs corresponding to the annotations to be used.

    project_dir
        Path of the project where the plots should be saved


    Returns
    -------
    For every celltypes, plot and save 2 figures :
    -The distribution of the marker genes supported by their shannon, gini,
    au specificity
    -The one_v_max specificity repartition of the celltype for which the marker
     genes are most specific.

    Also computes and save a summary of the analysis

    """
    celltypes = list(marker_genes.keys())
    celltypes_checks = list(
        get_average_celltype_counts(adata, partition_key=partition_key).index
    )
    if sorted(celltypes) != sorted(celltypes_checks):
        print("Celltypes from anndata and marker genes file don't match")
        celltypes = list(set(celltypes) & set(celltypes_checks))
        marker_genes = {keys: marker_genes[keys] for keys in celltypes}
    now = datetime.datetime.now()
    project_dir_path = Path(project_dir)
    if not project_dir_path.joinpath(Path("SpecAnalysis")).is_dir():
        os.mkdir(project_dir_path.joinpath("SpecAnalysis"))
    os.mkdir(
        project_dir_path.joinpath(
            f'SpecAnalysis/Analysis_{now.strftime("%y%m%d_%H%M")}'
        )
    )
    analysis_path = project_dir_path.joinpath(
        f'SpecAnalysis/Analysis_{now.strftime("%y%m%d_%H%M")}'
    )
    for celltype in celltypes:
        cell_path = analysis_path.joinpath(f"{celltype}")
        os.mkdir(cell_path)
        gene_list = list(set(marker_genes[celltype]) & set(adata.var.index))
        # in case some gene names are different/don't exist between
        # the marker gene file and adata file
        plot.marker_genes_distribution(
            adata=adata,
            gene_list=gene_list,
            celltype=celltype,
            partition_key=partition_key,
        )
        plt.savefig(cell_path.joinpath(f"gene_distribution_{celltype}"))
        plt.close()
        plot.one_v_max_genelist(
            adata=adata, gene_list=gene_list, partition_key=partition_key
        )
        plt.savefig(cell_path.joinpath(f"one_v_max_for_{celltype}_genes"))
        plt.close()
    summary = specificity_summary(
        adata=adata, marker_genes=marker_genes, partition_key=partition_key
    )
    summary.to_csv(
        path_or_buf=analysis_path.joinpath("spec_summary.tsv"), sep="\t"
    )
    plot.all_spec_distrib(adata=adata, partition_key=partition_key)
    plt.savefig(analysis_path.joinpath("spec_distrib.png"))
    plt.close()


# if __name__ == '__main__':
#     specificity_quality_control(adata=barbry,
#                                 marker_genes=marker_dict,
#                                 partition_key='CellType',
#                                 project_dir = r'
#                                 C:\Users\ipmc\PycharmProjects\test&tuto\Kobaye')

# Saves in an ad hoc 'result' directory :
# - gene_distribs, one for each celltype
# - one_v_max distrib : one_v_max distrib in the celltypes. Markers of
# different celltypes marked with different colors ? --> messy.
# For each gene, just the stripplot
# one_v_max in the most expressed celltype. (with legend pex 'FOXJ1,
# max in 'Basal')
# - an excel file with gene | expected celltype | most expressed celltype |
# one_v_max | Shannon | Tau | Gini


# adata = barbry
# genes = dict({celltype : [Gene1,Gene2,...]}. On gerera les
# autres types de fichier
# plus tard


# TODO
def specificity_gene_selection(adata: anndata, many, other, arguments):
    """

    Parameters
    ----------
    adata
    many
    other
    arguments

    Returns
    -------

    """
    return "do_nothing"
