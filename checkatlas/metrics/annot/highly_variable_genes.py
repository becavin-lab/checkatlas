import numpy as np
import scanpy as sc


def run(adata_before, adata_after, n_genes=2000, method='seurat'):
    """
    Calculate Highly Variable Genes conservation after batch correction.
    
    This metric measures the overlap of highly variable genes (HVGs) before and
    after batch correction or integration. High overlap indicates that biological
    variation is preserved.
    
    `Highly Variable Genes readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/hvg/>`__

    :param adata_before: AnnData object
        Annotated data matrix before integration/batch correction
    :param adata_after: AnnData object
        Annotated data matrix after integration/batch correction
    :param n_genes: int, default=2000
        Number of highly variable genes to select
    :param method: str, default='seurat'
        Method for identifying HVGs ('seurat', 'cell_ranger', 'seurat_v3')
    :return: float
        Jaccard index of HVG overlap. Range: [0, 1], where 1 indicates perfect overlap.
        
    """
    # Get HVG information from adata_before
    if 'highly_variable' not in adata_before.var.columns:
        # Compute HVGs if not already done
        import scanpy as sc
        sc.pp.highly_variable_genes(adata_before, n_top_genes=n_genes, flavor=method)
    
    hvgs_before = set(adata_before.var_names[adata_before.var['highly_variable']])
    
    # Get HVG information from adata_after
    if 'highly_variable' not in adata_after.var.columns:
        import scanpy as sc
        sc.pp.highly_variable_genes(adata_after, n_top_genes=n_genes, flavor=method)
    
    hvgs_after = set(adata_after.var_names[adata_after.var['highly_variable']])
    
    # Calculate Jaccard index
    intersection = len(hvgs_before & hvgs_after)
    union = len(hvgs_before | hvgs_after)
    
    if union == 0:
        return 0.0
    
    jaccard_index = intersection / union
    
    return jaccard_index
