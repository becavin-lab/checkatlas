import numpy as np
from scipy.stats import spearmanr


def run(adata_before, adata_after, cell_cycle_key='phase'):
    """
    Calculate Cell Cycle Conservation score after batch correction or integration.
    
    This metric measures how well cell cycle phase information is preserved
    after data integration or batch correction. It compares cell cycle scores
    or phase assignments before and after the integration process.
    
    `Cell Cycle Conservation readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/cell_cycle/>`__

    :param adata_before: AnnData object
        Annotated data matrix before integration/batch correction,
        should contain cell cycle information in .obs[cell_cycle_key]
    :param adata_after: AnnData object
        Annotated data matrix after integration/batch correction
    :param cell_cycle_key: str, default='phase'
        Key in .obs containing cell cycle phase or scores
    :return: float
        Cell cycle conservation score (Spearman correlation).
        Range: [-1, 1], where 1 indicates perfect conservation.
    """
    # Extract cell cycle information
    if cell_cycle_key not in adata_before.obs.columns:
        raise ValueError(f"Cell cycle key '{cell_cycle_key}' not found in adata_before.obs")
    
    if cell_cycle_key not in adata_after.obs.columns:
        raise ValueError(f"Cell cycle key '{cell_cycle_key}' not found in adata_after.obs")
    
    cc_before = adata_before.obs[cell_cycle_key]
    cc_after = adata_after.obs[cell_cycle_key]
    
    # Convert categorical to numeric if necessary
    if cc_before.dtype == 'object' or hasattr(cc_before, 'cat'):
        # Map phases to numeric values
        phase_map = {'G1': 0, 'S': 1, 'G2M': 2, 'G2/M': 2}
        cc_before_numeric = cc_before.map(phase_map).fillna(0).astype(float)
    else:
        cc_before_numeric = cc_before.values
    
    if cc_after.dtype == 'object' or hasattr(cc_after, 'cat'):
        phase_map = {'G1': 0, 'S': 1, 'G2M': 2, 'G2/M': 2}
        cc_after_numeric = cc_after.map(phase_map).fillna(0).astype(float)
    else:
        cc_after_numeric = cc_after.values
    
    # Calculate Spearman correlation
    correlation, _ = spearmanr(cc_before_numeric, cc_after_numeric)
    
    return correlation
