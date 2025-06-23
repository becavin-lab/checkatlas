# Cell-Cycle Conservation Metric Analysis

## Description 

The Cell-Cycle Conservation metric evaluates the ability of an integration method to preserve cell cycle effects after batch effect correction.
This metric is crucial because the cell cycle represents an important source of biological variation that must be preserved during single-cell data integration.
The metric assesses how well cell cycle effects can be captured before and after integration.
It specifically focuses on the S and G2/M phases of the cell cycle, which are the most easily identifiable through gene expression.
It calculates scores for S and G2/M cell cycle phases using [Scanpy's score_cell_cycle function](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html) and measures the difference between variance contribution before and after integration.
A good integration method should eliminate batch effects while preserving authentic biological signals related to the cell cycle. 
This metric is part of the "label-free" metrics that evaluate biological conservation beyond cell type annotations.


## Formulas 

The Cell-Cycle Conservation formula is based on calculating cell cycle scores and measuring their preservation :

*Cell cycle scores calculation* : 
- S score : $Score_S = \text{mean}(\text{expression}{S_genes}) - \text{mean}(\text{expression}{control_genes})$
- G2/M score : $Score_{G2M} = \text{mean}(\text{expression}{G2M_genes}) - \text{mean}(\text{expression}{control_genes})$

*Variance contribution before integration* :
- $Var_{before_S} = \frac{\text{var}(Score_{S_before})}{\text{var}_{total_before}}$
- $Var_{before_G2M} = \frac{\text{var}(Score_{G2M_before})}{\text{var}_{total_before}}$

*Variance contribution after integration* :
- $Var_{after_S} = \frac{\text{var}(Score_{S_after})}{\text{var}_{total_after}}$
- $Var_{after_G2M} = \frac{\text{var}(Score_{G2M_after})}{\text{var}_{total_after}}$

*Final conservation score* :
- $CC_{conservation} = 1 - \frac{|Var_{before} - Var_{after}|}{Var_{before}}$
- Where $Var_{before}$ and $Var_{after}$ represent the combined variance contributions of S and G2/M phases

A score close to 1 indicates excellent preservation of cell cycle signal, while a score close to 0 indicates significant loss of this biological signal.

## Sources 

[Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022).](https://doi.org/10.1038/s41592-021-01336-8)

[Reut Danino, Iftach Nachman, Roded Sharan, Batch correction of single-cell sequencing data via an autoencoder architecture, Bioinformatics Advances, Volume 4, Issue 1, 2024](https://doi.org/10.1093/bioadv/vbad186)

[Open problemns](https://openproblems.bio/results/batch_integration_embed/)

## Code 

[Sacnpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html)

