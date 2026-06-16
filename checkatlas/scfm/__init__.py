"""scFM QC layer 2 (diagnostics) and layer 3 (composite, reporting).

This package is the *interpretation* layer. It takes the long-format
metric table produced by ``checkatlas.metrics.scfm.run.cal_scfm`` and
the existing ``cal_annot/cluster/dimred`` and produces:

  * Nine ``ProblemVerdict`` objects (one per scFM failure mode).
  * Three composite scores (FMF, BF, PR).
  * MultiQC-friendly TSV outputs.

It is intentionally separate from ``checkatlas.metrics.scfm`` so
Layer 1 (metric collection) can be re-used without the diagnostic
overhead.
"""

from . import config, diagnostics, grades, report, rules

__all__ = ["config", "diagnostics", "grades", "report", "rules"]
