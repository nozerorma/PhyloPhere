# __init__.py — Package exports of the data models and the Ensembl gene loader.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/data/

"""
Public exports for data models (BiochemResults, CAASPosition, ConvergenceResult) and `load_ensembl_genes`.

Imported by: the submodules are imported directly (`src.data.<module>`) by src/convergence/disambiguate_single.py,
src/core/observed.py, src/utils/gene_wrapper.py, explain_positions.py and observed_b0_main.py
"""

from .models import BiochemResults, CAASPosition, ConvergenceResult
from .loaders import load_ensembl_genes

__all__ = [
    "BiochemResults",
    "CAASPosition",
    "ConvergenceResult",
    "load_ensembl_genes",
]
