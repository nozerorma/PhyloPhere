# __init__.py — Package exports of the amino acid grouping schemes.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/biochem/

"""
Public exports for biochemical grouping: the scheme tables (US, GS1 to GS4) and `get_grouping_scheme`.

Imported by: src/convergence/disambiguate_single.py, src/convergence/path_scores.py
"""

from .grouping import US, GS1, GS2, GS3, GS4, get_grouping_scheme

__all__ = [
    "US",
    "GS1",
    "GS2",
    "GS3",
    "GS4",
    "get_grouping_scheme",
]
