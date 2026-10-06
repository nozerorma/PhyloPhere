# __init__.py — Package exports of the convergence classification and disambiguation functions.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/convergence/

"""
Public exports for convergence classification and disambiguation: the node-state helpers of convergence.py and the
per-gene and per-position entry points of disambiguate_single.py.

Imported by: the submodules are imported directly (`src.convergence.<module>`) by src/core/driver.py,
src/core/observed.py and src/utils/gene_wrapper.py
"""

from .convergence import (
    NodeStates,
    extract_node_states_from_node_level,
    build_alignment_lookup,
    collect_tip_residues,
    format_amino_display,
)

from .disambiguate_single import (
    analyze_caas_position_disambiguation,
    analyze_gene_disambiguation,
)

__all__ = [
    "NodeStates",
    "analyze_caas_position_disambiguation",
    "analyze_gene_disambiguation",
    "build_alignment_lookup",
    "collect_tip_residues",
    "extract_node_states_from_node_level",
    "format_amino_display",
]
