"""
ASR Module
==========

Ancestral State Reconstruction using PAML.

Components:
    - reconstruct: Core ASR logic and PAML wrapper (ASRConfig, ASRReconstructor)
    - posterior: PAML posterior probability parsing and handling
    - tree_parser: Newick tree parsing and node mapping utilities
"""

from .reconstruct import ASRConfig, ASRReconstructor
from .posterior import (
    parse_paml_rst,
    parse_paml_rst_node_level,
)
from .tree_parser import (
    TreeNode,
    parse_newick,
    build_node_mapping,
    get_node_order,
    get_tip_labels,
    find_node_by_name,
    get_mrca,
)
from .node_identification import (
    identify_convergence_nodes,
    identify_convergence_nodes_from_file,
    validate_node_mapping,
)

__all__ = [
    "ASRConfig",
    "ASRReconstructor",
    "parse_paml_rst",
    "parse_paml_rst_node_level",
    "parse_newick",
    "build_node_mapping",
    "get_node_order",
    "get_tip_labels",
    "find_node_by_name",
    "get_mrca",
    "identify_convergence_nodes",
    "identify_convergence_nodes_from_file",
    "validate_node_mapping",
]
