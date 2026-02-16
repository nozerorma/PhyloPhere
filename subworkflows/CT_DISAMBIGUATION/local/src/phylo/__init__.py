"""
Phylo Module
============

Phylogenetic tree manipulation and utilities.

Components:
    - tree_utils: Tree parsing, pruning, node labeling, path traversal
"""

from .tree_utils import (
    load_tree,
    prune_tree,
    label_nodes,
    get_mrca,
    get_root_to_tip_paths,
    build_node_path_from_mapping,
    build_tip_node_lookup,
    build_tree_node_mapping,
    extract_tip_labels,
)
from .species_mapping import (
    read_taxid_mapping,
    validate_taxids_in_tree,
    match_tree_alignment_by_taxid,
)

__all__ = [
    "load_tree",
    "prune_tree",
    "label_nodes",
    "get_mrca",
    "get_root_to_tip_paths",
    "build_node_path_from_mapping",
    "build_tip_node_lookup",
    "build_tree_node_mapping",
    "extract_tip_labels",
    "read_taxid_mapping",
    "validate_taxids_in_tree",
    "match_tree_alignment_by_taxid",
]
