"""
Tree Utilities
===============

Phylogenetic tree manipulation and traversal utilities.

Key features:
- Load and parse Newick/Nexus trees
- Prune to species subset
- Node labeling (tax IDs for tips, auto-label or depth-based for internals)
- Polytomy detection
- Root-to-tip path traversal
- MRCA finding
"""

import logging
from pathlib import Path
from typing import List, Optional

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

logger = logging.getLogger(__name__)

from src.asr.tree_parser import build_node_mapping, get_tip_labels


def load_tree(tree_file: Path, format: str = "newick") -> Tree:
    """
    Load phylogenetic tree from file.

    Args:
        tree_file: Path to tree file (Path object or string)
        format: Tree format (newick or nexus, default: newick)

    Returns:
        BioPython Tree object

    Raises:
        FileNotFoundError: If tree file doesn't exist
        ValueError: If tree format is invalid
    """
    # Convert to Path if input is string
    if isinstance(tree_file, str):
        tree_file = Path(tree_file)

    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    try:
        tree = Phylo.read(tree_file, format)
        logger.info(f"Loaded tree from {tree_file}: {tree.count_terminals()} tips")
        return tree
    except Exception as e:
        raise ValueError(f"Failed to load tree from {tree_file}: {e}")


def prune_tree(tree: Tree, species_to_keep: List[str]) -> Tree:
    """
    Prune tree to keep only specified species.

    Args:
        tree: BioPython Tree object
        species_to_keep: List of species names to retain

    Returns:
        Pruned Tree object (new copy)

    Raises:
        ValueError: If no species match or all species would be removed
    """
    import copy

    # Create a deep copy to avoid modifying original
    pruned_tree = copy.deepcopy(tree)

    # Get all terminal names
    all_terminals = {term.name for term in pruned_tree.get_terminals()}
    species_set = set(species_to_keep)

    # Check overlap
    matching_species = all_terminals & species_set
    if not matching_species:
        raise ValueError(
            f"No species in tree match species_to_keep. "
            f"Tree has {len(all_terminals)} species, requested {len(species_set)}"
        )

    # Species to remove
    to_remove = all_terminals - species_set

    if not to_remove:
        logger.debug("No pruning needed, all species already in tree")
        return pruned_tree

    # Prune terminals
    for species in to_remove:
        try:
            pruned_tree.prune(species)
        except Exception as e:
            logger.warning(f"Failed to prune {species}: {e}")

    logger.info(
        f"Pruned tree from {len(all_terminals)} to {pruned_tree.count_terminals()} species "
        f"({len(matching_species)} kept, {len(to_remove)} removed)"
    )

    return pruned_tree


def get_mrca(tree: Tree, terminals: List[str]) -> Optional[Clade]:
    """
    Find most recent common ancestor (MRCA) of specified terminals.

    Args:
        tree: BioPython Tree object
        terminals: List of terminal names

    Returns:
        MRCA Clade or None if not found
    """
    if not terminals:
        return None

    if len(terminals) == 1:
        # Single terminal, return itself
        for term in tree.get_terminals():
            if term.name == terminals[0]:
                return term
        return None

    # Find MRCA using BioPython
    terminal_clades = []
    for term_name in terminals:
        for term in tree.get_terminals():
            if term.name == term_name:
                terminal_clades.append(term)
                break

    if len(terminal_clades) != len(terminals):
        logger.warning(f"Could not find all terminals in tree: {terminals}")
        return None

    mrca = tree.common_ancestor(terminal_clades)
    return mrca


def build_tree_node_mapping(tree_file: Path, rst_file: Path = None):
    """
    Expose ASR tree parser build_node_mapping via phylo module.

    Args:
        tree_file: Path to Newick tree file
        rst_file: Optional path to RST file (preferred source for PAML node IDs)

    Returns:
        Tuple of (ordered_nodes, id_mapping)
    """
    return build_node_mapping(tree_file=tree_file, rst_file=rst_file)


def extract_tip_labels(root_node) -> List[str]:
    """Expose get_tip_labels via phylo module for consistent access."""
    return get_tip_labels(root_node)
