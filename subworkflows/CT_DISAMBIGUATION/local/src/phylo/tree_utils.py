# tree_utils.py — Load and prune phylogenetic trees; thin access to the ASR tree parser.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/phylo/

"""
Tree utilities: load a tree with Biopython, prune it to a species subset, and expose the node-mapping and tip-label
helpers of src/asr/tree_parser.py through the phylo package.

Imported by: src/phylo/species_mapping.py (`prune_tree`), src/asr/asr_single.py (`load_tree`),
src/core/driver.py (`build_tree_node_mapping`, `extract_tip_labels`)
Inputs: a tree file (Newick by default)
Outputs: Bio.Phylo trees (in memory)
"""

import logging
from pathlib import Path
from typing import List

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

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
    # Accept a plain string path
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

    # Work on a deep copy: the input tree is left unchanged
    pruned_tree = copy.deepcopy(tree)

    # All terminal names
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


def build_tree_node_mapping(tree_file: Path, rst_file: Path = None):
    """
    PAML-aligned node order and id mapping of a tree (src.asr.tree_parser.build_node_mapping).

    Args:
        tree_file: Path to Newick tree file
        rst_file: Optional path to RST file (preferred source for PAML node IDs)

    Returns:
        Tuple of (ordered_nodes, id_mapping)
    """
    return build_node_mapping(tree_file=tree_file, rst_file=rst_file)


def extract_tip_labels(root_node) -> List[str]:
    """Tip labels of the tree below `root_node` (src.asr.tree_parser.get_tip_labels)."""
    return get_tip_labels(root_node)
