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
from typing import List, Dict, Optional

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

logger = logging.getLogger(__name__)

from src.asr.tree_parser import build_node_mapping, get_tip_labels


def has_polytomies(tree: Tree) -> bool:
    """Return True when any internal node has more than two children."""
    for clade in tree.get_nonterminals():
        if len(clade.clades) > 2:
            return True
    return False


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


def label_nodes(
    tree: Tree,
    tax_mapping: Optional[Dict[str, str]] = None,
    use_depth: bool = False,
) -> Dict[Clade, str]:
    """
    Label internal nodes in the tree.

    Strategy (from config):
    - Tips: Always use tax_id if available, else species name
    - Internal nodes:
      - If polytomies exist: auto-label (N1, N2, ...)
      - Else: use depth-based naming (depth_X.XXXX)
      - Override with use_depth=True to force depth-based

    Args:
        tree: BioPython Tree object
        tax_mapping: Dict mapping species -> tax_id (e.g., "Homo_sapiens" -> "txid9606")
        use_depth: Force depth-based naming even without polytomies

    Returns:
        Dict mapping Clade -> node_id (str)
    """
    node_labels = {}

    # Determine naming strategy for internal nodes
    has_poly = has_polytomies(tree)
    use_auto_label = has_poly or use_depth is False

    if use_auto_label:
        logger.info("Using auto-label naming for internal nodes (N1, N2, ...)")
    else:
        logger.info("Using depth-based naming for internal nodes (depth_X.XXXX)")

    # Label terminals (tips)
    for terminal in tree.get_terminals():
        if tax_mapping and terminal.name in tax_mapping:
            node_labels[terminal] = tax_mapping[terminal.name]
        else:
            node_labels[terminal] = terminal.name

    # Label internal nodes
    internal_counter = 1
    for clade in tree.get_nonterminals():
        if use_auto_label:
            node_labels[clade] = f"N{internal_counter}"
            internal_counter += 1
        else:
            # Depth-based naming
            depth = tree.distance(clade)
            node_labels[clade] = f"depth_{depth:.4f}"

    logger.info(
        f"Labeled {len(tree.get_terminals())} tips and "
        f"{len(tree.get_nonterminals())} internal nodes"
    )

    return node_labels


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


def get_root_to_tip_paths(
    tree: Tree,
    node_labels: Dict[Clade, str],
) -> Dict[str, List[str]]:
    """
    Get paths from root to all tips.

    Args:
        tree: BioPython Tree object
        node_labels: Dict mapping Clade -> node_id

    Returns:
        Dict mapping tip_node_id -> [node_ids from root to tip]
    """
    paths = {}

    for terminal in tree.get_terminals():
        tip_label = node_labels.get(terminal)
        if not tip_label:
            continue

        # Get path from root to this terminal
        path_clades = tree.get_path(terminal)
        if path_clades is None:
            logger.warning(f"Could not get path for terminal {tip_label}")
            continue
        path_labels = [node_labels.get(c, "UNKNOWN") for c in path_clades]

        paths[tip_label] = path_labels

    return paths


def get_path_between_nodes(
    tree: Tree,
    node_labels: Dict[Clade, str],
    node1_id: str,
    node2_id: str,
) -> Optional[List[str]]:
    """
    Get path between two nodes in the tree.

    Args:
        tree: BioPython Tree object
        node_labels: Dict mapping Clade -> node_id
        node1_id: First node ID
        node2_id: Second node ID

    Returns:
        List of node IDs in path, or None if nodes not found
    """
    # Reverse mapping: node_id -> Clade
    label_to_clade = {label: clade for clade, label in node_labels.items()}

    if node1_id not in label_to_clade or node2_id not in label_to_clade:
        return None

    clade1 = label_to_clade[node1_id]
    clade2 = label_to_clade[node2_id]

    # Get paths from root to each node
    path1 = tree.get_path(clade1)
    path2 = tree.get_path(clade2)

    # Check if paths are valid
    if path1 is None or path2 is None:
        return None

    # Find common ancestor (last common node in paths)
    common_ancestor = None
    for i in range(min(len(path1), len(path2))):
        if path1[i] == path2[i]:
            common_ancestor = path1[i]
        else:
            break

    if not common_ancestor:
        return None

    # Build path: node1 -> common_ancestor -> node2
    # Path from node1 to common ancestor (reverse)
    ca_idx_in_path1 = path1.index(common_ancestor)
    path_node1_to_ca = path1[ca_idx_in_path1:][::-1]  # Reverse

    # Path from common ancestor to node2
    ca_idx_in_path2 = path2.index(common_ancestor)
    path_ca_to_node2 = path2[ca_idx_in_path2 + 1 :]  # Skip CA itself

    # Combine
    full_path = path_node1_to_ca + path_ca_to_node2

    # Convert to labels
    path_labels = [node_labels.get(c, "UNKNOWN") for c in full_path]

    return path_labels


def build_node_path_from_mapping(
    node_mapping: Dict[int, any], root_id: Optional[int], target_id: Optional[int]
) -> List[int]:
    """
    Build path from root to target node (inclusive) using node mapping dict.

    Args:
        node_mapping: Dictionary mapping node IDs to node objects
        root_id: Root node ID
        target_id: Target node ID

    Returns:
        List of node IDs from root to target (inclusive)
        Empty list if path cannot be constructed
    """
    if not node_mapping or root_id is None or target_id is None:
        return []

    target = node_mapping.get(target_id)
    if target is None:
        return []

    path: List[int] = []
    node = target

    while node is not None:
        node_id = getattr(node, "node_id", None)
        if node_id is None:
            break
        path.insert(0, node_id)
        if node_id == root_id:
            return path
        # Get parent node (TreeNode object, not ID)
        parent = getattr(node, "parent", None)
        if parent is None:
            break
        # Get parent's node_id and look it up in mapping
        parent_node_id = getattr(parent, "node_id", None)
        node = node_mapping.get(parent_node_id) if parent_node_id is not None else None

    return []


def build_tip_node_lookup(node_mapping: Dict[int, any]) -> Dict[str, int]:
    """
    Build mapping from tip labels (taxids/species) to node IDs.

    Args:
        node_mapping: Dictionary mapping node IDs to node objects

    Returns:
        Dictionary mapping tip labels to node IDs
    """
    tip_lookup: Dict[str, int] = {}

    if not node_mapping:
        return tip_lookup

    for node_id, node in node_mapping.items():
        # Check both is_tip and is_leaf (different tree implementations)
        is_tip_node = getattr(node, "is_tip", False) or getattr(node, "is_leaf", False)

        # Also check if it's callable (some implementations use methods)
        if not is_tip_node and hasattr(node, "is_leaf") and callable(node.is_leaf):
            is_tip_node = node.is_leaf()

        if not is_tip_node:
            continue

        # Try multiple label attributes
        label = getattr(node, "taxid", None)
        if not label:
            label = getattr(node, "species", None)
        if not label:
            label = getattr(node, "name", None)

        if label:
            tip_lookup[str(label)] = node_id

    return tip_lookup


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
