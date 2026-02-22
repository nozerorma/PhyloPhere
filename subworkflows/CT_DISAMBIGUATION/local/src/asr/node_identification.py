"""
Node Identification for Convergence Analysis
=============================================

Helper functions to identify key phylogenetic nodes from tree structure
and species groupings for convergence analysis with dynamic pair support.

**REDESIGNED FOR DYNAMIC PAIRS (2025-12-05)**:
Supports variable numbers of focal pairs (1-N) rather than hard-coded focal_1/focal_2.

Key Nodes:
    X1 (root): Root node (ancestral state at tree root)
    X2 (mrca_contrast): MRCA of ALL focal nodes (encompasses all pairs)
    focal_nodes[0] through focal_nodes[N-1]: MRCA for each pair dynamically

Usage:
    >>> from asr.node_identification import identify_convergence_nodes
    >>> from asr.tree_parser import parse_newick
    >>>
    >>> # Parse tree
    >>> tree_str = Path("tree.nwk").read_text()
    >>> root = parse_newick(tree_str)
    >>>
    >>> # Identify nodes for variable pairs
    >>> node_mapping = identify_convergence_nodes(
    ...     root,
    ...     focal_pairs=[
    ...         ['taxid1', 'taxid2'],  # pair 1
    ...         ['taxid3', 'taxid4'],  # pair 2
    ...         ['taxid5', 'taxid6'],  # pair 3
    ...     ]
    ... )
    >>> # Returns: {'root': 116, 'mrca_contrast': 98, 'focal_nodes': [55, 78, 42], ...}

    Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-03
"""

from typing import List, Dict, Optional, Any
from pathlib import Path
import logging

from .tree_parser import (
    TreeNode,
    parse_newick,
    get_mrca,
    find_node_by_name,
    get_tip_labels,
)

logger = logging.getLogger(__name__)


def get_descendants(root: TreeNode) -> List[TreeNode]:
    """Return all descendant nodes under ``root`` (excluding ``root``)."""
    descendants: List[TreeNode] = []
    for child in root.children:
        descendants.append(child)
        descendants.extend(get_descendants(child))
    return descendants


def identify_convergence_nodes(
    root: TreeNode,
    focal_pairs: List[List[str]],
    contrast_species: Optional[List[str]] = None,
    require_all: bool = False,
) -> Dict[str, Any]:
    """
    Identify key phylogenetic nodes for convergence analysis with dynamic pair support.

    **REDESIGNED FOR DYNAMIC PAIRS**: This function now accepts a variable number of
    focal pairs (1-N) and computes focal nodes dynamically for each pair. The MRCA
    of the contrast is computed as the MRCA of ALL focal nodes.

    Args:
        root: Root TreeNode of parsed phylogeny
        focal_pairs: List of species ID lists, one per pair. Each list contains the
            species IDs for that pair (union of TOP and BOTTOM species).
            Example: [['human', 'chimp'], ['mouse', 'rat'], ['dog', 'cat']]
        contrast_species: Optional list of all species across all pairs (defaults to
            flattened focal_pairs)
        require_all: If True, raise error if any species not found in tree

    Returns:
        Dictionary with:
            - 'root': Tree root node ID
            - 'mrca_contrast': MRCA of all focal nodes (encompasses all pairs)
            - 'focal_nodes': List of focal node IDs, one per pair [focal_i_id ...]
            - 'mrca_all': MRCA of all tip species (optional)

    Raises:
        ValueError: If required species not found in tree or focal_pairs is empty

    Example:
        >>> root = parse_newick(tree_string)
        >>> nodes = identify_convergence_nodes(
        ...     root,
        ...     focal_pairs=[['human', 'chimp'], ['mouse', 'rat'], ['dog', 'cat']]
        ... )
        >>> print(nodes)
        {'root': 116, 'mrca_contrast': 98, 'focal_nodes': [78, 55, 42], 'mrca_all': 100}
    """
    if not focal_pairs:
        raise ValueError("focal_pairs cannot be empty")

    for idx, group in enumerate(focal_pairs, 1):
        if not group:
            raise ValueError(f"focal_pairs[{idx-1}] (pair {idx}) cannot be empty")

    # Node IDs should already be set (typically from PAML RST parsing)
    if root.node_id is None:
        logger.warning("Root node has no node_id set - node IDs may not be available")

    # Get all tip labels in tree (format typically lineage_taxid); also capture taxid suffixes
    all_tips = get_tip_labels(root)
    tip_set = set(all_tips)
    tip_taxid_set = set(label.split("_")[-1] for label in all_tips if label)

    logger.debug(f"Tree has {len(all_tips)} tips")
    logger.debug(f"Received {len(focal_pairs)} focal pairs")

    # Check which species are present (by full label or taxid suffix)
    def _present(species_list):
        present = []
        for sp in species_list:
            s = str(sp)
            if s in tip_set or s in tip_taxid_set:
                present.append(s)
        return present

    # Validate each focal pair
    focal_pairs_present = []
    for idx, group in enumerate(focal_pairs, 1):
        present = _present(group)
        if not present:
            msg = f"No species from focal pair {idx} found in tree (looked for {len(group)} species)"
            if require_all:
                raise ValueError(msg)
            logger.warning(msg)
        else:
            logger.info(
                f"Found {len(present)}/{len(group)} species for focal pair {idx}"
            )
        focal_pairs_present.append(present)

    # Determine contrast species (all species across all pairs)
    if contrast_species is None:
        contrast_species = [sp for group in focal_pairs for sp in group]

    contrast_present = _present(contrast_species)

    if not contrast_present:
        msg = f"No contrast species found in tree (looked for {len(contrast_species)} species)"
        if require_all:
            raise ValueError(msg)
        logger.warning(msg)

    logger.info(
        f"Found {len(contrast_present)}/{len(contrast_species)} total contrast species"
    )

    # Initialize node mapping
    node_mapping: Dict[str, Any] = {}

    # X1: Root node
    assert root.node_id is not None
    node_mapping["root"] = root.node_id
    logger.debug(f"Root (X1): node {root.node_id}")

    # MRCA of all tips
    if len(all_tips) > 1:
        mrca_all = get_mrca(root, all_tips)
        if mrca_all and mrca_all.node_id is not None:
            node_mapping["mrca_all"] = mrca_all.node_id
            logger.debug(f"MRCA all tips: node {mrca_all.node_id}")

    # Compute focal nodes for each pair dynamically
    focal_node_ids = []
    for idx, group_present in enumerate(focal_pairs_present, 1):
        if len(group_present) >= 2:
            mrca_focal = get_mrca(root, group_present)
            if mrca_focal and mrca_focal.node_id is not None:
                focal_node_ids.append(mrca_focal.node_id)
                logger.debug(f"MRCA focal pair {idx}: node {mrca_focal.node_id}")
            else:
                logger.warning(f"Could not find MRCA of focal pair {idx}")
                focal_node_ids.append(None)
        elif len(group_present) == 1:
            # Single species - use the tip node itself or its parent
            node = find_node_by_name(root, group_present[0])
            if node and node.node_id is not None:
                focal_node_ids.append(node.node_id)
                logger.debug(f"Focal pair {idx} (single species): node {node.node_id}")
            else:
                logger.warning(f"Could not find node for focal pair {idx}")
                focal_node_ids.append(None)
        else:
            logger.warning(f"Focal pair {idx} has no present species")
            focal_node_ids.append(None)

    node_mapping["focal_nodes"] = focal_node_ids

    # X2: MRCA of contrast (MRCA of ALL focal nodes, not just focal_1 and focal_2)
    # This ensures mrca_contrast encompasses all pairs, not just the first two
    valid_focal_ids = [fid for fid in focal_node_ids if fid is not None]

    if len(valid_focal_ids) >= 2:
        # Find nodes corresponding to focal IDs
        focal_nodes_objs = []
        for node in [root] + list(get_descendants(root)):
            if node.node_id in valid_focal_ids:
                focal_nodes_objs.append(node)

        if len(focal_nodes_objs) >= 2:
            # Get all tip labels under these focal nodes
            all_focal_tips = []
            for fnode in focal_nodes_objs:
                all_focal_tips.extend(get_tip_labels(fnode))

            # Compute MRCA of all these tips
            mrca_contrast = get_mrca(root, all_focal_tips)
            if mrca_contrast and mrca_contrast.node_id is not None:
                node_mapping["mrca_contrast"] = mrca_contrast.node_id
                logger.debug(
                    f"MRCA contrast (all {len(valid_focal_ids)} focal nodes): node {mrca_contrast.node_id}"
                )
            else:
                logger.warning("Could not find MRCA of contrast species")
        else:
            logger.warning("Could not find focal node objects for MRCA computation")
    elif len(valid_focal_ids) == 1:
        # Single focal node - use it as mrca_contrast
        node_mapping["mrca_contrast"] = valid_focal_ids[0]
        logger.debug(f"MRCA contrast (single focal node): node {valid_focal_ids[0]}")

    # Validate we have minimum required data
    if not focal_node_ids:
        msg = "No focal nodes could be identified"
        if require_all:
            raise ValueError(msg)
        logger.warning(msg)

    if "mrca_contrast" not in node_mapping and require_all:
        raise ValueError("Could not identify MRCA of contrast")

    return node_mapping


def identify_convergence_nodes_from_file(
    tree_file: Path,
    focal_pairs: List[List[str]],
    contrast_species: Optional[List[str]] = None,
    require_all: bool = False,
) -> Dict[str, Any]:
    """
    Convenience function to identify convergence nodes from tree file with dynamic pairs.

    Args:
        tree_file: Path to Newick tree file
        focal_pairs: List of species ID lists, one per pair
        contrast_species: Optional list of all contrast species
        require_all: If True, raise error if any species not found

    Returns:
        Dictionary with node mapping (root, mrca_contrast, focal_nodes list)

    Example:
        >>> nodes = identify_convergence_nodes_from_file(
        ...     Path("tree.nwk"),
        ...     focal_pairs=[['human', 'chimp'], ['mouse', 'rat']]
        ... )
    """
    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    tree_str = tree_file.read_text().strip()
    root = parse_newick(tree_str)

    # Note: For PAML trees, node IDs should be loaded via build_node_mapping with RST file
    if root.node_id is None:
        logger.warning(
            "Tree nodes have no IDs - use build_node_mapping with RST file for PAML trees"
        )

    return identify_convergence_nodes(root, focal_pairs, contrast_species, require_all)


def validate_node_mapping(
    node_mapping: Dict[str, int], required_keys: Optional[List[str]] = None
) -> bool:
    """
    Validate that node mapping has required keys and valid node IDs.

    Args:
        node_mapping: Dictionary of role names to node IDs
        required_keys: List of required keys (default: ['root', 'mrca_contrast', 'focal_nodes'])

    Returns:
        True if valid, False otherwise
    """
    if required_keys is None:
        required_keys = ["root", "mrca_contrast", "focal_nodes"]

    # Check all required keys present
    missing = [k for k in required_keys if k not in node_mapping]
    if missing:
        logger.error(f"Missing required keys in node mapping: {missing}")
        return False

    # Check all node IDs are non-negative integers
    for key, node_id in node_mapping.items():
        if not isinstance(node_id, int):
            logger.error(f"Node ID for '{key}' is not an integer: {type(node_id)}")
            return False
        if node_id < 0:
            logger.error(f"Node ID for '{key}' is negative: {node_id}")
            return False

    return True
