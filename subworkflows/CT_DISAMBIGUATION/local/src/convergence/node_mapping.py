"""
Node Mapping for Convergence Diagnostics
=========================================

Builds phylogenetic node mappings and path diagnostics for contrast-based
convergence analysis. Maps node IDs to states and traces paths from MRCA to tips.

Key Functions
-------------
**build_convergence_node_mapping**: Construct node→state map and path diagnostics

Path Building Logic
-------------------
For each contrast, the function:
    1. Maps the MRCA node to its modal ancestral state
    2. Maps all tip taxa (top/bottom) to their side labels
    3. Builds phylogenetic paths from MRCA to each tip
    4. Annotates internal nodes along paths with ASR states

Output Structure
----------------
Returns a dictionary with:
    - **node_map**: {node_id: state_or_label} for all nodes on paths
    - **mrca_node_id**: Node ID of the contrast MRCA
    - **top_paths**: List of paths (node ID sequences) from MRCA to top tips
    - **bottom_paths**: List of paths from MRCA to bottom tips

Usage Example
-------------
::

    from src.convergence.node_mapping import build_convergence_node_mapping

    mapping = build_convergence_node_mapping(
        tree, contrast, node_posteriors, paml_site=85
    )
    print(mapping['node_map'])  # {node_id: 'A', ...}
    print(mapping['top_paths'])  # [[mrca_id, internal_id, tip_id], ...]

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-07
"""

from typing import Any, Dict, List, Optional, Set
import logging

logger = logging.getLogger(__name__)


def build_convergence_node_mapping(
    tree,
    contrast,
    node_posteriors: Optional[Dict] = None,
    paml_site: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Build node → state mapping for a contrast, showing the phylogenetic context.

    Args:
        tree: Bio.Phylo tree object
        contrast: ContrastDefinition with mrca node_id and taxa
        node_posteriors: Optional ASR posterior probabilities
        paml_site: Optional 1-based PAML site index

    Returns:
        Dictionary with node_map (node_id → state/label) and path diagnostics
    """
    from src.phylo.tree_utils import find_node_by_id

    def get_modal_state(node_id: Optional[int]) -> Optional[str]:
        """Extract modal AA from ASR posteriors"""
        if not node_posteriors or node_id is None or paml_site is None:
            return None
        node_dict = node_posteriors.get(node_id)
        if not node_dict:
            return None
        site_post = node_dict.get(paml_site, {})
        if not site_post:
            return None
        try:
            return max(site_post.items(), key=lambda x: x[1])[0]
        except Exception:
            return None

    node_map: Dict[int, str] = {}

    # Map MRCA
    mrca_node = find_node_by_id(tree, contrast.node_id) if contrast.node_id else None
    if mrca_node and mrca_node.node_id is not None:
        mrca_state = get_modal_state(mrca_node.node_id)
        node_map[mrca_node.node_id] = mrca_state or f"MRCA_{contrast.pair_id}"
    else:
        logger.debug("No MRCA node found for contrast %s", getattr(contrast, "pair_id", "?"))

    # Map all taxa tips
    for taxid in contrast.all_taxa:
        node = find_node_by_id(tree, taxid)
        if node and node.node_id is not None:
            side = "top" if taxid in contrast.top_taxa else "bottom"
            node_map[node.node_id] = f"tip_{side}_{taxid}"

    # Build path from MRCA to each tip
    def build_path_to_node(target_id: int) -> List[int]:
        """Build path from MRCA to target node"""
        if mrca_node is None or target_id is None:
            return []

        target_node = find_node_by_id(tree, target_id)
        if not target_node:
            return []

        # Build path from target to root
        path_to_root = []
        current = target_node
        while current:
            if current.node_id is not None:
                path_to_root.append(current.node_id)
            if current == mrca_node:
                break
            current = tree.get_path(current)[-2] if len(tree.get_path(current)) > 1 else None

        return list(reversed(path_to_root))

    # Build paths for each side
    top_paths = [build_path_to_node(taxid) for taxid in contrast.top_taxa]
    bottom_paths = [build_path_to_node(taxid) for taxid in contrast.bottom_taxa]

    # Collect all internal nodes on paths
    all_path_nodes: Set[int] = set()
    for path in top_paths + bottom_paths:
        all_path_nodes.update(path)

    # Add internal node states to mapping
    for node_id in all_path_nodes:
        if node_id not in node_map:
            state = get_modal_state(node_id)
            node_map[node_id] = state or f"internal_{node_id}"

    return {
        "node_map": node_map,
        "mrca_node_id": contrast.node_id,
        "top_paths": top_paths,
        "bottom_paths": bottom_paths,
        "path_nodes": sorted(all_path_nodes),
    }
