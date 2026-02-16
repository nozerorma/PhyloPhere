"""
Tree parsing utilities for PAML RST file interpretation.

This module provides functions to parse Newick tree files with PAML node labels
from RST output files. PAML assigns specific node IDs (visible in the 'tree with
node labels for Rod Page's TreeView' section) that must be preserved for correct
posterior-to-node mapping.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import logging

logger = logging.getLogger(__name__)


class TreeNode:
    """Represents a node in a phylogenetic tree."""
    
    def __init__(self, name: Optional[str] = None, node_id: Optional[int] = None):
        self.name = name  # Tip label (e.g., 'taxid_speciesname') or None for internal nodes
        self.node_id = node_id  # PAML node ID (from RST labeled tree)
        self.children: List[TreeNode] = []
        self.parent: Optional[TreeNode] = None
        self.branch_length: float = 0.0
        self.paml_label: Optional[str] = None  # Original PAML label if present
    
    def is_leaf(self) -> bool:
        """Check if this is a leaf (tip) node."""
        return len(self.children) == 0
    
    def __repr__(self) -> str:
        if self.name:
            return f"TreeNode(name={self.name}, id={self.node_id})"
        return f"TreeNode(id={self.node_id})"


def parse_newick(tree_string: str) -> TreeNode:
    """
    Parse a Newick format tree string (with or without PAML node labels).
    
    For PAML trees with node IDs, use build_node_mapping() with rst_file
    which extracts node IDs from the connection map.
    
    Args:
        tree_string: Newick format tree string
        
    Returns:
        Root node of the parsed tree
    """
    # Preprocess: remove spaces to simplify parsing
    # PAML trees have spaces like ") 187 ," which we convert to ")187,"
    tree_string = tree_string.replace(' ', '')
    tree_string = tree_string.strip()
    if tree_string.endswith(';'):
        tree_string = tree_string[:-1]
    
    stack: List[TreeNode] = []
    current_node: Optional[TreeNode] = None
    current_label = ""
    current_length = ""
    reading_length = False
    reading_node_label = False  # True after ')' to read label for internal node
    
    for char in tree_string:
        if char == '(':
            # Before starting new internal node, finalize any pending label
            if reading_node_label and current_node and current_label:
                current_node.paml_label = current_label
                current_label = ""
                reading_node_label = False
            
            # Start a new internal node
            new_node = TreeNode()
            if current_node is not None:
                new_node.parent = current_node
                current_node.children.append(new_node)
            stack.append(new_node)
            current_node = new_node
            current_label = ""
            current_length = ""
            reading_length = False
            
        elif char == ')':
            # Before finalizing, check if we need to assign a pending label
            if reading_node_label and current_node and current_label:
                current_node.paml_label = current_label
                current_label = ""
                reading_node_label = False
            
            # Finalize any pending leaf node
            if current_label or current_length:
                child = TreeNode(name=current_label if current_label else None)
                if current_length:
                    try:
                        child.branch_length = float(current_length)
                    except ValueError:
                        logger.warning(f"Invalid branch length: {current_length}")
                if current_node:
                    child.parent = current_node
                    current_node.children.append(child)
            
            # Pop to parent - next text will be label for THIS node
            if stack:
                current_node = stack.pop()
            current_label = ""
            current_length = ""
            reading_length = False
            reading_node_label = True  # Next text is label for the node just closed
            
        elif char == ',':
            # Finalize leaf or assign internal node label
            if reading_node_label:
                # Assign label to the node we just closed (current_node)
                if current_node and current_label:
                    current_node.paml_label = current_label
                reading_node_label = False
                # After assigning label, move back to parent to handle next sibling
                if current_node and current_node.parent:
                    current_node = current_node.parent
            elif current_label or current_length:
                # Create a leaf child
                child = TreeNode(name=current_label if current_label else None)
                if current_length:
                    try:
                        child.branch_length = float(current_length)
                    except ValueError:
                        logger.warning(f"Invalid branch length: {current_length}")
                if current_node:
                    child.parent = current_node
                    current_node.children.append(child)
            
            current_label = ""
            current_length = ""
            reading_length = False
                
        elif char == ':':
            # Assign any accumulated label before starting branch length
            if reading_node_label and current_node and current_label:
                current_node.paml_label = current_label
                current_label = ""
                reading_node_label = False
            reading_length = True
            
        else:
            if reading_length:
                current_length += char
            else:
                current_label += char
    
    # Handle final label (root node label after last ')')
    if reading_node_label and current_node and current_label:
        current_node.paml_label = current_label
    elif current_label or current_length:
        # Trailing content - log debug (common in some tree formats)
        logger.debug(f"Unexpected trailing content after tree: label='{current_label}', length='{current_length}'")
    
    if current_node is None:
        raise ValueError("Failed to parse tree")
    
    return current_node
    
    return current_node


def get_node_order(root: TreeNode) -> List[TreeNode]:
    """
    Get list of nodes in post-order (children before parent).

    Note: PAML posterior tables are *not* in this order; they list nodes
    starting at the root. Use this helper only for structural traversals,
    not for mapping to PAML posterior ordering.

    Args:
        root: Root node of tree

    Returns:
        List of nodes in post-order
    """
    nodes = []
    for child in root.children:
        nodes.extend(get_node_order(child))
    nodes.append(root)
    return nodes


def assign_postorder_ids(root: TreeNode, start_id: int = 0) -> int:
    """
    Assign post-order node IDs to tree nodes.
    
    This is used for non-PAML trees (e.g., input phylogenies before ASR).
    For PAML output trees, use build_node_mapping() with rst_file instead,
    which extracts the actual PAML node IDs.
    
    Args:
        root: Root node of the tree
        start_id: Starting ID number (default 0)
        
    Returns:
        Next available ID after assignment
    """
    current_id = start_id
    for child in root.children:
        current_id = assign_postorder_ids(child, current_id)
    root.node_id = current_id
    current_id += 1
    return current_id


def extract_paml_node_range(rst_file: Path) -> Tuple[int, int]:
    """
    Extract PAML's internal node ID range from RST file.
    
    The RST file contains a line like "Nodes 187 to 371 are ancestral"
    which explicitly declares the internal node numbering.
    
    Args:
        rst_file: Path to PAML RST file
        
    Returns:
        Tuple of (min_node_id, max_node_id)
        
    Raises:
        FileNotFoundError: If RST file doesn't exist
        ValueError: If node range line not found
    """
    if not rst_file.exists():
        raise FileNotFoundError(f"RST file not found: {rst_file}")
    
    with open(rst_file, 'r') as f:
        content = f.read()
    
    # Find the node range declaration
    import re
    match = re.search(r'Nodes\s+(\d+)\s+to\s+(\d+)\s+are\s+ancestral', content)
    if not match:
        raise ValueError(f"Could not find 'Nodes X to Y are ancestral' in {rst_file}")
    
    min_id = int(match.group(1))
    max_id = int(match.group(2))
    
    logger.debug(f"PAML internal node range from RST: {min_id} to {max_id}")
    return min_id, max_id


def extract_paml_tree_structure(rst_file: Path) -> Dict[int, List[int]]:
    """
    Extract tree structure from PAML's connection map in RST file.
    
    The RST contains lines like: "187..188 188..189 189..190 ..."
    which define parent..child relationships.
    
    Args:
        rst_file: Path to PAML RST file
        
    Returns:
        Dict mapping parent node ID -> list of child node IDs
        
    Raises:
        FileNotFoundError: If RST file doesn't exist
        ValueError: If connection map not found
    """
    if not rst_file.exists():
        raise FileNotFoundError(f"RST file not found: {rst_file}")
    
    with open(rst_file, 'r') as f:
        content = f.read()
    
    # Find the connection map line (comes before "tree with node labels")
    # Look for pattern like "187..188 188..189 ..."
    import re
    
    # Find the section between the numeric tree and "tree with node labels"
    tree_marker = "tree with node labels for Rod Page's TreeView"
    tree_idx = content.find(tree_marker)
    if tree_idx == -1:
        raise ValueError(f"Could not find tree marker in {rst_file}")
    
    # Search backwards from tree marker to find the connection map
    search_section = content[max(0, tree_idx - 10000):tree_idx]
    
    # Find all parent..child pairs
    connections = re.findall(r'(\d+)\.\.(\d+)', search_section)
    if not connections:
        raise ValueError(f"Could not find connection map in {rst_file}")
    
    # Build parent -> children mapping
    tree_structure: Dict[int, List[int]] = {}
    for parent_str, child_str in connections:
        parent = int(parent_str)
        child = int(child_str)
        if parent not in tree_structure:
            tree_structure[parent] = []
        tree_structure[parent].append(child)
    
    logger.debug(f"Extracted tree structure with {len(tree_structure)} parent nodes")
    return tree_structure


def extract_rst_labeled_tree(rst_file: Path) -> str:
    """
    Extract the labeled Newick tree from PAML RST file.
    
    The RST file contains a section 'tree with node labels for Rod Page's TreeView'
    which includes PAML's actual node numbering. This is the source of truth.
    
    Args:
        rst_file: Path to PAML RST file
        
    Returns:
        Labeled Newick tree string
        
    Raises:
        FileNotFoundError: If RST file doesn't exist
        ValueError: If labeled tree section not found
    """
    if not rst_file.exists():
        raise FileNotFoundError(f"RST file not found: {rst_file}")
    
    with open(rst_file, 'r') as f:
        content = f.read()
    
    # Find the "tree with node labels for Rod Page's TreeView" section
    marker = "tree with node labels for Rod Page's TreeView"
    start_idx = content.find(marker)
    if start_idx == -1:
        raise ValueError(f"Could not find '{marker}' section in {rst_file}")
    
    # Extract tree string (starts after marker, ends at ';')
    tree_start = start_idx + len(marker)
    tree_section = content[tree_start:]
    
    # Find the semicolon that ends the tree
    semicolon_idx = tree_section.find(';')
    if semicolon_idx == -1:
        raise ValueError(f"Could not find tree terminator (;) in {rst_file}")
    
    tree_string = tree_section[:semicolon_idx + 1].strip()
    
    logger.debug(f"Extracted labeled tree from RST: {len(tree_string)} characters")
    return tree_string


def build_node_mapping(tree_file: Optional[Path] = None, rst_file: Optional[Path] = None) -> Tuple[List[TreeNode], Dict[int, TreeNode]]:
    """
    Build node ID to TreeNode mapping from PAML RST file.
    
    Parses the labeled Newick tree and assigns PAML node IDs to internal nodes.
    
    Args:
        tree_file: Path to Newick tree file (ignored, kept for API compatibility)
        rst_file: Path to RST file (required - contains labeled tree)
        
    Returns:
        Tuple of (ordered node list in postorder, node_id -> TreeNode dict)
        
    Raises:
        FileNotFoundError: If RST file doesn't exist
        ValueError: If tree parsing/building fails
    """
    if not rst_file or not rst_file.exists():
        raise FileNotFoundError(f"RST file required but not found: {rst_file}")
    
    logger.info(f"Building tree from RST labeled tree: {rst_file}")
    
    # Extract node range
    min_node, max_node = extract_paml_node_range(rst_file)
    expected_internal = max_node - min_node + 1
    logger.debug(f"Expected PAML node range: {min_node} to {max_node} ({expected_internal} nodes)")
    
    # Parse the labeled Newick tree
    tree_string = extract_rst_labeled_tree(rst_file)
    root = parse_newick(tree_string)
    
    # Assign IDs from PAML labels
    # Tips have format "pamlid_taxid" (with underscore)
    # Internal nodes have just the numeric ID (no underscore)
    def assign_ids_from_labels(node: TreeNode):
        """Recursively assign node_id from paml_label."""
        # Check if this is a tip (has underscore) or internal node (no underscore)
        if node.paml_label:
            if '_' in node.paml_label:
                # Tip: copy paml_label into name so downstream tip lookups work
                if node.name is None:
                    node.name = node.paml_label
            elif node.paml_label.isdigit():
                # This is an internal node - assign the ID
                node.node_id = int(node.paml_label)
        
        for child in node.children:
            assign_ids_from_labels(child)
    
    assign_ids_from_labels(root)
    
    # Build ordered list and ID mapping
    ordered_nodes = get_node_order(root)
    id_mapping = {n.node_id: n for n in ordered_nodes if n.node_id is not None}
    
    leaf_count = sum(1 for n in ordered_nodes if n.is_leaf())
    internal_count = len(id_mapping)
    
    logger.info(f"Built tree from RST labeled tree:")
    logger.info(f"  Total nodes: {len(ordered_nodes)}")
    logger.info(f"  Leaf nodes: {leaf_count}")
    logger.info(f"  Internal nodes with PAML IDs: {internal_count}")
    logger.info(f"  PAML node ID range: {min_node} to {max_node}")
    
    if internal_count != expected_internal:
        logger.warning(f"Expected {expected_internal} internal nodes, got {internal_count}")
    
    return ordered_nodes, id_mapping


def get_tip_labels(root: TreeNode) -> List[str]:
    """
    Extract all tip (leaf) labels from tree.
    
    Args:
        root: Root node of tree
        
    Returns:
        List of tip labels (e.g., taxids)
    """
    tips = []
    
    if root.is_leaf() and root.name:
        tips.append(root.name)
    else:
        for child in root.children:
            tips.extend(get_tip_labels(child))
    
    return tips


def find_node_by_name(root: TreeNode, name: str) -> Optional[TreeNode]:
    """
    Find a node by its name (tip label).
    
    Args:
        root: Root node of tree
        name: Name to search for
        
    Returns:
        TreeNode if found, None otherwise
    """
    if root.name == name:
        return root
    
    for child in root.children:
        result = find_node_by_name(child, name)
        if result:
            return result
    
    return None


def find_node_by_taxid(root: TreeNode, taxid: str) -> Optional[TreeNode]:
    """
    Find a tip node by taxid, where tips are labeled as 'lineage_taxid'.
    
    Args:
        root: Root node of tree
        taxid: Taxid string to search for (e.g., '9541')
        
    Returns:
        TreeNode if found, None otherwise
        
    Note:
        Tree tips from PAML have format 'lineage_taxid' (e.g., '53_9541').
        This function extracts the taxid part after the last underscore.
    """
    if root.is_leaf() and root.name:
        # Extract taxid from format 'lineage_taxid'
        tip_taxid = root.name.strip().split('_')[-1]
        if tip_taxid == taxid:
            return root
    
    for child in root.children:
        result = find_node_by_taxid(child, taxid)
        if result:
            return result
    
    return None


def get_mrca(root: TreeNode, tip_names: List[str]) -> Optional[TreeNode]:
    """
    Find most recent common ancestor (MRCA) of given tip nodes.
    
    Supports both exact name matching and taxid-based matching.
    For taxid matching, specify taxids as strings (e.g., '9541').
    Tree tips are expected to have format 'lineage_taxid'.
    
    Args:
        root: Root node of tree
        tip_names: List of tip labels or taxids to find MRCA for
        
    Returns:
        MRCA node, or None if not found
    """
    # Find all tip nodes - try exact match first, then taxid match
    tip_nodes = []
    for name in tip_names:
        name_str = str(name).strip()
        # Try exact match first
        node = find_node_by_name(root, name_str)
        if not node:
            # Try taxid match (for PAML trees with 'lineage_taxid' format)
            node = find_node_by_taxid(root, name_str)
        if node:
            tip_nodes.append(node)
    
    if len(tip_nodes) < 2:
        return None
    
    def _ancestors(node: Optional[TreeNode]) -> List[TreeNode]:
        anc = []
        while node:
            anc.append(node)
            node = node.parent
        return anc

    # Intersect ancestor sets by object identity
    common: Optional[set] = None
    for tip in tip_nodes:
        anc_set = set(_ancestors(tip))
        if common is None:
            common = anc_set
        else:
            common &= anc_set
        if not common:
            return None

    # Compute depth from root (root depth = 0)
    depth_cache: Dict[TreeNode, int] = {}
    def depth_from_root(node: TreeNode) -> int:
        if node in depth_cache:
            return depth_cache[node]
        d = 0
        cur = node
        while cur.parent is not None:
            d += 1
            cur = cur.parent
        depth_cache[node] = d
        return d

    # Choose the deepest common ancestor (max depth from root)
    if common is None:
        return None
    deepest = max(common, key=depth_from_root)
    return deepest
