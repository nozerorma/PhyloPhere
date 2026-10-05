"""
Posterior Probability Parsing
==============================

Parse ancestral state posteriors from PAML (.rst) files.
"""

import logging
import re
from pathlib import Path
from typing import Dict, Tuple, Union, Optional

from .tree_parser import build_node_mapping

logger = logging.getLogger(__name__)


_FULL_MARGINAL_SECTION_RE = re.compile(r"Prob distribution at node (\d+), by site")
_FULL_MARGINAL_SITE_LINE_RE = re.compile(r"^\s*(\d+)\s+\d+\s+\S+:\s*(.+)$", re.MULTILINE)
_FULL_MARGINAL_TOKEN_RE = re.compile(r"([A-Z*-])\(([\d.]+)\)")


def parse_paml_rst_full_marginal(
    content: str,
    threshold: float = 0.0,
    positions: Optional[set] = None,
) -> Dict[int, Dict[int, Dict[str, float]]]:
    """
    Parse PAML's per-node full marginal amino-acid distributions.

    Reads the "Prob distribution at node N, by site" sections that codeml
    writes when run with ``verbose = 2``: one section per internal node,
    each listing every site with its complete 20-amino-acid posterior (as
    opposed to "Prob of best state at each node", which keeps only the
    single most likely residue per site). Each section is self-labeled with
    its PAML node id, so no positional node-ordering is needed.

    Args:
        content: full text of a PAML rst file.
        threshold: minimum posterior probability for a residue to be kept.
            The residue with the highest probability at a site is always
            kept even if below threshold, so a site's distribution is never
            emptied out entirely.
        positions: optional set of 1-based site positions to keep; other
            sites are skipped.

    Returns:
        Dict mapping node_id -> position -> {AA: posterior}. Empty if the
        rst file has no "Prob distribution at node" sections (i.e. it was
        produced with ``verbose = 1``).
    """
    posteriors: Dict[int, Dict[int, Dict[str, float]]] = {}
    section_starts = list(_FULL_MARGINAL_SECTION_RE.finditer(content))
    for idx, section_marker in enumerate(section_starts):
        node_id = int(section_marker.group(1))
        section_end = (
            section_starts[idx + 1].start()
            if idx + 1 < len(section_starts)
            else len(content)
        )
        section = content[section_marker.end() : section_end]

        node_posteriors: Dict[int, Dict[str, float]] = {}
        for site_match in _FULL_MARGINAL_SITE_LINE_RE.finditer(section):
            site_num = int(site_match.group(1))
            if positions is not None and site_num not in positions:
                continue

            dist: Dict[str, float] = {}
            best_aa, best_prob = None, -1.0
            for aa, prob_str in _FULL_MARGINAL_TOKEN_RE.findall(site_match.group(2)):
                if aa in ("-", "*"):
                    continue
                prob = float(prob_str)
                if prob > best_prob:
                    best_aa, best_prob = aa, prob
                if prob >= threshold:
                    dist[aa] = prob
            if not dist and best_aa is not None:
                dist[best_aa] = best_prob
            if dist:
                node_posteriors[site_num] = dist

        posteriors[node_id] = node_posteriors

    return posteriors


def parse_paml_rst_node_level(
    rst_file: Path,
    tree_file: Path,
    threshold: float = 0.0,
    return_mapping: bool = False,
    positions: Optional[set] = None,
) -> Union[
    Dict[int, Dict[int, Dict[str, float]]],
    Tuple[Dict[int, Dict[int, Dict[str, float]]], Dict[int, int]],
]:
    """
    Parse PAML rst file and extract per-node posteriors.

    This function extracts node-level ancestral state posteriors from PAML
    RST output files. PAML lists nodes starting from the root and then outward
    (ascending node IDs). When the rst file contains full per-site marginal
    distributions (``verbose = 2`` in codeml.ctl), those are used directly via
    :func:`parse_paml_rst_full_marginal`; otherwise this falls back to parsing
    the single-best-state table.

    Args:
        rst_file: Path to PAML rst file
        tree_file: Path to Newick tree file (needed for node ordering)
        threshold: Minimum posterior probability to include (default 0.0)

    Args:
        rst_file: Path to PAML rst file
        tree_file: Path to Newick tree file (needed for node ordering)
        threshold: Minimum posterior probability to include (default 0.0)
        return_mapping: If True, also return mapping of local tree node IDs to
            PAML node numbers (as reported in the RST file)
        positions: Optional set of 1-based site positions to filter posteriors for.

    Returns:
        Dict mapping node_id -> position -> {AA: posterior}; if return_mapping=True
        returns a tuple of (posteriors, paml_to_internal_map) where the map is
        PAML node id -> internal tree node id used in posteriors.

    Raises:
        FileNotFoundError: If rst or tree file doesn't exist
        ValueError: If file format is invalid or no sites parsed

    Example:
        >>> posteriors = parse_paml_rst_node_level(
        ...     Path("output/rst"),
        ...     Path("tree.nwk")
        ... )
        >>> # Get state at node 5, position 10
        >>> node5_pos10 = posteriors[5][10]  # {'A': 0.95, 'R': 0.03}
        >>> most_likely = max(node5_pos10, key=node5_pos10.get)  # 'A'
    """
    if not rst_file.exists():
        raise FileNotFoundError(f"RST file not found: {rst_file}")

    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    # Parse tree to get node IDs. PAML marginal reconstruction lists nodes
    # starting at the root and then proceeding outward; IDs are contiguous
    # ascending (e.g., TP53: root=187, children 188/322, max=371). The order
    # of posterior tokens follows this ascending PAML node_id ordering.
    logger.info(f"Extracting labeled tree from RST file: {rst_file}")
    ordered_nodes, id_mapping = build_node_mapping(
        tree_file=tree_file, rst_file=rst_file
    )
    internal_nodes = [node for node in ordered_nodes if not node.is_leaf()]
    num_internal = len(internal_nodes)
    logger.info(
        f"Tree has {num_internal} internal nodes (out of {len(ordered_nodes)} total)"
    )

    # Order internal nodes the way PAML reports posteriors: ascending node_id.
    if not all(node.node_id is not None for node in internal_nodes):
        raise ValueError(
            "All internal nodes must have PAML node IDs for posterior parsing."
        )

    # Only consider nodes that have PAML node IDs assigned
    internal_nodes_for_posteriors = sorted(
        [n for n in internal_nodes if getattr(n, "node_id", None) is not None],
        key=lambda n: int(n.node_id) if n.node_id is not None else -1,
    )

    posteriors: Dict[int, Dict[int, Dict[str, float]]] = {
        node.node_id: {}
        for node in internal_nodes_for_posteriors
        if node.node_id is not None
    }
    # Map local sequential index (PAML order) to PAML node ID from the parsed tree
    node_id_map: Dict[int, int] = {
        idx: node.node_id
        for idx, node in enumerate(internal_nodes_for_posteriors)
        if node.node_id is not None
    }

    # Read RST file
    with open(rst_file, "r") as f:
        content = f.read()

    if "Prob distribution at node" in content:
        posteriors = parse_paml_rst_full_marginal(
            content, threshold=threshold, positions=positions
        )
        logger.info(
            f"Parsed full marginal posteriors for {len(posteriors)} nodes from {rst_file}"
        )
        if return_mapping:
            return posteriors, node_id_map
        return posteriors

    # PAML node IDs are already extracted from the labeled tree in build_node_mapping
    logger.debug(
        f"Using PAML node IDs from labeled tree: {sorted(id_mapping.keys()) if id_mapping else 'none'}"
    )

    # Find the marginal reconstruction section
    start_marker = "(1) Marginal reconstruction"
    start_idx = content.find(start_marker)
    if start_idx == -1:
        raise ValueError(
            f"Could not find '{start_marker}' in {rst_file}. "
            f"This may not be a valid PAML rst file."
        )

    # Find the "Prob of best state at each node" section
    prob_marker = "Prob of best state at each node"
    prob_idx = content.find(prob_marker, start_idx)
    if prob_idx == -1:
        raise ValueError(f"Could not find '{prob_marker}' section in {rst_file}")

    # Find data section start
    data_header_marker = "site   Freq   Data:"
    data_start_idx = content.find(data_header_marker, prob_idx)
    if data_start_idx == -1:
        raise ValueError(
            f"Could not find data header '{data_header_marker}' in {rst_file}"
        )

    # Extract data section
    data_section_start = data_start_idx + len(data_header_marker)

    # Find section end
    end_markers = ["(2)", "List of extant", "TREE #", "Nodes", "tree with node"]
    data_section_end = len(content)
    for marker in end_markers:
        idx = content.find(marker, data_section_start)
        if idx != -1:
            data_section_end = min(data_section_end, idx)

    data_section = content[data_section_start:data_section_end]
    lines = data_section.split("\n")

    # Parse site entries
    i = 0
    sites_parsed = 0

    while i < len(lines):
        line = lines[i]
        line_stripped = line.strip()

        # Skip empty lines
        if not line_stripped:
            i += 1
            continue

        # Check for site number pattern
        match = re.match(r"^\s*(\d+)\s+(\d+)\s+(.*)", line)
        if match:
            site_num = int(match.group(1))
            rest = match.group(3)

            # Collect full site data (may span multiple lines)
            site_data = rest
            i += 1
            while i < len(lines) and ":" not in site_data:
                site_data += lines[i]
                i += 1

            # If positions filter is provided, skip site entries not in set
            if positions is not None and site_num not in positions:
                # Skip parsing the tokens for this site
                # but continue scanning the file so line pointer is advanced
                while i < len(lines) and ":" not in site_data:
                    site_data += lines[i]
                    i += 1
                # Continue main loop - do not store data for this site
                continue

            # Parse if colon found
            if ":" in site_data:
                parts = site_data.split(":", 1)
                consensus_str = parts[0].strip()
                posterior_str = parts[1].strip() if len(parts) > 1 else ""

                # Remove trailing markers
                posterior_str = posterior_str.rstrip("$").strip()

                # Parse posterior tokens
                posterior_tokens = posterior_str.split()

                # Validate lengths
                if len(posterior_tokens) != len(consensus_str):
                    min_len = min(len(consensus_str), len(posterior_tokens))
                    posterior_tokens = posterior_tokens[:min_len]
                    consensus_str = consensus_str[:min_len]

                usable_tokens = posterior_tokens
                if len(usable_tokens) != num_internal:
                    logger.debug(
                        f"Site {site_num}: Expected {num_internal} internal nodes, "
                        f"got {len(usable_tokens)} posteriors. Using min length."
                    )
                    usable_tokens = usable_tokens[
                        : min(len(usable_tokens), num_internal)
                    ]

                for posterior_idx, token in enumerate(usable_tokens):
                    node = internal_nodes_for_posteriors[posterior_idx]
                    node_id = node.node_id
                    if node_id is None:
                        continue
                    record_id = node_id

                    # Parse token: "F(0.995)"
                    if not token or "(" not in token or ")" not in token:
                        logger.debug(
                            f"Skipping malformed token '{token}' at site {site_num} posterior {posterior_idx}"
                        )
                        continue

                    try:
                        # Extract AA and probability
                        paren_idx = token.find("(")
                        aa = token[:paren_idx]
                        prob_str = token[paren_idx + 1 : -1]
                        prob = float(prob_str)

                        # Always keep the modal/best state; threshold only prunes secondary AAs
                        node_entry = posteriors.setdefault(record_id, {})
                        site_entry = node_entry.setdefault(site_num, {})
                        if prob >= threshold or not site_entry:
                            site_entry[aa] = prob

                    except (ValueError, IndexError) as e:
                        logger.debug(
                            f"Failed to parse token '{token}' at site {site_num} node {node_id}: {e}"
                        )
                        continue

                sites_parsed += 1
            else:
                logger.warning(f"Site {site_num}: No colon separator found")
                i += 1
        else:
            i += 1

    if sites_parsed == 0:
        raise ValueError(
            f"No sites parsed from {rst_file}. "
            f"Ensure this is a PAML rst file with marginal reconstruction."
        )

    # Count nodes with data
    nodes_with_data = sum(1 for node_id in posteriors if posteriors[node_id])

    logger.info(
        f"Parsed node-level posteriors for {sites_parsed} sites from {rst_file}"
    )
    logger.info(f"Tree has {len(ordered_nodes)} total nodes")
    logger.info(
        f"Extracted posteriors for {nodes_with_data} internal (ancestral) nodes"
    )
    logger.info("NOTE: Tip nodes (observed data) do not have posteriors in PAML output")

    if return_mapping:
        return posteriors, node_id_map
    return posteriors
