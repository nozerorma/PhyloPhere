"""
Convergence Classification Module
=================================

Classifies amino acid substitution patterns as divergent, convergent,
or parallel evolution using ancestral state reconstruction (ASR) data.
Supports dynamic multi-pair focal logic for contrast-based analyses.

Key Concepts
------------
**Convergent**: Multiple lineages evolve the same derived state from
    a shared ancestral state (e.g., A→V in both groups).
**Divergent**: Lineages evolve different derived states from the same
    ancestral state (e.g., A→V vs A→I).
**Parallel**: Both groups change but to different states, with different
    ancestral states (e.g., A→V and T→I).
**Codivergent**: Mixed pattern where one group converges while the other diverges.

Classification Hierarchy
-----------------------
1. **no_change**: No substitutions detected in any group
2. **convergent**: One or both groups converge to a single derived state
3. **divergent**: One group diverges to multiple states, other unchanged/insufficient
4. **parallel_convergence**: Both groups converge to different states
5. **parallel_divergence**: Both groups diverge to different state sets
6. **codivergent**: Mixed convergence + divergence patterns
7. **ambiguous**: Complex pattern requiring manual inspection

Data Contracts
--------------
- **PairDetail**: TypedDict for tip-level pair data (pair_id, focal_state, tip residues)
- **NodeStates**: Dataclass container for key phylogenetic node states
- **ConvergenceClassification**: Dataclass for classification results

Usage Example
-------------
::

    from src.convergence.convergence import classify_tip_level_pattern

    pair_details = [
        {'pair_id': 'pair_1', 'focal_state': 'A',
         'top_tip_mode': 'V', 'bottom_tip_mode': 'V'},
        {'pair_id': 'pair_2', 'focal_state': 'A',
         'top_tip_mode': 'V', 'bottom_tip_mode': 'V'}
    ]
    result = classify_tip_level_pattern(pair_details, convergence_mode='focal_clade')
    print(result['pattern'])  # 'convergent'
    print(result['description'])  # 'Convergent (TOP only): → V, BOTTOM unchanged'

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-07
"""

from typing import Dict, List, Optional, Tuple, Set, Any, Sequence, TypedDict
from dataclasses import dataclass, field
import logging

from ..biochem.grouping import get_grouping_scheme

logger = logging.getLogger(__name__)


class PairDetail(TypedDict, total=False):
    """Tip-level pair detail contract used for convergence classification."""

    pair_id: Optional[str]
    focal_state: Optional[str]
    mrca_contrast: Optional[str]
    top_tip_mode: Optional[str]
    bottom_tip_mode: Optional[str]
    top_tip_residue: Optional[str]
    bottom_tip_residue: Optional[str]
    # Additional diagnostic or probability fields are tolerated


@dataclass
class NodeStates:
    """
    Container for amino acid states at key phylogenetic nodes.

    **REDESIGNED FOR DYNAMIC PAIRS (2025-12-05)**:
    focal_states and focal_probs are now lists supporting variable pairs (1-N).

    **FALLBACK TRACKING (2025-12-07)**:
    fallback_depths tracks how many ancestor levels were traversed for each focal node.

    Attributes:
        root (Optional[str]): State at root
        pre_split (Optional[str]): State at internal node before first split
        mrca_contrast (str): State at MRCA of focal contrast
        position (int): Alignment position (1-based)
        gene (str): Gene identifier
        focal_states (List[str]): States at focal lineage MRCAs [focal_1..focal_N]
        focal_probs (List[float]): Posterior probabilities for focal states
        fallback_depths (Dict[str, int]): Depth of ancestor fallback for each node (focal_1, focal_2, root, etc.)
    """

    root: Optional[str]
    pre_split: Optional[str]
    mrca_contrast: str
    position: int
    gene: str
    focal_states: List[str] = field(default_factory=list)
    focal_probs: List[float] = field(default_factory=list)
    root_prob: Optional[float] = None
    pre_split_prob: Optional[float] = None
    mrca_contrast_prob: Optional[float] = None
    low_confidence_nodes: Optional[List[str]] = None
    fallback_depths: Optional[Dict[str, int]] = None


@dataclass
class ConvergenceClassification:
    """
    Result of convergence pattern classification.

    Attributes:
        is_divergent (bool): True if pattern shows divergent evolution
        is_convergent (bool): True if pattern shows convergent evolution
        is_parallel (bool): True if pattern shows parallel evolution
        is_codivergent (bool): True if one group converges, other diverges
        is_caap (str): Comma-separated list of converging schemes (e.g., "GS3,GS4")
        confidence (float): Classification confidence (0.0-1.0)
        ancestral_state (str): Reconstructed ancestral amino acid
        derived_states (List[str]): Derived amino acids in focal lineages
        pattern_description (str): Human-readable pattern description
        ambiguous (bool): True if pattern cannot be unambiguously classified
        individual_transitions (Optional[Dict]): Per-transition metrics for divergent/parallel
        top_convergent (bool): True if top group shows convergence (tip-level)
        bottom_convergent (bool): True if bottom group shows convergence (tip-level)
        top_changed (List[str]): Descendant states in top group (tip-level)
        bottom_changed (List[str]): Descendant states in bottom group (tip-level)
    """

    is_divergent: bool
    is_convergent: bool
    is_parallel: bool
    is_codivergent: bool
    is_caap: str
    confidence: float
    ancestral_state: str
    derived_states: List[str]
    pattern_description: str
    ambiguous: bool = False
    individual_transitions: Optional[Dict[str, Any]] = None
    top_convergent: bool = False
    bottom_convergent: bool = False
    top_changed: List[str] = field(default_factory=list)
    bottom_changed: List[str] = field(default_factory=list)


# =============================================================================
# TIP-LEVEL ANALYSIS FUNCTIONS
# =============================================================================


def build_alignment_lookup(
    alignment, taxid_to_species: Optional[Dict[str, str]] = None
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Build lookup dictionaries for fast alignment access by taxid or species name.

    Creates two complementary dictionaries:
    - seq_by_id: Maps alignment record IDs (typically taxids) to sequences
    - seq_by_species: Maps species names to sequences

    **Dual-key lookup strategy:**
    This function supports flexible species identification by indexing sequences
    under both taxid (record ID) and species name. When taxid_to_species mapping
    is provided, each sequence is indexed under both keys. When unavailable,
    the record ID is used as both the taxid and species name.

    This dual strategy handles:
    - Alignment files where record IDs are taxids (common in phylogenetic data)
    - Alignment files where record IDs are species names (common in MSA tools)
    - Trait files that reference species by name or taxid

    The lookup functions (_fetch_residue_from_lookup, collect_tip_residues) try
    taxid lookup first, then fall back to species name lookup. This is NOT a
    failure fallback - it's an intentional dual-key strategy to handle different
    identifier formats without requiring normalization.

    Args:
        alignment: BioPython MultipleSeqAlignment object
        taxid_to_species: Optional mapping from taxid (str) -> species name (str)

    Returns:
        Tuple of (seq_by_id, seq_by_species) dictionaries

    Example:
        >>> seq_by_id, seq_by_species = build_alignment_lookup(alignment, taxid_map)
        >>> residue = seq_by_id['9606'][position]  # Access by taxid
        >>> residue = seq_by_species['Homo_sapiens'][position]  # Access by species
    """
    seq_by_id = {}
    seq_by_species = {}
    taxid_to_species = taxid_to_species or {}

    for rec in alignment:
        rec_id = str(rec.id)
        seq = str(rec.seq)
        seq_by_id[rec_id] = seq

        # Also index by species name if available
        species_name = taxid_to_species.get(rec_id)
        if species_name:
            seq_by_species[species_name] = seq
        else:
            # Also index by record ID (handles alignments where ID is species name)
            seq_by_species[rec_id] = seq

    return seq_by_id, seq_by_species


def _fetch_residue_from_lookup(
    label: str, msa_pos: int, seq_by_id: Dict[str, str], seq_by_species: Dict[str, str]
) -> Optional[str]:
    """
    Fetch residue from alignment lookup dictionaries using dual-key strategy.

    Tries taxid lookup first (seq_by_id), then species name lookup (seq_by_species).
    This handles trait files that reference species by either identifier format.

    Args:
        label: Taxid or species name to look up
        msa_pos: Alignment position (0-based)
        seq_by_id: Dictionary mapping taxid -> sequence
        seq_by_species: Dictionary mapping species name -> sequence

    Returns:
        Single-letter amino acid code, or None if not found
    """
    if label is None:
        return None

    key = str(label)

    # Try ID lookup first (most common: taxid in trait file matches record ID)
    seq = seq_by_id.get(key)
    if seq is None:
        # Try species name lookup (handles trait files with species names instead of taxids)
        seq = seq_by_species.get(key)
        if seq is not None:
            logger.debug(f"Resolved '{key}' via species name lookup (not found in taxid lookup)")

    if seq and 0 <= msa_pos < len(seq):
        return seq[msa_pos]

    return None


def collect_tip_residues(
    taxa_list: List[str],
    species_list: List[str],
    msa_pos: int,
    seq_by_id: Dict[str, str],
    seq_by_species: Dict[str, str],
    taxid_to_species: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """
    Collect residues observed at tip species for a trait-specific group.

    Args:
        taxa_list: List of taxids for this group
        species_list: List of species names for this group
        msa_pos: Alignment position (0-based)
        seq_by_id: Dictionary mapping taxid -> sequence
        seq_by_species: Dictionary mapping species name -> sequence
        taxid_to_species: Optional mapping from taxid to species name

    Returns:
        List of dicts with keys: 'taxid', 'species', 'residue'

    Example:
        >>> tip_records = collect_tip_residues(
        ...     ['9606', '9598'],
        ...     ['Homo_sapiens', 'Pan_troglodytes'],
        ...     85, seq_by_id, seq_by_species
        ... )
        >>> print(tip_records)
        [{'taxid': '9606', 'species': 'Homo_sapiens', 'residue': 'A'},
         {'taxid': '9598', 'species': 'Pan_troglodytes', 'residue': 'A'}]
    """
    tip_records = []
    taxid_to_species = taxid_to_species or {}
    species_list = species_list or []

    # Process taxids
    for idx, taxid in enumerate(taxa_list or []):
        species_name = taxid_to_species.get(str(taxid))
        if not species_name and idx < len(species_list):
            species_name = species_list[idx]

        # Try to fetch residue by taxid first
        residue = _fetch_residue_from_lookup(taxid, msa_pos, seq_by_id, seq_by_species)

        # Retry with species name if taxid lookup failed (handles trait files with species names)
        if residue is None and species_name:
            residue = _fetch_residue_from_lookup(species_name, msa_pos, seq_by_id, seq_by_species)
            if residue is not None:
                logger.debug(f"Resolved taxid '{taxid}' via species name '{species_name}'")

        tip_records.append({"taxid": str(taxid), "species": species_name, "residue": residue})

    # If only species names provided (no taxids)
    if not (taxa_list or []) and species_list:
        for species_name in species_list:
            residue = _fetch_residue_from_lookup(species_name, msa_pos, seq_by_id, seq_by_species)
            tip_records.append({"taxid": None, "species": species_name, "residue": residue})

    return tip_records


def extract_tip_residue(tip_records: List[Dict[str, Any]]) -> Optional[str]:
    """
    Get the residue among tip species.

    Args:
        tip_records: List of tip residue records from collect_tip_residues()

    Returns:
        Amino acid residue among tip species, or None if no valid residues

    Example:
        >>> tips = [
        ...     {'taxid': '9606', 'species': 'Homo_sapiens', 'residue': 'A'},
        ...     {'taxid': '9598', 'species': 'Pan_troglodytes', 'residue': 'A'},
        ...     {'taxid': '9544', 'species': 'Macaca_mulatta', 'residue': 'T'}
        ... ]
        >>> extract_tip_residue(tips)
        'A'
    """
    # Filter out missing/gap/unknown residues
    residues = [tip["residue"] for tip in tip_records if tip.get("residue") not in {None, "-", "X"}]

    if not residues:
        return None

    # Return the residue
    return residues[0] if residues else None


def normalize_amino_list(values: List[Optional[str]]) -> List[str]:
    """
    Clean and deduplicate a list of amino acids.

    Args:
        values: List of amino acid codes (may contain None, gaps, etc.)

    Returns:
        Cleaned list of unique amino acids (uppercase, no gaps/unknowns)

    Example:
        >>> normalize_amino_list(['A', 'a', 'T', None, '-', 'X', 'A'])
        ['A', 'T']
    """
    cleaned = []
    seen = set()

    for value in values or []:
        if value is None:
            continue

        aa = str(value).strip().upper()

        # Skip empty, gaps, unknowns, and duplicates
        if aa in {"", "-", "X", "?"} or aa in seen:
            continue

        seen.add(aa)
        cleaned.append(aa)

    return cleaned


def format_amino_display(amino_list: List[str]) -> str:
    """
    Format amino acid list for display.

    Args:
        amino_list: List of amino acid codes

    Returns:
        Formatted string (e.g., "A+T", "H+R+K")

    Example:
        >>> format_amino_display(['A', 'T'])
        'A+T'
        >>> format_amino_display([])
        '?'
    """
    if not amino_list:
        return "?"
    return "+".join(amino_list)


def compute_derived_state_similarity(derived_states: List[str], grouping_getter) -> Dict[str, Any]:
    """
    Check if derived states belong to the same biochemical group.

    For divergent/parallel patterns, this measures whether the
    different derived states are biochemically similar.

    Args:
        derived_states: List of derived amino acid states

    Returns:
        Dictionary with per-scheme similarity analysis

    Example:
        {
            'GS1': {'derived_groups': ['AGPS', 'CV'], 'convergent': False, ...},
            ...
        }
    """
    summary = {}

    if grouping_getter is None or not derived_states:
        return summary

    for idx in range(5):
        scheme = f"GS{idx}"

        # Get groups for all derived states
        groups = []
        seen = set()
        for aa in derived_states:
            group = grouping_getter(aa, scheme)
            if group and group not in seen:
                seen.add(group)
                groups.append(group)

        # Convergent if all derived states in same group
        convergent = len(groups) == 1 if groups else False

        # Description
        if convergent:
            description = f"All in {groups[0]}"
        else:
            description = f"Across {len(groups)} groups: {'/'.join(groups)}"

        summary[scheme] = {
            "derived_groups": groups,
            "convergent": convergent,
            "description": description,
        }

    return summary


# =============================================================================
# NODE-LEVEL ANALYSIS FUNCTIONS (MRCA-based)
# =============================================================================


def extract_node_states_from_node_level(
    asr_posteriors: Dict[int, Dict[int, Dict[str, float]]],
    node_mapping: Optional[Dict[str, Any]],
    position: int,
    gene: str,
    posterior_threshold: float = 0.7,
    tree_node_lookup: Optional[Dict[int, Any]] = None,
) -> Optional[NodeStates]:
    """
    Extract amino acid states at key phylogenetic nodes from node-level ASR posteriors.

    **UPDATED FOR DYNAMIC PAIRS (2025-12-05)**:
    Now extracts variable numbers of focal states from node_mapping['focal_nodes'] list.

    Args:
        asr_posteriors: Dict mapping node_id (int) -> position -> amino_acid -> probability
        node_mapping: Dict with 'root', 'mrca_contrast', 'focal_nodes' (list of node IDs)
        position: Alignment position (1-based)
        gene: Gene identifier
        posterior_threshold: Minimum posterior probability to accept state
        tree_node_lookup: Optional mapping of node_id -> TreeNode (with parent links) to
            enable fallback to ancestral nodes when focal nodes are missing or have
            sub-threshold posteriors.

    Returns:
        NodeStates object or None if extraction fails or confidence too low
    """

    # Guard against missing or incomplete node mappings
    if not node_mapping:
        logger.debug(
            f"No node mapping provided for {gene}:{position}; skipping node-level extraction"
        )
        return None

    required_roles = {"root", "mrca_contrast", "focal_nodes"}
    missing_roles = required_roles - set(node_mapping.keys())
    if missing_roles:
        logger.debug(
            f"Incomplete node mapping for {gene}:{position}; missing roles {missing_roles}"
        )
        return None

    def get_modal_state(node_id: int) -> Optional[Tuple[str, float]]:
        """Get most probable amino acid and its posterior probability."""
        if node_id not in asr_posteriors:
            logger.debug(f"Node {node_id} not found in ASR posteriors")
            return None

        if position not in asr_posteriors[node_id]:
            logger.debug(f"Position {position} not found in posteriors for node {node_id}")
            return None

        posteriors = asr_posteriors[node_id][position]
        if not posteriors:
            return None

        modal_aa = max(posteriors.items(), key=lambda x: x[1])
        return modal_aa[0], modal_aa[1]

    def resolve_with_parent_fallback(
        start_id: Optional[int], idx: int
    ) -> Tuple[Optional[Tuple[str, float]], bool, int]:
        """Return modal state for node, falling back up the tree if needed (max 3 levels).

        Returns:
            (state, fallback_used, fallback_depth) where fallback_depth is number of ancestors traversed
        """
        if start_id is None:
            return None, False, 0

        current_id = start_id
        visited: Set[int] = set()
        fallback_used = False
        fallback_depth = 0
        max_fallback_depth = 5
        fallback_trace: List[str] = []  # Track all attempts for logging

        while (
            current_id is not None
            and current_id not in visited
            and fallback_depth <= max_fallback_depth
        ):
            visited.add(current_id)
            modal = get_modal_state(current_id)

            # Record this level for trace logging
            if modal:
                fallback_trace.append(
                    f"node={current_id}:residue={modal[0]}:posterior={modal[1]:.4f}"
                )
            else:
                fallback_trace.append(f"node={current_id}:residue=None:posterior=0.0000")

            if modal and modal[1] >= posterior_threshold:
                if fallback_used and start_id != current_id:
                    logger.info(
                        "Focal node %d missing/low confidence; fell back to ancestor %d (depth=%d). "
                        "Trace: %s",
                        idx + 1,
                        current_id,
                        fallback_depth,
                        " -> ".join(fallback_trace),
                    )
                return modal, fallback_used, fallback_depth

            # Move to parent if available (but respect max fallback depth)
            if not tree_node_lookup or fallback_depth >= max_fallback_depth:
                break
            node_obj = tree_node_lookup.get(current_id)
            parent_id = node_obj.parent.node_id if node_obj and node_obj.parent else None
            if parent_id is None:
                break
            fallback_used = True
            fallback_depth += 1
            current_id = parent_id

        # Log final unresolved trace
        if fallback_trace:
            logger.warning(
                "Focal node %d unresolved after fallback (depth=%d). Trace: %s",
                idx + 1,
                fallback_depth,
                " -> ".join(fallback_trace),
            )

        return None, fallback_used, fallback_depth

    low_conf_nodes: List[str] = []
    fallback_depths: Dict[str, int] = {}  # Track fallback depth for each node

    try:
        root_data = get_modal_state(node_mapping["root"])
        mrca_data = get_modal_state(node_mapping["mrca_contrast"])

        # Extract focal states dynamically from focal_nodes list
        focal_node_ids = node_mapping.get("focal_nodes", [])
        focal_states_data = []
        for idx, focal_id in enumerate(focal_node_ids):
            focal_data, used_fallback, fallback_depth = resolve_with_parent_fallback(focal_id, idx)
            fallback_depths[f"focal_{idx+1}"] = fallback_depth
            if focal_data is None:
                if focal_id is None:
                    logger.warning(f"Focal node {idx+1} is None for {gene}:{position}")
                else:
                    logger.warning(
                        "Focal node %d unresolved for %s:%s after ancestor fallback",
                        idx + 1,
                        gene,
                        position,
                    )
            elif used_fallback:
                low_conf_nodes.append(f"focal_{idx+1}")
            focal_states_data.append(focal_data)

        if not root_data or not mrca_data or not any(focal_states_data):
            logger.debug(f"Missing required node states for {gene}:{position}")
            return None

        if root_data and root_data[1] < posterior_threshold:
            low_conf_nodes.append("root")
        if mrca_data and mrca_data[1] < posterior_threshold:
            low_conf_nodes.append("mrca_contrast")

        for idx, data in enumerate(focal_states_data, 1):
            if data and data[1] < posterior_threshold:
                low_conf_nodes.append(f"focal_{idx}")

        # Extract optional nodes
        pre_split_data = (
            get_modal_state(node_mapping.get("pre_split", -1))
            if "pre_split" in node_mapping
            else None
        )

        # Build lists of focal states and probabilities
        focal_states = [d[0] if d else "" for d in focal_states_data]
        focal_probs = [d[1] if d else 0.0 for d in focal_states_data]

        return NodeStates(
            root=root_data[0] if root_data else None,
            pre_split=pre_split_data[0] if pre_split_data else None,
            mrca_contrast=mrca_data[0] if mrca_data else "",
            focal_states=focal_states,
            focal_probs=focal_probs,
            position=position,
            gene=gene,
            root_prob=root_data[1] if root_data else None,
            pre_split_prob=pre_split_data[1] if pre_split_data else None,
            mrca_contrast_prob=mrca_data[1] if mrca_data else None,
            low_confidence_nodes=low_conf_nodes or None,
            fallback_depths=fallback_depths if fallback_depths else None,
        )

    except Exception as e:
        logger.error(f"Error extracting node states for {gene}:{position}: {e}")
        return None


def classify_tip_level_pattern(
    pair_details: Sequence[PairDetail],
    convergence_mode: str = "focal_clade",
    grouping_scheme: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Classify evolutionary pattern using tip-level transitions and ancestral states.

    :param pair_details: Iterable of tip-level pair detail dicts containing
        ``pair_id``, ``focal_state`` or ``mrca_contrast``, and tip residues
        (``top_tip_mode``/``top_tip_residue`` and ``bottom_tip_mode``/``bottom_tip_residue``).
    :param convergence_mode: ``'mrca'`` to use MRCA contrast as ancestor for both sides,
        ``'focal_clade'`` (default) to use focal-specific ancestors.
    :returns: Classification dictionary with pattern label, description, change flags,
        derived state lists, and per-transition diagnostics.

    The function validates input shape, sorts pairs deterministically by ``pair_id``,
    and requires at least two usable pairs.

    Patterns:
        - convergent: One group converges, other unchanged
        - parallel: Both groups converge to different states
        - codivergent_top_converges: TOP converges, BOTTOM diverges
        - codivergent_bottom_converges: BOTTOM converges, TOP diverges
        - divergent: One group diverges, other unchanged
        - ambiguous: Complex pattern requiring manual inspection
        - no_change: No substitutions detected

    Args:
        pair_details: List of pair dicts with keys:
            - focal_k: Ancestral state at focal node k (dynamic list)
            - mrca_contrast: Ancestral state at contrast MRCA
            - top_tip_mode: Residue in top group
            - bottom_tip_mode: Residue in bottom group
        convergence_mode: 'mrca' to use mrca_contrast as ancestral state,
                         'focal_clade' to use per-focal-node ancestry (dynamic, default)

    Returns:
        Dictionary with classification results

    Example:
        >>> pair_details = [
        ...     {'focal_1': 'A', 'top_tip_mode': 'A', 'bottom_tip_mode': 'T'},
        ...     {'focal_2': 'A', 'top_tip_mode': 'A', 'bottom_tip_mode': 'T'}
        ... ]
        >>> result = classify_tip_level_pattern(pair_details)
        >>> print(result['pattern'])
        'convergent'
    """

    def _validate_pair_details(entries: Sequence[PairDetail]) -> Optional[List[PairDetail]]:
        """Validate and order pair details; return None when insufficient."""
        if not entries:
            return None

        usable: List[PairDetail] = []
        for entry in entries:
            top_tip = entry.get("top_tip_mode") or entry.get("top_tip_residue")
            bottom_tip = entry.get("bottom_tip_mode") or entry.get("bottom_tip_residue")
            if top_tip is None and bottom_tip is None:
                logger.debug(f"Skipping pair with no tip residues: {entry}")
                continue
            usable.append(entry)

        if len(usable) < 2:
            return None

        # Deterministic ordering: sort by pair_id (case-insensitive) with stable index fallback
        indexed = list(enumerate(usable))
        indexed.sort(key=lambda item: (str(item[1].get("pair_id", "")).lower(), item[0]))
        return [item[1] for item in indexed]

    validated_pairs = _validate_pair_details(pair_details)
    if validated_pairs is None:
        return {
            "pattern": "insufficient_data",
            "description": "Need at least 2 pairs with tip residues for classification",
            "top_convergent": False,
            "bottom_convergent": False,
        }

    invalid_states = {None, "-", "X", "?"}

    def is_valid(state):
        return state not in invalid_states

    top_transitions: List[Dict[str, Any]] = []
    bottom_transitions: List[Dict[str, Any]] = []
    ancestor_states: Set[str] = set()

    # Log which convergence mode is being used
    logger.info(
        f"=== classify_tip_level_pattern CALLED with convergence_mode='{convergence_mode}' ==="
    )
    logger.debug(f"Number of pairs: {len(validated_pairs)}")

    def _map_state(state: Optional[str]) -> Optional[str]:
        if state is None:
            return None
        aa = str(state).strip().upper()
        if not aa:
            return None
        if not grouping_scheme or grouping_scheme == "US":
            return aa
        mapped = get_grouping_scheme(aa, grouping_scheme)
        return mapped if mapped else aa

    for pair in validated_pairs:
        # Choose ancestral state based on convergence_mode
        pair_id = pair.get("pair_id", "unknown")
        focal_state = pair.get("focal_state")
        mrca_contrast = pair.get("mrca_contrast")

        logger.debug(
            f"Pair {pair_id} available states: focal_state={focal_state}, mrca_contrast={mrca_contrast}"
        )

        if convergence_mode == "mrca":
            # Use mrca_contrast as the common ancestor for both focal groups
            ancestor = mrca_contrast
            logger.info(f"Pair {pair_id}: MODE='mrca' -> USING ancestor={ancestor} (mrca_contrast)")
        else:  # focal_clade mode (default)
            # Use focal-specific ancestors (focal_i aligned to pair order) for top/bottom
            ancestor = focal_state
            logger.info(
                f"Pair {pair_id}: MODE='focal_clade' -> USING ancestor={ancestor} (focal_state)"
            )

        top_tip = pair.get("top_tip_mode") or pair.get("top_tip_residue")
        bottom_tip = pair.get("bottom_tip_mode") or pair.get("bottom_tip_residue")

        mapped_ancestor = _map_state(ancestor)
        mapped_top_tip = _map_state(top_tip)
        mapped_bottom_tip = _map_state(bottom_tip)

        if is_valid(mapped_ancestor):
            if ancestor is not None:
                ancestor_states.add(str(mapped_ancestor))

        top_changed_flag = False
        bottom_changed_flag = False

        if is_valid(mapped_ancestor) and is_valid(mapped_top_tip):
            top_changed_flag = mapped_ancestor != mapped_top_tip
            top_transitions.append(
                {
                    "pair_id": pair.get("pair_id"),
                    "ancestor": mapped_ancestor,
                    "descendant": mapped_top_tip,
                    "changed": top_changed_flag,
                }
            )

        if is_valid(mapped_ancestor) and is_valid(mapped_bottom_tip):
            bottom_changed_flag = mapped_ancestor != mapped_bottom_tip
            bottom_transitions.append(
                {
                    "pair_id": pair.get("pair_id"),
                    "ancestor": mapped_ancestor,
                    "descendant": mapped_bottom_tip,
                    "changed": bottom_changed_flag,
                }
            )

    logger.debug(f"Top transitions: {top_transitions}")
    logger.debug(f"Bottom transitions: {bottom_transitions}")
    logger.debug(f"Ancestor states: {ancestor_states}")

    top_changed = [t["descendant"] for t in top_transitions if t["changed"]]
    bottom_changed = [t["descendant"] for t in bottom_transitions if t["changed"]]

    # Count actual amino acid changes (number of substitution events, not unique states)
    top_change_count = len(top_changed)
    bottom_change_count = len(bottom_changed)

    # Require ≥2 changes to classify as convergent/divergent
    # With <2 changes, we cannot reliably assess the pattern
    top_sufficient_changes = top_change_count >= 2
    bottom_sufficient_changes = bottom_change_count >= 2

    # Determine change types based on unique derived states and change count
    top_unique_states = set(top_changed) if top_changed else set()
    bottom_unique_states = set(bottom_changed) if bottom_changed else set()

    # Classify change types (convergent, divergent, codivergent, or no_change)
    top_convergent = len(top_unique_states) == 1 and top_sufficient_changes
    bottom_convergent = len(bottom_unique_states) == 1 and bottom_sufficient_changes

    top_divergent = len(top_unique_states) > 1 and top_sufficient_changes
    bottom_divergent = len(bottom_unique_states) > 1 and bottom_sufficient_changes

    # Codivergent: both convergence and divergence in same group
    # This happens when some pairs converge to one state while others converge to another
    top_codivergent = (
        len(top_unique_states) > 1
        and top_sufficient_changes
        and any(top_changed.count(state) > 1 for state in top_unique_states)
    )
    bottom_codivergent = (
        len(bottom_unique_states) > 1
        and bottom_sufficient_changes
        and any(bottom_changed.count(state) > 1 for state in bottom_unique_states)
    )

    top_unchanged = top_change_count == 0
    bottom_unchanged = bottom_change_count == 0

    # Insufficient changes: cannot assess pattern
    top_insufficient = top_change_count > 0 and not top_sufficient_changes
    bottom_insufficient = bottom_change_count > 0 and not bottom_sufficient_changes

    top_descendants = sorted(top_unique_states) if top_unique_states else []
    bottom_descendants = sorted(bottom_unique_states) if bottom_unique_states else []

    # Determine top_change_type and bottom_change_type for later use
    if top_codivergent:
        top_change_type = "codivergent"
    elif top_convergent:
        top_change_type = "convergent"
    elif top_divergent:
        top_change_type = "divergent"
    elif top_insufficient:
        top_change_type = "no_change"  # Insufficient evidence
    else:
        top_change_type = "none"

    if bottom_codivergent:
        bottom_change_type = "codivergent"
    elif bottom_convergent:
        bottom_change_type = "convergent"
    elif bottom_divergent:
        bottom_change_type = "divergent"
    elif bottom_insufficient:
        bottom_change_type = "no_change"  # Insufficient evidence
    else:
        bottom_change_type = "none"

    pattern = "ambiguous"
    description = "Complex pattern requiring manual inspection"

    # Pattern classification with new hierarchy
    if not top_changed and not bottom_changed:
        pattern = "no_change"
        description = "No substitutions detected in any group"

    elif top_insufficient and bottom_insufficient:
        pattern = "no_change"
        description = "Insufficient changes (<2 per group) to assess pattern"

    elif top_insufficient and bottom_insufficient:
        pattern = "no_change"
        description = "Insufficient changes (<2 per group) to assess pattern"

    elif (top_sufficient_changes and bottom_insufficient) or (
        top_insufficient and bottom_sufficient_changes
    ):
        # Only one side has sufficient changes
        if top_convergent and bottom_insufficient:
            top_state = top_descendants[0] if top_descendants else "?"
            pattern = "convergent"
            description = f"Convergent (TOP only): → {top_state}, BOTTOM insufficient changes"
        elif bottom_convergent and top_insufficient:
            bottom_state = bottom_descendants[0] if bottom_descendants else "?"
            pattern = "convergent"
            description = f"Convergent (BOTTOM only): → {bottom_state}, TOP insufficient changes"
        elif top_divergent and bottom_insufficient:
            all_descendants = sorted(top_descendants)
            pattern = "divergent"
            description = (
                f"Divergent (TOP only): {'/'.join(all_descendants)}, BOTTOM insufficient changes"
            )
        elif bottom_divergent and top_insufficient:
            all_descendants = sorted(bottom_descendants)
            pattern = "divergent"
            description = (
                f"Divergent (BOTTOM only): {'/'.join(all_descendants)}, TOP insufficient changes"
            )
        elif top_codivergent and bottom_insufficient:
            pattern = "codivergent_top"
            description = (
                f"Codivergent TOP: {'/'.join(top_descendants)}, BOTTOM insufficient changes"
            )
        elif bottom_codivergent and top_insufficient:
            pattern = "codivergent_bottom"
            description = (
                f"Codivergent BOTTOM: {'/'.join(bottom_descendants)}, TOP insufficient changes"
            )

    elif (top_convergent or top_divergent or top_codivergent) and bottom_unchanged:
        # TOP changed, BOTTOM unchanged
        if top_convergent:
            top_state = top_descendants[0] if top_descendants else "?"
            pattern = "convergent"
            description = f"Convergent (TOP only): → {top_state}, BOTTOM unchanged"
        elif top_divergent:
            all_descendants = sorted(top_descendants)
            pattern = "divergent"
            description = f"Divergent (TOP only): {'/'.join(all_descendants)}, BOTTOM unchanged"
        elif top_codivergent:
            pattern = "codivergent_top"
            description = f"Codivergent TOP: {'/'.join(top_descendants)}, BOTTOM unchanged"

    elif (bottom_convergent or bottom_divergent or bottom_codivergent) and top_unchanged:
        # BOTTOM changed, TOP unchanged
        if bottom_convergent:
            bottom_state = bottom_descendants[0] if bottom_descendants else "?"
            pattern = "convergent"
            description = f"Convergent (BOTTOM only): → {bottom_state}, TOP unchanged"
        elif bottom_divergent:
            all_descendants = sorted(bottom_descendants)
            pattern = "divergent"
            description = f"Divergent (BOTTOM only): {'/'.join(all_descendants)}, TOP unchanged"
        elif bottom_codivergent:
            pattern = "codivergent_bottom"
            description = f"Codivergent BOTTOM: {'/'.join(bottom_descendants)}, TOP unchanged"

    elif top_sufficient_changes and bottom_sufficient_changes:
        # Both sides have sufficient changes - determine parallel subtype
        if top_convergent and bottom_convergent:
            if set(top_descendants) == set(bottom_descendants):
                pattern = "no_convergence"
                description = (
                    f"Both groups converge to the same state: TOP/BOTTOM → {top_descendants[0]}"
                )
            else:
                pattern = "parallel_convergence"
                description = (
                    f"Parallel convergence: TOP → {top_descendants[0]}, "
                    f"BOTTOM → {bottom_descendants[0]}"
                )
        elif top_divergent and bottom_divergent:
            pattern = "parallel_divergence"
            top_states = "/".join(top_descendants)
            bottom_states = "/".join(bottom_descendants)
            description = f"Parallel divergence: TOP → {top_states}, BOTTOM → {bottom_states}"
        elif top_codivergent and bottom_codivergent:
            pattern = "parallel_codivergent"
            description = f"Parallel codivergence: both sides show convergence+divergence"
        elif (top_convergent or top_codivergent) and (bottom_divergent or bottom_codivergent):
            # One side converges/codiverges, other diverges/codiverges = mixed
            pattern = "parallel_mixed"
            description = f"Parallel mixed: TOP {top_change_type}, BOTTOM {bottom_change_type}"
        elif (bottom_convergent or bottom_codivergent) and (top_divergent or top_codivergent):
            pattern = "parallel_mixed"
            description = f"Parallel mixed: TOP {top_change_type}, BOTTOM {bottom_change_type}"
        else:
            # Fallback for unexpected combinations
            pattern = "parallel_mixed"
            description = f"Parallel mixed: TOP {top_change_type}, BOTTOM {bottom_change_type}"

    result = {
        "pattern": pattern,
        "description": description,
        "top_convergent": top_convergent,
        "bottom_convergent": bottom_convergent,
        "top_divergent": top_divergent,
        "bottom_divergent": bottom_divergent,
        "top_codivergent": top_codivergent,
        "bottom_codivergent": bottom_codivergent,
        "top_descendants": top_descendants,
        "bottom_descendants": bottom_descendants,
        "top_changed": top_descendants,
        "bottom_changed": bottom_descendants,
        "ancestors": sorted(ancestor_states),
        "descendants": sorted(set(top_descendants + bottom_descendants)),
        "top_change_type": top_change_type,
        "bottom_change_type": bottom_change_type,
        "top_change_count": top_change_count,
        "bottom_change_count": bottom_change_count,
        # Add diagnostic info
        "diagnostic": {
            "top_transitions": top_transitions,
            "bottom_transitions": bottom_transitions,
            "top_unchanged": top_unchanged,
            "bottom_unchanged": bottom_unchanged,
            "top_sufficient_changes": top_sufficient_changes,
            "bottom_sufficient_changes": bottom_sufficient_changes,
        },
    }

    return result


def describe_transition(trans: Dict[str, Any]) -> str:
    """
    Format a transition dict as human-readable string.

    Args:
        trans: Transition dict with 'ancestor', 'descendant', 'status' keys

    Returns:
        Formatted string (e.g., "A→V (changed)")

    Example:
        >>> trans = {'ancestor': 'A', 'descendant': 'V', 'status': 'changed'}
        >>> describe_transition(trans)
        'A→V (changed)'
    """
    ancestor = trans.get("ancestor") or "?"
    descendant = trans.get("descendant") or "?"
    status = trans.get("status")
    return f"{ancestor}→{descendant} ({status})"


def build_consolidated_multiset(
    alignment,
    msa_pos: int,
    trait_pairs: List[Tuple[str, str]],
    hypothesis_caas: str,
    seq_by_id: Optional[Dict[str, str]] = None,
    seq_by_species: Optional[Dict[str, str]] = None,
    taxid_to_species: Optional[Dict[str, str]] = None,
) -> str:
    """
    Build consolidated amino acid multiset from actual tip states.

    Extracts residues from ALL species in trait file (involved + uninvolved pairs),
    builds multiset with one letter per pair (grouped by side), and ensures
    hypothesis residues are present.

    Args:
        alignment: BioPython alignment object
        msa_pos: MSA position (0-based)
        trait_pairs: List of (species1, species2) tuples from trait file
        hypothesis_caas: Expected pattern from CAAS discovery (e.g., "A/V")
        seq_by_id: Optional prebuilt sequence lookup by taxid
        seq_by_species: Optional prebuilt sequence lookup by species name
        taxid_to_species: Optional mapping from taxid to species name

    Returns:
        new_convAA string with multiplicities (e.g., "AAV/VVV" for 3 pairs)

    Example:
        >>> # 3 pairs: (Sp1,Sp2), (Sp3,Sp4), (Sp5,Sp6)
        >>> # Discovery pattern: A/V
        >>> # Observed: Sp1=A, Sp2=V; Sp3=A, Sp4=V; Sp5=A, Sp6=I
        >>> build_consolidated_multiset(
        ...     alignment, 85,
        ...     [('Sp1','Sp2'), ('Sp3','Sp4'), ('Sp5','Sp6')],
        ...     'A/V'
        ... )
        'AAA/VVI'  # Left side: A,A,A; Right side: V,V,I
    """
    logger.info(f"=== build_consolidated_multiset START ===")
    logger.info(f"Gene: {getattr(alignment, 'gene_name', 'unknown')}")
    logger.info(f"MSA position (0-based): {msa_pos}")
    logger.info(f"Number of trait pairs: {len(trait_pairs)}")
    logger.debug(f"Trait pairs: {trait_pairs}")

    # Build lookups if not provided
    if seq_by_id is None or seq_by_species is None:
        logger.debug("Building alignment lookups...")
        seq_by_id, seq_by_species = build_alignment_lookup(alignment, taxid_to_species)
        logger.debug(
            f"Alignment contains {len(seq_by_id)} sequences by ID, {len(seq_by_species)} by species"
        )
        logger.debug(f"Sample seq_by_id keys: {list(seq_by_id.keys())[:5]}")
        logger.debug(f"Sample seq_by_species keys: {list(seq_by_species.keys())[:5]}")
        logger.debug(f"taxid_to_species mapping: {len(taxid_to_species or {})} entries")

    taxid_to_species = taxid_to_species or {}

    # Parse hypothesis pattern
    if "/" in hypothesis_caas:
        hyp_left, hyp_right = hypothesis_caas.split("/", 1)
        logger.debug(f"Parsed hypothesis: LEFT='{hyp_left}', RIGHT='{hyp_right}'")
    else:
        logger.warning(
            f"Invalid caas pattern '{hypothesis_caas}', expected 'AA/BB' format"
        )
        hyp_left, hyp_right = hypothesis_caas, hypothesis_caas

    left_residues = []
    right_residues = []

    # Collect residues from all pairs
    logger.info(f"Collecting residues from {len(trait_pairs)} pairs...")
    skipped_pairs = []
    for pair_idx, (species1, species2) in enumerate(trait_pairs, 1):
        logger.debug(f"  Pair {pair_idx}: ({species1}, {species2})")

        # Fetch residue for left species (first in pair)
        res1 = _fetch_residue_from_lookup(species1, msa_pos, seq_by_id, seq_by_species)
        logger.debug(f"    LEFT species '{species1}' -> residue: {res1}")

        if res1 is None:
            # Try with underscore replacement (common species name format)
            species1_alt = species1.replace(" ", "_")
            res1 = _fetch_residue_from_lookup(species1_alt, msa_pos, seq_by_id, seq_by_species)
            logger.debug(f"    LEFT retry with '{species1_alt}' -> residue: {res1}")

        if res1 is None:
            logger.warning(
                f"Skipping pair {pair_idx}: LEFT species '{species1}' not found in alignment at position {msa_pos}. "
                f"This pair will not contribute to the multiset."
            )
            skipped_pairs.append((species1, species2, "LEFT species not in alignment"))
            continue

        if res1 == "-":
            logger.warning(
                f"Skipping pair {pair_idx}: LEFT species '{species1}' has gap at position {msa_pos}. "
                f"This pair will not contribute to the multiset."
            )
            skipped_pairs.append((species1, species2, "LEFT species has gap"))
            continue

        # Fetch residue for right species (second in pair)
        res2 = _fetch_residue_from_lookup(species2, msa_pos, seq_by_id, seq_by_species)
        logger.debug(f"    RIGHT species '{species2}' -> residue: {res2}")

        if res2 is None:
            species2_alt = species2.replace(" ", "_")
            res2 = _fetch_residue_from_lookup(species2_alt, msa_pos, seq_by_id, seq_by_species)
            logger.debug(f"    RIGHT retry with '{species2_alt}' -> residue: {res2}")

        if res2 is None:
            logger.warning(
                f"Skipping pair {pair_idx}: RIGHT species '{species2}' not found in alignment at position {msa_pos}. "
                f"This pair will not contribute to the multiset."
            )
            skipped_pairs.append((species1, species2, "RIGHT species not in alignment"))
            continue

        if res2 == "-":
            logger.warning(
                f"Skipping pair {pair_idx}: RIGHT species '{species2}' has gap at position {msa_pos}. "
                f"This pair will not contribute to the multiset."
            )
            skipped_pairs.append((species1, species2, "RIGHT species has gap"))
            continue

        logger.info(f"    Pair {pair_idx} ({species1} vs {species2}): LEFT={res1}, RIGHT={res2}")
        left_residues.append(res1)
        right_residues.append(res2)

    if skipped_pairs:
        logger.warning(
            f"Skipped {len(skipped_pairs)} pairs due to missing species or gaps: {skipped_pairs}"
        )

    logger.info(f"Collected LEFT residues: {left_residues}")
    logger.info(f"Collected RIGHT residues: {right_residues}")

    # Build multiset strings from ONLY observed residues (no injection)
    # The multiset represents actual tip-level states in trait pair order
    left_multiset = "".join(left_residues)
    right_multiset = "".join(right_residues)

    # Keep ordering exactly as pairs appear in trait file (no reordering)
    new_convAA = f"{left_multiset}/{right_multiset}"

    logger.info(f"Final multiset: {new_convAA}")
    logger.info(f"  LEFT multiset (sorted): {left_multiset} (from {left_residues})")
    logger.info(f"  RIGHT multiset (sorted): {right_multiset} (from {right_residues})")
    logger.info(f"  Hypothesis pattern: {hypothesis_caas}")
    logger.info(f"=== build_consolidated_multiset END ===")

    return new_convAA
