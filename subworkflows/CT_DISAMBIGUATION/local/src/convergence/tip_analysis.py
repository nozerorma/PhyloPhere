"""
Tip-Level Convergence Analysis
===============================

Collects tip residues for contrasts and classifies tip-level convergence patterns.
Pure helpers expecting alignment data, ASR posteriors, and MRCA lookups from callers.

Workflow
--------
1. **collect_contrast_tip_residues**: Extract observed residues at tip species for a contrast
2. **analyze_tip_convergence_patterns**: Classify pattern and populate diagnostics

Diagnostics Structure
---------------------
The diagnostics dictionary populated by `analyze_tip_convergence_patterns` includes:
    - pair_details: Enriched pair detail dicts with tip/focal states
    - pair_transition_summary: Per-pair ancestor→descendant transitions
    - transition_focus_side: Focus side for classification ('top' or 'bottom')
    - focus_transition_detail: Human-readable transition summary

Usage Example
-------------
::

    from src.convergence.tip_analysis import collect_contrast_tip_residues, analyze_tip_convergence_patterns

    # Collect tip residues
    enriched = collect_contrast_tip_residues(
        contrast, position_info, alignment_data,
        seq_by_id, seq_by_species, mrca_func, node_posteriors
    )

    # Classify pattern
    diagnostics = {}
    pattern = analyze_tip_convergence_patterns([enriched], diagnostics)
    print(pattern['pattern'])  # 'convergent' or 'divergent'
    print(diagnostics['focus_transition_detail'])  # 'Ancestors=['A'], Descendants=['V']'

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-07
"""

from typing import Any, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


def collect_contrast_tip_residues(
    contrast: Any,
    position_info: Any,
    alignment_data: Any,
    seq_by_id: Dict[str, str],
    seq_by_species: Dict[str, str],
    mrca_func: Any,
    node_posteriors: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
) -> Optional[Any]:
    """
    Collect tip residues for a single contrast at a position.

    :param contrast: ContrastDefinition-like object (expects top_taxa, bottom_taxa,
                     top_species, bottom_species, all_taxa, pair_id fields)
    :param position_info: Object with ``position_zero_based`` and ``position_one_based``
    :param alignment_data: Alignment metadata (expects ``taxid_to_species`` map)
    :param seq_by_id: Alignment lookup keyed by taxid/record id
    :param seq_by_species: Alignment lookup keyed by species name
    :param mrca_func: Callable returning MRCA node with ``node_id`` for a taxid list
    :param node_posteriors: Optional ASR posterior probabilities ``{node: {site: {aa: prob}}}``

    :returns: Enriched contrast with tip residues and focal ASR state, or ``None``
              when MRCA or ASR data are unavailable
    :rtype: Optional[Any]
    """
    from src.convergence.convergence import collect_tip_residues, extract_tip_residue

    # Find MRCA for this pair
    mrca_node = mrca_func(contrast.all_taxa)
    if not mrca_node or mrca_node.node_id is None:
        logger.debug(f"No MRCA found for pair {contrast.pair_id}")
        return None

    # Build enriched contrast
    contrast_dict = contrast.__dict__.copy()
    if "node_id" in contrast_dict:
        contrast_dict.pop("node_id")

    # Import ContrastDefinition from the data models location
    from dataclasses import replace

    enriched_contrast = replace(contrast, node_id=mrca_node.node_id)

    position_zero_based = position_info.position_zero_based

    # Collect tip residues for this position
    top_tip_records = collect_tip_residues(
        enriched_contrast.top_taxa,
        enriched_contrast.top_species,
        position_zero_based,
        seq_by_id,
        seq_by_species,
        alignment_data.taxid_to_species,
    )

    bottom_tip_records = collect_tip_residues(
        enriched_contrast.bottom_taxa,
        enriched_contrast.bottom_species,
        position_zero_based,
        seq_by_id,
        seq_by_species,
        alignment_data.taxid_to_species,
    )

    # Get modal residues
    enriched_contrast.top_tip_residues = top_tip_records
    enriched_contrast.bottom_tip_residues = bottom_tip_records
    enriched_contrast.top_tip_residue = extract_tip_residue(top_tip_records)
    enriched_contrast.bottom_tip_residue = extract_tip_residue(bottom_tip_records)
    enriched_contrast.top_tip_mode = enriched_contrast.top_tip_residue
    enriched_contrast.bottom_tip_mode = enriched_contrast.bottom_tip_residue

    # Get focal state from ASR if available
    def get_modal_state(node_id: Optional[int], paml_site: int) -> Optional[str]:
        if not node_posteriors or node_id is None:
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

    paml_site = position_info.position_one_based
    enriched_contrast.focal_state = get_modal_state(mrca_node.node_id, paml_site)

    # Return None if ASR not available - position cannot be assessed
    if enriched_contrast.focal_state is None:
        logger.debug(
            f"No ASR available for pair {contrast.pair_id} at position "
            f"{position_info.position_one_based}, skipping position"
        )
        return None

    enriched_contrast.mrca_modal_aa = enriched_contrast.focal_state

    return enriched_contrast


def analyze_tip_convergence_patterns(
    pair_details: List[Dict[str, Any]], diagnostics: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """
    Classify and annotate tip-level convergence pattern.

    :param pair_details: Enriched pair detail dictionaries (expects tip/focal states)
    :param diagnostics: Diagnostics dictionary to populate in-place
    :returns: Pattern dict or ``None`` if insufficient data
    :rtype: Optional[Dict[str, Any]]
    """
    from .patterns import (
        summarize_pair_transitions,
        classify_focus_transitions,
    )
    from .convergence import classify_change_and_parallelism

    if len(pair_details) < 2:
        return None

    # Cast to Sequence for type checking (dicts are compatible at runtime)
    tip_level_pattern = classify_change_and_parallelism(pair_details)  # type: ignore[arg-type]
    logger.info(f"✓ Tip-level pattern: {tip_level_pattern.get('pattern_type', 'unknown')}")

    diagnostics["pair_details"] = pair_details
    diagnostics["pair_transition_summary"] = summarize_pair_transitions(pair_details)

    # Determine focus side
    def _has_changes(side: str) -> bool:
        for item in diagnostics["pair_transition_summary"]:
            if (item.get("transitions") or {}).get(side, {}).get("status") == "changed":
                return True
        return False

    focus_side = "top"
    if not _has_changes("top") and _has_changes("bottom"):
        focus_side = "bottom"
    diagnostics["transition_focus_side"] = focus_side

    classification, detail = classify_focus_transitions(
        diagnostics["pair_transition_summary"], focus=focus_side
    )
    diagnostics["focus_transition_detail"] = detail

    return tip_level_pattern
