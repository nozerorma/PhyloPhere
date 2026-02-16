"""Convergence exports for ASR-driven pattern classification."""

from .convergence import (
    NodeStates,
    ConvergenceClassification,
    extract_node_states_from_node_level,
    build_alignment_lookup,
    collect_tip_residues,
    normalize_amino_list,
    format_amino_display,
    compute_derived_state_similarity,
    classify_tip_level_pattern,
    describe_transition,
)

from .patterns import (
    transition_status,
    summarize_pair_transitions,
    classify_focus_transitions,
)
from .node_mapping import (
    build_convergence_node_mapping,
)

__all__ = [
    # Pattern classification
    "classify_tip_level_pattern",
    "transition_status",
    "summarize_pair_transitions",
    "classify_focus_transitions",
    # Node mapping
    "build_convergence_node_mapping",
    # Convergence
    "NodeStates",
    "ConvergenceClassification",
    "extract_node_states_from_node_level",
    # Tip-level analysis
    "build_alignment_lookup",
    "collect_tip_residues",
    "normalize_amino_list",
    "format_amino_display",
    "compute_derived_state_similarity",
    "describe_transition",
]
