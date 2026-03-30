"""Public exports for convergence classification and disambiguation."""

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
from .disambiguate_single import (
    analyze_caas_position_disambiguation,
    analyze_gene_disambiguation,
)

__all__ = [
    "ConvergenceClassification",
    "NodeStates",
    "analyze_caas_position_disambiguation",
    "analyze_gene_disambiguation",
    "build_alignment_lookup",
    "build_convergence_node_mapping",
    "classify_focus_transitions",
    "classify_tip_level_pattern",
    "collect_tip_residues",
    "compute_derived_state_similarity",
    "describe_transition",
    "extract_node_states_from_node_level",
    "format_amino_display",
    "normalize_amino_list",
    "summarize_pair_transitions",
    "transition_status",
]
