"""
Convergence Pattern Classification
==================================

Pure helpers for transition summarization, pattern classification, and stability assessment.
Operates on in-memory structures with no hidden I/O. Delegates tip-level classification
to the core convergence module while providing transition-level analysis utilities.

Key Functions
-------------
**transition_status**: Classify ancestor→descendant transition as maintained/changed/unknown
**summarize_pair_transitions**: Build per-pair transition summaries with ancestor/descendant states
**classify_focus_transitions**: Label focused side as convergent/divergent/ambiguous
**assess_convergence_stability**: Determine stability based on within-side conservation

Stability Criteria
------------------
A convergence is **STABLE** if at least one side (top or bottom) is fully conserved
across all pairs. Conservation means a single unique residue for that side.

Examples:
    - 'AAA/TSS' → STABLE (top conserved as A)
    - 'AAS/TTT' → STABLE (bottom conserved as T)
    - 'AAA/TTT' → STABLE (both sides conserved)
    - 'AAS/TST' → AMBIGUOUS (neither side fully conserved)

Data Contracts
--------------
- **Transition**: TypedDict for ancestor→descendant state with metrics
- **TransitionSummary**: TypedDict for per-pair transition pairs (top/bottom)

Usage Example
-------------
::

    from src.convergence.patterns import summarize_pair_transitions, classify_focus_transitions

    pair_details = [
        {'pair_id': 'pair_1', 'focal_state': 'A',
         'top_tip_mode': 'V', 'bottom_tip_mode': 'T'},
        {'pair_id': 'pair_2', 'focal_state': 'A',
         'top_tip_mode': 'V', 'bottom_tip_mode': 'I'}
    ]
    summary = summarize_pair_transitions(pair_details)
    label, detail = classify_focus_transitions(summary, focus='top')
    print(label)  # 'convergent'
    print(detail)  # 'Ancestors=['A'], Descendants=['V']'

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-07
"""

from typing import Any, Dict, List, Optional, Tuple, TypedDict

from .convergence import classify_tip_level_pattern

__all__ = [
    "transition_status",
    "summarize_pair_transitions",
    "classify_focus_transitions",
    "classify_tip_level_pattern",
    "assess_convergence_stability",
]


class Transition(TypedDict, total=False):
    """Transition descriptor between ancestor and descendant states."""

    ancestor: Optional[str]
    descendant: Optional[str]
    status: str


class TransitionSummary(TypedDict, total=False):
    """Summary per pair with top/bottom transitions."""

    pair_id: Optional[str]
    transitions: Dict[str, Transition]


def transition_status(ancestor: Optional[str], descendant: Optional[str]) -> str:
    """
    Classify transition between ancestral and descendant states.

    Args:
        ancestor: Ancestral amino acid state
        descendant: Descendant amino acid state

    Returns:
        Status string: 'maintained', 'changed', or 'unknown'
    """
    if not ancestor or not descendant:
        return "unknown"
    if ancestor == descendant:
        return "maintained"
    return "changed"


def summarize_pair_transitions(
    pair_details: Optional[List[Dict[str, Any]]],
) -> List[TransitionSummary]:
    """
    Summarize transitions for each pair (top and bottom).

    Args:
        pair_details: List of pair detail dictionaries with focal_state and tip states

    Returns:
        List of transition summaries with pair_id and transitions dict
    """
    summary: List[TransitionSummary] = []
    for pair in pair_details or []:
        ancestor = pair.get("focal_state") or pair.get("mrca_modal_aa")
        transitions = {}
        for side in ("top", "bottom"):
            tip_state = pair.get(f"{side}_tip_mode") or pair.get(f"{side}_tip_residue")
            transitions[side] = {
                "ancestor": ancestor,
                "descendant": tip_state,
                "status": transition_status(ancestor, tip_state),
            }
        summary.append({"pair_id": pair.get("pair_id"), "transitions": transitions})
    return summary


def classify_focus_transitions(
    transition_summary: List[TransitionSummary], focus: str = "top"
) -> Tuple[str, str]:
    """
    Classify transitions for a focused side (top or bottom) as convergent/divergent.

    Args:
        transition_summary: List of transition summaries from summarize_pair_transitions
        focus: Which side to focus on ('top' or 'bottom')

    Returns:
        Tuple of (classification label, detail description)
        - Label: 'convergent', 'divergent', or 'ambiguous'
        - Detail: Human-readable description
    """
    focus_transitions: List[Transition] = []
    for item in transition_summary:
        transitions = item.get("transitions") or {}
        trans_focus = transitions.get(focus)
        if not trans_focus:
            continue
        if trans_focus.get("status") != "unknown":
            focus_transitions.append(trans_focus)

    if len(focus_transitions) < 2:
        return "ambiguous", "Not enough resolved transitions"

    ancestors = {t.get("ancestor") for t in focus_transitions if t.get("ancestor")}
    descendants = {t.get("descendant") for t in focus_transitions if t.get("descendant")}

    ancestor_list = sorted([a for a in ancestors if a is not None])
    descendant_list = sorted([d for d in descendants if d is not None])

    if not descendants:
        return "ambiguous", "Missing descendant states"

    if len(descendants) == 1:
        label = "convergent"
    else:
        label = "divergent"

    detail = f"Ancestors={ancestor_list if ancestor_list else ['?']}, Descendants={descendant_list}"
    return label, detail


def assess_convergence_stability(multi_caas: str) -> Dict[str, Any]:
    """
    Assess convergence stability based on conservation pattern in multi_caas string.

    **CONVERGENCE STABILITY CRITERIA**:
    A convergence is considered STABLE if either the top or bottom side
    (or both) is fully conserved across all pairs.

    Examples:
        - "AAA/TSS" → STABLE (top fully conserved as A)
        - "AAS/TTT" → STABLE (bottom fully conserved as T)
        - "AAA/TTT" → STABLE (both sides fully conserved)
        - "AAS/TST" → AMBIGUOUS (neither side fully conserved)
        - "ABC/XYZ" → AMBIGUOUS (high variation on both sides)

    This is independent of antiCAAS flags. AntiCAAS tracks pattern contradictions
    in uninvolved pairs, while stability tracks within-side conservation.

    Args:
        multi_caas: Multi-species CAAS pattern string (e.g., "AAA/TSS" or "AAS/TTT")

    Returns:
        Dictionary with:
            - 'is_stable' (bool): True if convergence is stable
            - 'stability_pattern' (str): 'top_conserved', 'bottom_conserved', 'both_conserved', or 'ambiguous'
            - 'top_residues' (set): Unique residues on top side
            - 'bottom_residues' (set): Unique residues on bottom side
            - 'top_conserved' (bool): True if top side fully conserved
            - 'bottom_conserved' (bool): True if bottom side fully conserved

    Example:
        >>> assess_convergence_stability("AAA/TSS")
        {
            'is_stable': True,
            'stability_pattern': 'top_conserved',
            'top_residues': {'A'},
            'bottom_residues': {'S', 'T'},
            'top_conserved': True,
            'bottom_conserved': False
        }
    """
    if not multi_caas or "/" not in multi_caas:
        return {
            "is_stable": False,
            "stability_pattern": "invalid",
            "top_residues": set(),
            "bottom_residues": set(),
            "top_conserved": False,
            "bottom_conserved": False,
        }

    top_str, bottom_str = multi_caas.split("/", 1)

    # Extract unique residues (excluding gaps, X, ?)
    top_residues = {r for r in top_str if r not in {"-", "X", "?", ""}}
    bottom_residues = {r for r in bottom_str if r not in {"-", "X", "?", ""}}

    # Check conservation (single unique residue = fully conserved)
    top_conserved = len(top_residues) == 1
    bottom_conserved = len(bottom_residues) == 1

    # Determine stability
    if top_conserved and bottom_conserved:
        is_stable = True
        stability_pattern = "both_conserved"
    elif top_conserved:
        is_stable = True
        stability_pattern = "top_conserved"
    elif bottom_conserved:
        is_stable = True
        stability_pattern = "bottom_conserved"
    else:
        is_stable = False
        stability_pattern = "ambiguous"

    return {
        "is_stable": is_stable,
        "stability_pattern": stability_pattern,
        "top_residues": top_residues,
        "bottom_residues": bottom_residues,
        "top_conserved": top_conserved,
        "bottom_conserved": bottom_conserved,
    }
