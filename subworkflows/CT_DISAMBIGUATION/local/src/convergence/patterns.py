"""
Convergence Pattern Classification
==================================

Pure helpers for transition summarization and pattern classification.
Operates on in-memory structures with no hidden I/O. Delegates tip-level classification
to the core convergence module while providing transition-level analysis utilities.

Key Functions
-------------
**transition_status**: Classify ancestor→descendant transition as maintained/changed/unknown
**summarize_pair_transitions**: Build per-pair transition summaries with ancestor/descendant states
**classify_focus_transitions**: Label focused side as convergent/divergent/ambiguous

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

__all__ = [
    "transition_status",
    "summarize_pair_transitions",
    "classify_focus_transitions",
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
    descendants = {
        t.get("descendant") for t in focus_transitions if t.get("descendant")
    }

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
