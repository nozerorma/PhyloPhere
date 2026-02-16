#!/usr/bin/env python3
"""
Amino Acid State Inference and Change Classification
=====================================================

Functions for inferring ancestral and derived amino acid states from ASR
(Ancestral State Reconstruction) node states and evolutionary pattern
classifications. Determines which phylogenetic sides (TOP/BOTTOM trait groups)
experienced changes based on pattern type and change classifications.

This module handles pattern-aware logic for:

- **Convergent/parallel patterns**: Both sides change toward same derived state
- **Divergent/codivergent patterns**: Sides change in opposite directions
- **Asymmetric codivergent**: One side converges while other diverges

Core Functions
--------------
:func:`compute_change_side`: Determine which sides experienced changes

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-06
"""

from typing import Optional
import logging

logger = logging.getLogger(__name__)

#: Accepted change types for TOP/BOTTOM groups
VALID_CHANGE_TYPES = {"none", "convergent", "divergent", "changed", "conserved"}


def compute_change_side(
    top_change_type: str, bottom_change_type: str, pattern_type: Optional[str] = None
) -> str:
    """
    Determine which side(s) experienced amino acid changes.

    Analyzes change types for TOP and BOTTOM trait groups to determine which
    phylogenetic sides experienced changes, with pattern-aware logic for
    different evolutionary scenarios.

    :param top_change_type: Change type in top group
                            ('none', 'convergent', 'divergent', 'changed')
    :type top_change_type: str
    :param bottom_change_type: Change type in bottom group
                               ('none', 'convergent', 'divergent', 'changed')
    :type bottom_change_type: str
    :param pattern_type: Overall pattern classification (optional)
                         ('convergent', 'divergent', 'parallel',
                         'codivergent', 'codivergent_top_converges', etc.)
    :type pattern_type: Optional[str]

    :returns: String indicating change side - 'both', 'top', 'bottom', or 'none'
    :rtype: str

    :raises ValueError: If pattern_type is provided but not recognized

    Note
    ----
    Pattern-specific logic:

    - **Convergent/parallel**: Both sides typically change toward same state
    - **Divergent/codivergent**: Sides change in opposite directions
    - **Codivergent with asymmetric convergence**: One side converges while
      other diverges

    Example
    -------
    ::

        >>> compute_change_side('changed', 'changed', 'convergent')
        'both'
        >>> compute_change_side('changed', 'none', 'parallel')
        'top'
        >>> compute_change_side('none', 'none')
        'none'
    """
    # Validate change types and log unexpected values without failing hard
    if top_change_type not in VALID_CHANGE_TYPES:
        logger.warning(
            "Unexpected top_change_type '%s'. Valid values: %s",
            top_change_type,
            VALID_CHANGE_TYPES,
        )
    if bottom_change_type not in VALID_CHANGE_TYPES:
        logger.warning(
            "Unexpected bottom_change_type '%s'. Valid values: %s",
            bottom_change_type,
            VALID_CHANGE_TYPES,
        )

    logger.debug(
        "Resolving change side with top=%s, bottom=%s, pattern=%s",
        top_change_type,
        bottom_change_type,
        pattern_type,
    )

    # No pattern provided: determine based on change types alone
    if pattern_type is None:
        if top_change_type not in ("none", "insufficient") and bottom_change_type not in ("none", "insufficient"):
            return "both"
        elif top_change_type not in ("none", "insufficient"):
            return "top"
        elif bottom_change_type not in ("none", "insufficient"):
            return "bottom"
        else:
            return "none"

    # Codivergent patterns (one side only)
    if pattern_type in ("codivergent_top", "codivergent_bottom"):
        return "top" if pattern_type == "codivergent_top" else "bottom"

    # No change pattern
    if pattern_type == "no_change":
        return "none"

    # Divergent patterns: changes occur in one or both directions
    if pattern_type == "divergent":
        if top_change_type != "none" and bottom_change_type != "none":
            return "both"
        elif top_change_type != "none":
            return "top"
        elif bottom_change_type != "none":
            return "bottom"
        else:
            return "none"

    # Convergent/parallel patterns: changes toward same or different derived states
    if pattern_type in (
        "convergent",
        "parallel_convergence",
        "parallel_divergence",
        "parallel_mixed",
        "parallel_codivergent",
    ):
        if top_change_type != "none" and bottom_change_type != "none":
            return "both"
        elif top_change_type != "none":
            return "top"
        elif bottom_change_type != "none":
            return "bottom"
        else:
            return "none"

    # Edge cases: ambiguous, insufficient data, no convergence, unknown
    # Treat as neutral/no-change patterns to avoid breaking downstream analysis
    if pattern_type in ("ambiguous", "insufficient_data", "no_convergence", "unknown"):
        if top_change_type != "none" and bottom_change_type != "none":
            return "both"
        elif top_change_type != "none":
            return "top"
        elif bottom_change_type != "none":
            return "bottom"
        else:
            return "none"

    # Unrecognized pattern type
    raise ValueError(f"Unrecognized pattern_type: {pattern_type}")
