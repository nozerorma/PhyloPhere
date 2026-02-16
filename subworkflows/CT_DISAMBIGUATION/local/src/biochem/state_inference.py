#!/usr/bin/env python3
"""Infer which trait side(s) changed at a CAAS position.

This module centralizes side-level change inference (``top``, ``bottom``,
``both``, or ``none``) from per-side change labels and optional
pattern-level classifications.
"""

from typing import Optional
import logging

logger = logging.getLogger(__name__)

VALID_CHANGE_TYPES = {
    "none",
    "convergent",
    "divergent",
    "changed",
    "conserved",
    "insufficient",
}


def _has_change(change_type: str) -> bool:
    """Return ``True`` when a change label represents a non-neutral change."""
    return change_type not in ("none", "insufficient")


def compute_change_side(
    top_change_type: str, bottom_change_type: str, pattern_type: Optional[str] = None
) -> str:
    """
    Determine which side(s) experienced amino-acid changes.

    :param top_change_type: Change type label for the TOP trait group.
    :type top_change_type: str
    :param bottom_change_type: Change type label for the BOTTOM trait group.
    :type bottom_change_type: str
    :param pattern_type: Optional global pattern class for the position.
    :type pattern_type: Optional[str]

    :returns: ``"both"``, ``"top"``, ``"bottom"``, or ``"none"``.
    :rtype: str
    :raises ValueError: If ``pattern_type`` is not recognized.
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

    # No pattern provided: infer directly from side-level change labels.
    if pattern_type is None:
        if _has_change(top_change_type) and _has_change(bottom_change_type):
            return "both"
        if _has_change(top_change_type):
            return "top"
        if _has_change(bottom_change_type):
            return "bottom"
        return "none"

    # Codivergent patterns (one side only).
    if pattern_type in ("codivergent_top", "codivergent_bottom"):
        return "top" if pattern_type == "codivergent_top" else "bottom"

    # No-change pattern.
    if pattern_type == "no_change":
        return "none"

    # Divergent patterns: changes can occur on one side or both.
    if pattern_type == "divergent":
        if _has_change(top_change_type) and _has_change(bottom_change_type):
            return "both"
        if _has_change(top_change_type):
            return "top"
        if _has_change(bottom_change_type):
            return "bottom"
        return "none"

    # Convergent/parallel patterns.
    if pattern_type in (
        "convergent",
        "parallel_convergence",
        "parallel_divergence",
        "parallel_mixed",
        "parallel_codivergent",
    ):
        if _has_change(top_change_type) and _has_change(bottom_change_type):
            return "both"
        if _has_change(top_change_type):
            return "top"
        if _has_change(bottom_change_type):
            return "bottom"
        return "none"

    # Neutral/edge pattern classes used in downstream summaries.
    if pattern_type in ("ambiguous", "insufficient_data", "no_convergence", "unknown"):
        if _has_change(top_change_type) and _has_change(bottom_change_type):
            return "both"
        if _has_change(top_change_type):
            return "top"
        if _has_change(bottom_change_type):
            return "bottom"
        return "none"

    # Unrecognized pattern type
    raise ValueError(f"Unrecognized pattern_type: {pattern_type}")
