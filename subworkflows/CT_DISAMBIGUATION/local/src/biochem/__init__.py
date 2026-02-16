"""Biochem module: amino-acid grouping + change-side inference only."""

from .grouping import GS1, GS2, GS3, GS4, get_grouping_scheme

from .state_inference import (
    compute_change_side,
)

__all__ = [
    # Grouping schemes
    "GS1",
    "GS2",
    "GS3",
    "GS4",
    "get_grouping_scheme",
    # State inference
    "compute_change_side",
]
