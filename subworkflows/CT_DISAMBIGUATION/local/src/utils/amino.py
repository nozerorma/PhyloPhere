"""Amino-acid utility helpers shared across pipeline layers.

This module intentionally sits outside the ``data`` and ``convergence``
packages so both can depend on it without creating circular imports.
"""

from typing import List, Optional


def normalize_amino_list(values: List[Optional[str]]) -> List[str]:
    """Clean and deduplicate a list of amino-acid codes.

    Normalization rules:
    - cast to uppercase strings
    - drop empty/unknown/gap tokens: ``""``, ``-``, ``X``, ``?``
    - preserve first-seen order while removing duplicates
    """
    cleaned: List[str] = []
    seen = set()

    for value in values or []:
        if value is None:
            continue

        aa = str(value).strip().upper()
        if aa in {"", "-", "X", "?"} or aa in seen:
            continue

        seen.add(aa)
        cleaned.append(aa)

    return cleaned
