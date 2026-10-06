# scores.py — Position and gene CAAS scores from per-scheme path scores.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Position and gene CAAS scores, shared by the observed chain (labeling b_0) and the null.

Imported by: src/utils/gene_wrapper.py
Inputs: per-scheme (or per-side) score mappings and the reference pools of position scores
Outputs: position scores, per-direction values and gene scores (in memory)

* The score of a row (``caas_row``) is its ``asr_path_score`` (identity), so a position's score is the
  mean of ``asr_path_score`` over the caap_group schemes that detected it, per side.
* Directions: ``top`` and ``bottom`` use only that side's rows; ``all`` keeps one entry per
  position, its best side, so a position detected on both sides is not counted twice.
* Gene score: ``size_adj_max(x) = F(max(x)) ** len(x)`` with F the ECDF of the reference
  pool of the same direction (and, in the null, of the same labeling). It is None when the
  gene has no scored position in that direction or the pool is empty: no positions, no score.
  F counts pool values up to ``max + TIE_TOL``: position scores are means of a few values in
  [0, 1], so means that are equal in exact arithmetic can differ by rounding noise (~1e-16),
  and the pool is heavily tied. The tolerance makes those ties deterministic.

Pure Python. Sums are correctly rounded (``math.fsum``), so the result does not depend on
the order the schemes were filled in and the same inputs give the same bits wherever the
function is called.
"""

from __future__ import annotations

import bisect
import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

__all__ = ["DIRECTIONS", "TIE_TOL", "position_sum", "position_score", "collapse_sides",
           "direction_values", "size_adj_max", "gene_scores"]

DIRECTIONS = ("all", "top", "bottom")

# Absolute; scores lie in [0, 1]. Rounding noise of a mean of <= 5 values is ~1e-16, and real
# differences between position scores are orders of magnitude larger.
TIE_TOL = 1e-12


def _is_value(x) -> bool:
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def position_sum(scheme_scores: Mapping[str, float]) -> Tuple[float, int]:
    """(sum, n) of a position's per-scheme scores, skipping missing values.

    The sum is correctly rounded, hence independent of the order of the mapping.
    """
    vals = [v for v in scheme_scores.values() if _is_value(v)]
    return math.fsum(vals), len(vals)


def position_score(scheme_scores: Mapping[str, float]) -> Optional[float]:
    """Mean over the schemes that scored the position; None when none did."""
    total, n = position_sum(scheme_scores)
    return total / n if n else None


def collapse_sides(side_scores: Mapping[str, float]) -> Dict[str, float]:
    """Per-direction entries of one position: ``all`` = best side, plus ``top`` / ``bottom`` when present."""
    vals = {s: v for s, v in side_scores.items() if _is_value(v)}
    if not vals:
        return {}
    out = {"all": max(vals.values())}
    for d in ("top", "bottom"):
        if d in vals:
            out[d] = vals[d]
    return out


def direction_values(positions: Iterable[Mapping[str, float]]) -> Dict[str, List[float]]:
    """Lists of position scores per direction, from an iterable of ``{side: score}`` (one per position)."""
    out: Dict[str, List[float]] = {d: [] for d in DIRECTIONS}
    for sides in positions:
        for d, v in collapse_sides(sides).items():
            out[d].append(v)
    return out


def size_adj_max(values: Sequence[float], pool_sorted: Sequence[float]) -> Optional[float]:
    """``(#{pool <= max(values)} / len(pool)) ** len(values)``; None if either input is empty.

    ``pool_sorted`` must be ascending. Pool values within ``TIE_TOL`` above the maximum are
    ties and count as below-or-equal.
    """
    xs = [v for v in values if _is_value(v)]
    if not xs or len(pool_sorted) == 0:
        return None
    return (bisect.bisect_right(pool_sorted, max(xs) + TIE_TOL) / len(pool_sorted)) ** len(xs)


def gene_scores(
    positions: Iterable[Mapping[str, float]],
    pools: Mapping[str, Sequence[float]],
) -> Dict[str, Optional[float]]:
    """``size_adj_max`` per direction for one gene, from its positions (``{side: score}`` each)."""
    vals = direction_values(positions)
    return {d: size_adj_max(vals[d], pools[d]) for d in DIRECTIONS}
