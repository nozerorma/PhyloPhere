# scores.py — Position and gene CAAS scores from the path scores of the detections.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Position and gene CAAS scores, shared by the observed chain (labeling b_0) and the null.

Imported by: src/utils/null_tables.py, src/utils/null_merge.py, reduce_null_blocks.py (`cycle_order_key`), the SCORING scripts observed_core_scores.py and null_position_stats.py (staged next to them)
Inputs: per-side score mappings and the reference pools of position scores
Outputs: position scores, per-direction values and gene scores (in memory)

* The score of a row (``caas_row``) is its ``asr_path_score`` (identity). A position's score, per side,
  is the score of its row, in (0, 1]. A position that is not detected, or detected with a score of 0, has no
  score; a side without a score does not enter the direction values nor the positive-position count of the gene.
* Directions: ``top`` and ``bottom`` use only that side's rows; ``all`` keeps one entry per
  position, its best side, so a position detected on both sides is not counted twice.
* Gene score: the largest positive position score in each direction. No ECDF
  or adjustment for the number of positions is applied. None means no score.

Pure Python. The same inputs give the same bits wherever the functions are called.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

__all__ = ["DIRECTIONS", "GENE_STATS", "TIE_TOL", "SCORE_RULE", "position_score", "collapse_sides", "direction_values",
           "gene_stat_value", "gene_scores", "cycle_order_key"]

DIRECTIONS = ("all", "top", "bottom")

# Only the raw maximum is supported.
GENE_STATS = ("max",)

# Rule of the position score. It is stamped on the null (counts.json, score_rule) so that
# null_position_stats.py refuses a null scored with another rule. "caas_pos_species_pss": the CAAS score (mean over the domains of the PSS-weighted mean, over the species of the domain, of their hypothesis scores), defined only when positive.
SCORE_RULE = "caas_pos_species_pss_v4"

# Absolute; scores lie in [0, 1]. Rounding noise of a weighted mean is ~1e-16, and real
# differences between position scores are orders of magnitude larger.
TIE_TOL = 1e-12


def cycle_order_key(cycle: str):
    """Replay order of a cycle: the numeric suffix ("b_2" before "b_10"), the label for cycles without one.

    The order fixes the checkpoints of the saturation curve of the null and which cycle is the first to score a position;
    the cycles are exchangeable otherwise.
    """
    tail = cycle.rsplit("_", 1)[-1]
    return (0, int(tail), cycle) if tail.isdigit() else (1, 0, cycle)


def _is_value(x) -> bool:
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def position_score(caas_row: Optional[float]) -> Optional[float]:
    """Score of a position side: its ``caas_row``; None when it has none or it is <= 0.0."""
    return float(caas_row) if (_is_value(caas_row) and float(caas_row) > 0.0) else None


def collapse_sides(side_scores: Mapping[str, float]) -> Dict[str, float]:
    """Per-direction entries of one position: ``all`` = best side, plus ``top`` / ``bottom`` when present."""
    vals = {s: v for s, v in side_scores.items() if _is_value(v) and v > 0.0}
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


def gene_stat_value(values: Sequence[float], pool_sorted: Sequence[float], stat: str = "max") -> Optional[float]:
    """Largest positive position score. pool_sorted is ignored for legacy callers.

    None means no scored position; it is distinct from the zero used to fill a
    canonical null matrix when a gene has no detection in a cycle.
    """
    if stat != "max":
        raise ValueError("Only gene_stat='max' is supported; rebuild legacy gene statistics")
    xs = [v for v in values if _is_value(v) and v > 0.0]
    return max(xs) if xs else None


def gene_scores(
    positions: Iterable[Mapping[str, float]],
    pools: Mapping[str, Sequence[float]],
    stat: str = "max",
) -> Dict[str, Optional[float]]:
    """The gene statistic ``stat`` per direction for one gene, from its positions (``{side: score}`` each)."""
    vals = direction_values(positions)
    return {d: gene_stat_value(vals[d], pools[d], stat) for d in DIRECTIONS}
