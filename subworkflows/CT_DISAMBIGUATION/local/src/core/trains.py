# trains.py — Cluster-train detection shared by the observed chain and the null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Cluster trains, shared by the observed chain and the permulation null.

Imported by: src/utils/null_tables.py (`ctrain`), CT_POSTPROC/local/src/filter_caas_clusters-param.py (`ctrain`)
Inputs: positions per labeling, an optional position -> column map
Outputs: flagged positions (in memory)

:func:`ctrain` flags every position that lies in an interval whose density ``count / span`` reaches
``maxcaas`` (``span >= minlen``), with the span measured in the positions or, given a column map, in
untrimmed alignment columns. :func:`train_flags` is the single entry point that decides over which
positions a train is computed. The gene-removal rule that consumes the trains is in :mod:`core.gene_rules`.

Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import math
from typing import Dict, Hashable, Iterable, List, Mapping, Optional, Sequence, Set

__all__ = ["ctrain", "train_flags"]


# ── Cluster trains ────────────────────────────────────────────────────────────

def ctrain(
    positions: Sequence[int],
    maxcaas: float = 0.7,
    minlen: int = 3,
    columns: Optional[Mapping[int, int]] = None,
    min_positions: Optional[int] = None,
) -> List[int]:
    """Sorted positions inside a high-density interval.

    For the sorted unique positions, every interval ``[l, r]`` (indices) with
    ``span = end - start + 1 >= minlen`` and ``count / span >= maxcaas`` flags all of its
    positions. ``columns`` maps each position to the column the span is measured in (see
    :mod:`core.columns`): with it, columns that lie between two positions and hold none (the ones
    the trimmer removed) count towards the span; without it the positions themselves are the
    coordinates. An explicit ``min_positions`` applies to each interval. None retains
    the legacy rule: at least ``minlen`` candidates in the gene, with no separate
    interval-count requirement.
    """
    legacy = min_positions is None
    minimum_gene_count = minlen if legacy else min_positions
    min_positions = 1 if legacy else min_positions
    if (not isinstance(minlen, int) or minlen < 1 or not isinstance(min_positions, int)
            or min_positions < 1 or not math.isfinite(maxcaas) or not 0 <= maxcaas <= 1):
        raise ValueError('Invalid train span, position count or density')
    if any(not math.isfinite(p) or p < 0 or p != int(p) for p in positions):
        raise ValueError('Invalid train position')
    uniq = sorted(set(positions))
    if columns is None:
        coords = uniq
    else:
        missing = [p for p in uniq if p not in columns]
        if missing:
            raise ValueError(f'position {missing[0]} has no column in the map')
        uniq = sorted(uniq, key=columns.__getitem__)
        coords = [columns[p] for p in uniq]
        if len(set(coords)) != len(coords):
            raise ValueError('two positions share a column')
    n = len(uniq)
    if n < minimum_gene_count:
        return []
    full_span = coords[-1] - coords[0] + 1
    if full_span >= minlen and n / full_span >= maxcaas:
        return sorted(uniq)
    # Mark intervals with a difference array: no repeated O(interval length)
    # set updates. O(n²) time and O(n) additional memory, identical decisions.
    marks = [0] * (n + 1)
    for r in range(min_positions - 1, n):
        for l in range(r - min_positions + 2):
            span = coords[r] - coords[l] + 1
            if span < minlen:
                break
            if (r - l + 1) / span >= maxcaas:
                marks[l] += 1
                marks[r + 1] -= 1
    active, bad = 0, []
    for i, pos in enumerate(uniq):
        active += marks[i]
        if active:
            bad.append(pos)
    return sorted(bad)


def train_flags(
    pos_by_key: Mapping[Hashable, Iterable[int]],
    maxcaas: float = 0.7,
    minlen: int = 3,
    columns: Optional[Mapping[int, int]] = None,
    min_positions: Optional[int] = None,
) -> Dict[Hashable, Set[int]]:
    """Train positions per key, for one gene.

    The key names what a train is computed over, e.g. the labeling: the
    positions filed under one key form one train universe. Callers choose the grain by
    how they file positions; this is the only place trains are computed. ``columns`` is the
    gene's position -> untrimmed column map (:func:`ctrain`), or None to measure in the
    positions themselves.
    """
    return {key: set(ctrain(list(pos), maxcaas, minlen, columns, min_positions)) for key, pos in pos_by_key.items()}
