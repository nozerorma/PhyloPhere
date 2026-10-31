# postproc.py — Cluster-train and gene-removal filters shared by the observed chain and the null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Post-processing filters shared by the observed chain and the permulation null.

Imported by: src/utils/null_tables.py (`ctrain`), src/utils/null_merge.py (`gene_removal`), merge_null_block.py (`load_gene_lengths`)
Inputs: positions per labeling, an optional position -> column map, gene lengths (TSV)
Outputs: flagged positions and units to remove (in memory)

One implementation of the two CT_POSTPROC filters, applied per labeling to the pooled
scored rows (observed = labeling ``b_0``):

* Cluster trains: :func:`ctrain` flags every position that lies in an interval whose
  density ``count / span`` reaches ``maxcaas`` (``span >= minlen``), with the span measured
  in the positions or, given a column map, in untrimmed alignment columns. :func:`train_flags`
  is the single entry point that decides over which positions a train is computed.
* Gene removal: :func:`gene_removal` drops (labeling, gene) units whose count of distinct positions exceeds the
  fences of their gene-length stratum (:func:`stratified_fences`), calibrated within each labeling pool:
  ``extreme`` above the far fence, ``dubious`` above the near fence with at least one train position.

Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import csv
import os
import math
from array import array
from typing import Dict, Hashable, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple

__all__ = ["ctrain", "train_flags", "GeneUnit", "gene_unit_stats", "gene_removal", "load_gene_lengths",
           "stratified_fences", "STRATIFIED"]

STRATIFIED = "stratified"  # gene_filter_mode of the length-stratified fences
Unit = Tuple[str, str]  # (labeling, gene)


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


# ── Gene removal ──────────────────────────────────────────────────────────────

class GeneUnit(NamedTuple):
    """One (labeling, gene) unit: its distinct scored positions and whether
    any of them lies in a train."""
    labeling: str
    gene: str
    n_caas: int
    has_clustered: bool


def _percentile_linear(values: Sequence[float], q: float) -> float:
    """Linear-interpolation percentile (numpy / pandas ``quantile`` default), ``q`` in [0, 100]."""
    xs = sorted(float(v) for v in values)
    if not xs:
        return float("nan")
    if len(xs) == 1:
        return xs[0]
    rank = (q / 100.0) * (len(xs) - 1)
    lo = int(rank)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (rank - lo)


def stratified_fences(
    genes: Sequence[str],
    counts: Sequence[int],
    gene_lengths: Mapping[str, float],
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> List[Tuple[Optional[float], Optional[float]]]:
    """(far fence, near fence) of every unit of one pool; (None, None) for a gene without a length.

    The units with a length are ordered by (length, gene) and cut into ``strata`` strata of equal size, so a unit
    is compared with units of its own length. The fence of a stratum is ``Q3 + k * scale`` with
    ``scale = max(IQR of the stratum, IQR of all the units with a length)``: the global IQR keeps the fence from
    collapsing onto Q3 where the counts of a stratum hardly vary.
    """
    fences: List[Tuple[Optional[float], Optional[float]]] = [(None, None)] * len(counts)
    with_length = [i for i, g in enumerate(genes) if (gene_lengths.get(g) or 0) > 0]
    if not with_length:
        return fences
    with_length.sort(key=lambda i: (gene_lengths[genes[i]], genes[i]))
    n = len(with_length)
    pool_iqr = (_percentile_linear([counts[i] for i in with_length], 75.0)
                - _percentile_linear([counts[i] for i in with_length], 25.0))
    members: List[List[int]] = [[] for _ in range(max(1, strata))]
    for rank, i in enumerate(with_length):
        members[rank * len(members) // n].append(i)
    for idx in members:
        if not idx:
            continue
        values = [counts[i] for i in idx]
        q1, q3 = _percentile_linear(values, 25.0), _percentile_linear(values, 75.0)
        scale = max(q3 - q1, pool_iqr)
        fence = (q3 + k_far * scale, q3 + k_near * scale)
        for i in idx:
            fences[i] = fence
    return fences


def gene_unit_stats(
    units: Iterable[GeneUnit],
    gene_lengths: Mapping[str, float],
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> List[dict]:
    """Per-unit statistics and outlier flags, calibrated within each labeling.

    The thresholds are the fences of :func:`stratified_fences`: ``extreme`` is a count above the far fence
    (``threshold_extreme``) and ``dubious`` a count above the near fence (``threshold_dubious``) with a train
    position. Units without a length have ``length``, ``density`` and both thresholds as None.
    """
    pools: Dict[str, List[GeneUnit]] = {}
    for u in units:
        pools.setdefault(u.labeling, []).append(u)

    out: List[dict] = []
    for labeling, members in pools.items():
        fences = stratified_fences([u.gene for u in members], [u.n_caas for u in members], gene_lengths,
                                   strata, k_far, k_near)
        for u, (far, near) in zip(members, fences):
            length = gene_lengths.get(u.gene) if far is not None else None
            out.append({
                "labeling": labeling, "Gene": u.gene,
                "n_caas": u.n_caas,
                "length": length,
                "density": u.n_caas / length * 100.0 if length else None,
                "threshold_dubious": near,
                "threshold_extreme": far,
                "has_clustered": u.has_clustered,
                "dubious": near is not None and u.n_caas > near and u.has_clustered,
                "extreme": far is not None and u.n_caas > far,
            })
    return out


def gene_removal(
    units: Iterable[GeneUnit],
    gene_lengths: Mapping[str, float],
    mode: str = STRATIFIED,
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> Dict[Unit, str]:
    """Units to remove, mapped to their category (``Dubious``, ``Extreme`` or ``Both``).

    ``mode``: ``stratified`` or ``none`` (empty result); any other value is an error, so that a setting of a
    retired mode cannot turn the filter off silently. A unit is removed when its count exceeds the far fence of its
    length stratum (``Extreme``) or the near fence with a train position (``Dubious``); both is ``Both`` (see
    :func:`stratified_fences`).
    """
    mode = (mode or "none").lower()
    if mode == "none":
        return {}
    if mode != STRATIFIED:
        raise ValueError(f"unknown gene filter mode {mode!r}: use 'none' or '{STRATIFIED}'")
    # Same decision rule as gene_unit_stats, but the stream is held in compact per-pool columns
    # (gene, count, train flag) instead of one stats dict per unit: the null has tens of millions
    # of units, and only the removed ones are kept.
    pools: Dict[str, Tuple[List[str], array, bytearray]] = {}
    for u in units:
        genes, counts, flags = pools.setdefault(u.labeling, ([], array("l"), bytearray()))
        genes.append(u.gene)
        counts.append(u.n_caas)
        flags.append(1 if u.has_clustered else 0)

    removed: Dict[Unit, str] = {}
    for labeling, (genes, counts, flags) in pools.items():
        for g, n, fl, (far, near) in zip(genes, counts, flags, stratified_fences(genes, counts, gene_lengths, strata, k_far, k_near)):
            if far is None:
                continue
            ext, dub = n > far, n > near and bool(fl)
            if dub or ext:
                removed[(labeling, g)] = "Both" if dub and ext else "Dubious" if dub else "Extreme"
    return removed


def load_gene_lengths(
    path: Optional[str] = None,
    map_dir: Optional[str] = None,
    map_suffix: str = ".map.tsv",
    genes: Optional[Iterable[str]] = None,
    sources: Optional[Dict[str, str]] = None,
) -> Dict[str, float]:
    """TSV or MAP directory to ``{gene: length}``.

    Prioritizes untrimmed alignment columns from each gene's MAP table in ``map_dir`` (the same
    coordinates used by cluster trains). Falls back to the Ensembl TSV at ``path`` (columns gene, length)
    for any gene missing from ``map_dir`` or when ``map_dir`` is None.
    """
    annotation: Dict[str, float] = {}
    known_genes = set(genes) if genes is not None else None
    if path:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        with open(path, newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            names = {key.strip().lower(): key for key in (reader.fieldnames or [])}
            if not {'gene', 'length'} <= names.keys():
                raise ValueError('Gene annotation needs gene and length columns')
            seen = set()
            for row in reader:
                gene = (row[names['gene']] or '').strip()
                if not gene or gene in seen:
                    raise ValueError(f'Invalid or duplicate annotated gene: {gene}')
                seen.add(gene)
                text = (row[names['length']] or '').strip()
                if text in ('', 'NA'):
                    continue  # An absent optional length is not a malformed number.
                length = float(text)
                if not math.isfinite(length) or length <= 0:
                    raise ValueError(f'Invalid length for {gene}')
                annotation[gene] = length
                if sources is not None:
                    source = 'annotation_units'
                    if {'start','end'} <= names.keys():
                        try:
                            if float(row[names['end']]) - float(row[names['start']]) + 1 == length:
                                source = 'genomic_bp'
                        except (TypeError,ValueError):
                            pass
                    sources[gene] = source
            if known_genes is None:
                known_genes = seen
    out = dict(annotation)
    if map_dir:
        if not os.path.isdir(map_dir):
            raise FileNotFoundError(map_dir)
        from .columns import index_files, find_file, read_map
        index = index_files(map_dir, map_suffix, known_genes)
        for gene in (known_genes if known_genes is not None else index):
            if gene not in index:
                continue
            map_path = find_file(index, gene, map_suffix)
            length = len(read_map(map_path)[0])
            if not length:
                raise ValueError(f'Empty MAP table: {map_path}')
            out[gene] = float(length)
            if sources is not None:
                sources[gene] = 'untrimmed_alignment_columns'
    return out
