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
* Gene removal: :func:`gene_removal` drops (labeling, gene) units by
  ``dubious`` (IQR outlier in distinct positions AND at least one train position) and
  ``extreme`` (density above a percentile), each calibrated within its own
  labeling pool. In ``stratified`` mode both are replaced by one rule on the count against the fences of
  the unit's length stratum (:func:`stratified_fences`): ``extreme`` above the far fence, ``dubious`` above the
  near fence with a train position.

Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import csv
import os
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
) -> List[int]:
    """Sorted positions inside a high-density interval.

    For the sorted unique positions, every interval ``[l, r]`` (indices) with
    ``span = end - start + 1 >= minlen`` and ``count / span >= maxcaas`` flags all of its
    positions. ``columns`` maps each position to the column the span is measured in (see
    :mod:`core.columns`): with it, columns that lie between two positions and hold none (the ones
    the trimmer removed) count towards the span; without it the positions themselves are the
    coordinates. A unit with fewer than ``minlen`` positions has no train either way.
    """
    uniq = sorted({int(p) for p in positions})
    if columns is None:
        coords = uniq
    else:
        missing = [p for p in uniq if p not in columns]
        if missing:
            raise ValueError(f"position {missing[0]} has no column in the map")
        uniq = sorted(uniq, key=columns.__getitem__)
        coords = [columns[p] for p in uniq]
        if len(set(coords)) != len(coords):
            raise ValueError("two positions share a column")
    n = len(uniq)
    if n < minlen:
        return []
    bad: Set[int] = set()
    for r in range(n):
        for l in range(r + 1):
            span = coords[r] - coords[l] + 1
            if span >= minlen and (r - l + 1) / span >= maxcaas:
                bad.update(uniq[l:r + 1])
    return sorted(bad)


def train_flags(
    pos_by_key: Mapping[Hashable, Iterable[int]],
    maxcaas: float = 0.7,
    minlen: int = 3,
    columns: Optional[Mapping[int, int]] = None,
) -> Dict[Hashable, Set[int]]:
    """Train positions per key, for one gene.

    The key names what a train is computed over, e.g. the labeling: the
    positions filed under one key form one train universe. Callers choose the grain by
    how they file positions; this is the only place trains are computed. ``columns`` is the
    gene's position -> untrimmed column map (:func:`ctrain`), or None to measure in the
    positions themselves.
    """
    return {key: set(ctrain(list(pos), maxcaas, minlen, columns)) for key, pos in pos_by_key.items()}


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
    iqr_multiplier: float = 3.0,
    extreme_percentile: float = 0.99,
    mode: str = "both",
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> List[dict]:
    """Per-unit statistics and outlier flags, calibrated within each labeling.

    ``dubious`` pool: every unit (it needs no length); threshold ``Q3 + k * IQR`` of the
    distinct-position counts; flagged when ``n_caas > threshold`` and the unit has a train
    position. ``extreme`` pool: units whose gene has a positive length; density
    ``n_caas / length * 100``; flagged when it exceeds the ``extreme_percentile`` quantile.
    Units without a length have ``length``, ``density`` and ``threshold_extreme`` as None.

    With ``mode == "stratified"`` the thresholds are the fences of :func:`stratified_fences`: ``extreme`` is the
    far fence (``threshold_extreme``) and ``dubious`` the near fence (``threshold_dubious``) with a train position.
    """
    units = list(units)
    pools: Dict[str, List[GeneUnit]] = {}
    for u in units:
        pools.setdefault(u.labeling, []).append(u)

    out: List[dict] = []
    for labeling, members in pools.items():
        if (mode or "").lower() == STRATIFIED:
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
            continue
        counts = [u.n_caas for u in members]
        q1, q3 = _percentile_linear(counts, 25.0), _percentile_linear(counts, 75.0)
        thr_dub = q3 + iqr_multiplier * (q3 - q1)

        dens = {}
        for u in members:
            length = gene_lengths.get(u.gene)
            if length and length > 0:
                dens[u.gene] = u.n_caas / length * 100.0
        thr_ext = _percentile_linear(list(dens.values()), extreme_percentile * 100.0) if dens else None

        for u in members:
            d = dens.get(u.gene)
            out.append({
                "labeling": labeling, "Gene": u.gene,
                "n_caas": u.n_caas,
                "length": gene_lengths.get(u.gene) if d is not None else None,
                "density": d,
                "threshold_dubious": thr_dub,
                "threshold_extreme": thr_ext,
                "has_clustered": u.has_clustered,
                "dubious": u.n_caas > thr_dub and u.has_clustered,
                "extreme": d is not None and d > thr_ext,
            })
    return out


def gene_removal(
    units: Iterable[GeneUnit],
    gene_lengths: Mapping[str, float],
    mode: str = "dubious",
    iqr_multiplier: float = 3.0,
    extreme_percentile: float = 0.99,
    thresholds: Optional[Dict[str, Tuple[float, Optional[float]]]] = None,
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> Dict[Unit, str]:
    """Units to remove, mapped to their category (``Dubious``, ``Extreme`` or ``Both``).

    ``mode``: ``dubious``, ``extreme``, ``both``, ``stratified`` or ``none`` (empty result).
    ``stratified`` removes a unit whose count exceeds the far fence of its length stratum (``Extreme``), or the near
    fence with a train position (``Dubious``); both is ``Both`` (see :func:`stratified_fences`; the fences vary with
    the gene, so ``thresholds`` stays empty).
    ``thresholds``, when given, is filled with (count threshold, density threshold) of every labeling; it does not
    change the decisions.
    """
    mode = (mode or "none").lower()
    if mode not in ("dubious", "extreme", "both", STRATIFIED):
        return {}
    want_dub, want_ext = mode in ("dubious", "both"), mode in ("extreme", "both")
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
        if mode == STRATIFIED:
            for g, n, fl, (far, near) in zip(genes, counts, flags, stratified_fences(genes, counts, gene_lengths, strata, k_far, k_near)):
                if far is None:
                    continue
                ext, dub = n > far, n > near and bool(fl)
                if dub or ext:
                    removed[(labeling, g)] = "Both" if dub and ext else "Dubious" if dub else "Extreme"
            continue
        q1, q3 = _percentile_linear(counts, 25.0), _percentile_linear(counts, 75.0)
        thr_dub = q3 + iqr_multiplier * (q3 - q1)
        dens = []
        for g, n in zip(genes, counts):
            length = gene_lengths.get(g)
            dens.append(n / length * 100.0 if length and length > 0 else None)
        valid = [d for d in dens if d is not None]
        thr_ext = _percentile_linear(valid, extreme_percentile * 100.0) if valid else None
        if thresholds is not None:
            thresholds[labeling] = (thr_dub, thr_ext)
        for g, n, fl, d in zip(genes, counts, flags, dens):
            dub = want_dub and n > thr_dub and bool(fl)
            ext = want_ext and d is not None and d > thr_ext
            if dub or ext:
                removed[(labeling, g)] = "Both" if dub and ext else "Dubious" if dub else "Extreme"
    return removed


def load_gene_lengths(
    path: Optional[str] = None,
    map_dir: Optional[str] = None,
    map_suffix: str = ".map.tsv",
    genes: Optional[Iterable[str]] = None,
) -> Dict[str, float]:
    """TSV or MAP directory to ``{gene: length}``.

    Prioritizes untrimmed alignment columns from each gene's MAP table in ``map_dir`` (the same
    coordinates used by cluster trains). Falls back to the Ensembl TSV at ``path`` (columns gene, length)
    for any gene missing from ``map_dir`` or when ``map_dir`` is None.
    """
    out: Dict[str, float] = {}

    # 1. Primary source: untrimmed alignment lengths from MAP directory
    if map_dir and os.path.isdir(map_dir):
        from .columns import index_files, map_untrimmed_length
        m_index = index_files(map_dir, map_suffix)
        target_genes = genes if genes is not None else m_index.keys()
        for g in target_genes:
            map_path = m_index.get(g)
            if map_path and os.path.isfile(map_path):
                length = map_untrimmed_length(map_path)
                if length is not None and length > 0:
                    out[g] = float(length)

    # 2. Fallback source: Ensembl gene annotation TSV
    if path and os.path.isfile(path):
        with open(path, "r", newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header:
                lc = [h.strip().lower() for h in header]
                if "gene" in lc and "length" in lc:
                    gi, li = lc.index("gene"), lc.index("length")
                    for row in reader:
                        if len(row) <= max(gi, li):
                            continue
                        gene_name = row[gi].strip()
                        if gene_name not in out:
                            try:
                                val = float(row[li])
                                if val > 0:
                                    out[gene_name] = val
                            except (TypeError, ValueError):
                                continue
    return out
