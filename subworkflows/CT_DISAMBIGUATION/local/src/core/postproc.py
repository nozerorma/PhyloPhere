"""Post-processing filters shared by the observed chain and the permulation null.

One implementation of the two CT_POSTPROC filters, applied per labeling to the pooled
scored rows (observed = labeling ``b_0``):

* Cluster trains: :func:`ctrain` flags every position that lies in an interval whose
  density ``count / span`` reaches ``maxcaas`` (``span >= minlen``), with the span measured
  in the positions or, given a column map, in untrimmed alignment columns. :func:`train_flags`
  is the single entry point that decides over which positions a train is computed.
* Gene removal: :func:`gene_removal` drops (labeling, caap_group, gene) units by
  ``dubious`` (IQR outlier in distinct positions AND at least one train position) and
  ``extreme`` (density above a percentile), each calibrated within its own
  (labeling, caap_group) pool.

Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import csv
from typing import Dict, Hashable, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple

__all__ = ["ctrain", "train_flags", "GeneUnit", "gene_unit_stats", "gene_removal", "load_gene_lengths"]

Unit = Tuple[str, str, str]  # (labeling, caap_group, gene)


# ── Cluster trains ───────────────────────────────────────────────────────────

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

    The key names what a train is computed over, e.g. ``(labeling, caap_group)``: the
    positions filed under one key form one train universe. Callers choose the grain by
    how they file positions; this is the only place trains are computed. ``columns`` is the
    gene's position -> untrimmed column map (:func:`ctrain`), or None to measure in the
    positions themselves.
    """
    return {key: set(ctrain(list(pos), maxcaas, minlen, columns)) for key, pos in pos_by_key.items()}


# ── Gene removal ─────────────────────────────────────────────────────────────

class GeneUnit(NamedTuple):
    """One (labeling, caap_group, gene) unit: its distinct scored positions and whether
    any of them lies in a train."""
    labeling: str
    caap_group: str
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


def gene_unit_stats(
    units: Iterable[GeneUnit],
    gene_lengths: Mapping[str, float],
    iqr_multiplier: float = 3.0,
    extreme_percentile: float = 0.99,
) -> List[dict]:
    """Per-unit statistics and outlier flags, calibrated within (labeling, caap_group).

    ``dubious`` pool: every unit (it needs no length); threshold ``Q3 + k * IQR`` of the
    distinct-position counts; flagged when ``n_caas > threshold`` and the unit has a train
    position. ``extreme`` pool: units whose gene has a positive length; density
    ``n_caas / length * 100``; flagged when it exceeds the ``extreme_percentile`` quantile.
    Units without a length have ``length``, ``density`` and ``threshold_extreme`` as None.
    """
    units = list(units)
    pools: Dict[Tuple[str, str], List[GeneUnit]] = {}
    for u in units:
        pools.setdefault((u.labeling, u.caap_group), []).append(u)

    out: List[dict] = []
    for (labeling, group), members in pools.items():
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
                "labeling": labeling, "caap_group": group, "Gene": u.gene,
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
) -> Dict[Unit, str]:
    """Units to remove, mapped to their category (``Dubious``, ``Extreme`` or ``Both``).

    ``mode``: ``dubious``, ``extreme``, ``both`` or ``none`` (empty result).
    """
    mode = (mode or "none").lower()
    if mode not in ("dubious", "extreme", "both"):
        return {}
    want_dub, want_ext = mode in ("dubious", "both"), mode in ("extreme", "both")
    removed: Dict[Unit, str] = {}
    for s in gene_unit_stats(units, gene_lengths, iqr_multiplier, extreme_percentile):
        dub, ext = want_dub and s["dubious"], want_ext and s["extreme"]
        if dub or ext:
            removed[(s["labeling"], s["caap_group"], s["Gene"])] = (
                "Both" if dub and ext else "Dubious" if dub else "Extreme")
    return removed


def load_gene_lengths(path: str) -> Dict[str, float]:
    """TSV to ``{gene: length}``; case-insensitive ``gene`` and ``length`` columns."""
    out: Dict[str, float] = {}
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            return out
        lc = [h.strip().lower() for h in header]
        if "gene" not in lc or "length" not in lc:
            return out
        gi, li = lc.index("gene"), lc.index("length")
        for row in reader:
            if len(row) <= max(gi, li):
                continue
            try:
                out[row[gi].strip()] = float(row[li])
            except (TypeError, ValueError):
                continue
    return out
