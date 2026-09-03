#!/usr/bin/env python3
"""Cycle-aware CT_POSTPROC filtering for the CAAS permulation null (Gap B).

Pure, dependency-light ports of the observed post-processing steps so the
per-cycle null CAAS candidate pool is filtered the SAME way the observed
``filtered_discovery.tsv`` is before ``scoring_compute.R`` scores it:

* :func:`ctrain` / :func:`clustering_discards` — verbatim algorithm of
  ``CT_POSTPROC/local/src/filter_caas_clusters-param.py`` (density "trains").
* :func:`cycle_gene_removal` — the dubious + extreme gene logic of
  ``CT_POSTPROC/local/src/filter_caas_genes.py``, with the base-cycle tag playing
  the role the observed side gives to ``trait`` (one calibration per cycle,
  matching ``build_percent_rank_lookup`` / ``_build_cycle_score_pools``).

None of these touch ASR. They operate on the stacked
``perm_pos_detail`` table (+ per-cycle position lists inside ``_perms_worker``).

Author: PhyloPhere pipeline (Gap B)
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Set, Tuple

__all__ = ["ctrain", "clustering_discards", "cycle_gene_removal"]


# ── Cluster ("train") detection ──────────────────────────────────────────────

def ctrain(positions: Sequence[int], maxcaas: float = 0.7, minlen: int = 3) -> List[int]:
    """Positions inside a high-density CAAS interval.

    Verbatim port of ``filter_caas_clusters-param.py::ctrain``: for the sorted
    unique position list, every interval ``[l, r]`` with
    ``span = end - start + 1 >= minlen`` and ``count / span >= maxcaas`` flags
    all of its positions.  Returns the sorted list of flagged positions.
    """
    uniq = sorted({int(p) for p in positions})
    n = len(uniq)
    if n < minlen:
        return []
    bad: Set[int] = set()
    for r in range(n):
        for l in range(r + 1):
            start, end = uniq[l], uniq[r]
            span = end - start + 1
            if span < minlen:
                continue
            count = r - l + 1
            if count / span >= maxcaas:
                bad.update(uniq[l : r + 1])
    return sorted(bad)


def clustering_discards(
    pos_by_group: Dict[str, Iterable[int]],
    maxcaas: float = 0.7,
    minlen: int = 3,
) -> Dict[str, Set[int]]:
    """Per caap_group set of clustered ("Discarded") positions for one cycle/gene.

    Mirrors the group-aware branch of ``filterCAAS`` (each caap_group processed
    independently).
    """
    return {
        grp: set(ctrain(list(pos), maxcaas, minlen))
        for grp, pos in pos_by_group.items()
    }


# ── Gene-level filter (dubious + extreme), per base cycle ────────────────────

def _percentile_linear(values: Sequence[float], q: float) -> float:
    """``numpy.percentile(values, q)`` default (linear) without importing numpy.

    ``q`` in [0, 100].  Matches pandas ``Series.quantile`` default interpolation,
    which is what ``filter_caas_genes.py`` uses.
    """
    xs = sorted(float(v) for v in values)
    if not xs:
        return float("nan")
    if len(xs) == 1:
        return xs[0]
    rank = (q / 100.0) * (len(xs) - 1)
    lo = int(rank)
    hi = min(lo + 1, len(xs) - 1)
    frac = rank - lo
    return xs[lo] + (xs[hi] - xs[lo]) * frac


def _iqr_threshold(values: Sequence[float], multiplier: float) -> float:
    q1 = _percentile_linear(values, 25.0)
    q3 = _percentile_linear(values, 75.0)
    return q3 + multiplier * (q3 - q1)


def cycle_gene_removal(
    rows: Iterable[Tuple[str, str, str, int, bool]],
    gene_lengths: Dict[str, float],
    mode: str = "dubious",
    iqr_multiplier: float = 3.0,
    extreme_percentile: float = 0.99,
) -> Set[Tuple[str, str, str]]:
    """Gene/group units to drop from each base cycle's null CAAS pool.

    ``rows`` — iterable of ``(cycle, caap_group, gene, n_caas, has_clustered)``
    where ``n_caas`` = distinct detected Positions for that (cycle, group, gene)
    and ``has_clustered`` = the unit carries ≥1 ``ctrain``-flagged position.

    Returns the set of ``(cycle, caap_group, gene)`` to remove.  Mirrors
    ``filter_caas_genes.py``:

    * **dubious** — per (cycle, group): ``n_caas > Q3 + k*IQR`` AND
      ``has_clustered``.
    * **extreme** — per (cycle, group): ``density > quantile(extreme_percentile)``
      with ``density = n_caas / length * 100`` (genes lacking a length are
      excluded, matching ``dropna(subset=['length'])``).
    * **both** — union.  **none** — empty.
    """
    mode = (mode or "none").lower()
    if mode not in ("dubious", "extreme", "both"):
        return set()

    # group rows by (cycle, caap_group)
    by_cg: Dict[Tuple[str, str], List[Tuple[str, int, bool]]] = {}
    for cycle, grp, gene, n_caas, has_clustered in rows:
        by_cg.setdefault((cycle, grp), []).append((gene, int(n_caas), bool(has_clustered)))

    remove: Set[Tuple[str, str, str]] = set()
    want_dub = mode in ("dubious", "both")
    want_ext = mode in ("extreme", "both")

    for (cycle, grp), units in by_cg.items():
        counts = [n for _g, n, _c in units]
        if want_dub and len(counts) >= 2:
            thr = _iqr_threshold(counts, iqr_multiplier)
            for gene, n_caas, has_clustered in units:
                if n_caas > thr and has_clustered:
                    remove.add((cycle, grp, gene))
        if want_ext:
            dens = [
                (gene, n / gene_lengths[gene] * 100.0)
                for gene, n, _c in units
                if gene in gene_lengths and gene_lengths[gene]
            ]
            if len(dens) >= 2:
                ethr = _percentile_linear([d for _g, d in dens], extreme_percentile * 100.0)
                for gene, d in dens:
                    if d > ethr:
                        remove.add((cycle, grp, gene))
    return remove


def load_gene_lengths(path: str) -> Dict[str, float]:
    """TSV → ``{gene: length}``.  Case-insensitive ``gene`` / ``length`` columns
    (ensembl annotation format used by ``filter_caas_genes.py``)."""
    import csv

    out: Dict[str, float] = {}
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            return out
        lc = [h.strip().lower() for h in header]
        try:
            gi = lc.index("gene")
            li = lc.index("length")
        except ValueError:
            return out
        for row in reader:
            if len(row) <= max(gi, li):
                continue
            try:
                out[row[gi].strip()] = float(row[li])
            except (TypeError, ValueError):
                continue
    return out
