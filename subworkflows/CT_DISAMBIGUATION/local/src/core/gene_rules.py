# gene_rules.py — Gene-removal rule and gene-length loading shared by the observed chain and the null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Gene-removal rule shared by the observed chain and the permulation null.

Imported by: src/utils/null_merge.py (`gene_removal`), merge_null_block.py (`read_gene_lengths_table`),
CT_POSTPROC/local/src/filter_caas_genes.py (`read_gene_lengths_table` and the rule),
CT_POSTPROC/local/src/gene_lengths.py (`load_gene_lengths`)
Inputs: (labeling, gene) units with their count of CAAS positions and train flag, gene lengths (TSV or MAP tables)
Outputs: units to remove and per-unit statistics (in memory)

:func:`gene_removal` drops (labeling, gene) units whose count of distinct positions exceeds the fences of their
gene-length stratum (:func:`stratified_fences`), calibrated within each labeling pool (observed = labeling ``b_0``):
``extreme`` above the far fence, ``dubious`` above the near fence with at least one train position. The train
flags come from :mod:`core.trains`.

Pure Python, no numpy/pandas: the null workers import it too.
"""

from __future__ import annotations

import csv
import os
import math
from array import array
from typing import Dict, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Tuple

__all__ = ["GeneUnit", "gene_unit_stats", "gene_removal", "load_gene_lengths", "read_gene_lengths_table",
           "stratified_fences", "LENGTH_TABLE_COLUMNS"]

Unit = Tuple[str, str]  # (labeling, gene)
LENGTH_TABLE_COLUMNS = ("gene", "length", "length_source")  # gene_lengths.tsv; length is NA for a gene without one


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


def _strata(
    genes: Sequence[str],
    counts: Sequence[int],
    gene_lengths: Mapping[str, float],
    strata: int,
    k_far: float,
    k_near: float,
) -> List[dict]:
    """The strata of one pool, in order of length, with the statistics their fences come from.

    The units with a length are ordered by (length, gene) and cut into ``strata`` strata of equal size, so a unit is
    compared with units of its own length. The fence of a stratum is ``Q3 + k * scale`` with
    ``scale = max(IQR of the stratum, IQR of all the units with a length)``: the global IQR keeps the fence from
    collapsing onto Q3 where the counts of a stratum hardly vary. Each entry holds the unit indices (``members``),
    ``stratum`` (1-based), the length range, ``q1``, ``median``, ``q3``, ``iqr``, ``pool_iqr``, ``scale``,
    ``scale_source`` (``stratum`` or ``pool``) and the two fences. Empty without a unit that has a length.
    """
    with_length = [i for i, g in enumerate(genes) if (gene_lengths.get(g) or 0) > 0]
    if not with_length:
        return []
    with_length.sort(key=lambda i: (gene_lengths[genes[i]], genes[i]))
    n = len(with_length)
    pool_iqr = (_percentile_linear([counts[i] for i in with_length], 75.0)
                - _percentile_linear([counts[i] for i in with_length], 25.0))
    members: List[List[int]] = [[] for _ in range(max(1, strata))]
    for rank, i in enumerate(with_length):
        members[rank * len(members) // n].append(i)
    out = []
    for number, idx in enumerate(members, start=1):
        if not idx:
            continue
        values = [counts[i] for i in idx]
        q1, q3 = _percentile_linear(values, 25.0), _percentile_linear(values, 75.0)
        scale = max(q3 - q1, pool_iqr)
        lengths = [gene_lengths[genes[i]] for i in idx]
        out.append({"members": idx, "stratum": number, "length_min": min(lengths), "length_max": max(lengths),
                    "q1": q1, "median": _percentile_linear(values, 50.0), "q3": q3, "iqr": q3 - q1, "pool_iqr": pool_iqr,
                    "scale": scale, "scale_source": "stratum" if q3 - q1 >= pool_iqr else "pool",
                    "far": q3 + k_far * scale, "near": q3 + k_near * scale})
    return out


def stratified_fences(
    genes: Sequence[str],
    counts: Sequence[int],
    gene_lengths: Mapping[str, float],
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> List[Tuple[Optional[float], Optional[float]]]:
    """(far fence, near fence) of every unit of one pool; (None, None) for a gene without a length (see :func:`_strata`)."""
    fences: List[Tuple[Optional[float], Optional[float]]] = [(None, None)] * len(counts)
    for stratum in _strata(genes, counts, gene_lengths, strata, k_far, k_near):
        for i in stratum["members"]:
            fences[i] = (stratum["far"], stratum["near"])
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
        genes = [u.gene for u in members]
        counts = [u.n_caas for u in members]
        detail: Dict[int, dict] = {}
        for stratum in _strata(genes, counts, gene_lengths, strata, k_far, k_near):
            for i in stratum["members"]:
                detail[i] = stratum
        fences = stratified_fences(genes, counts, gene_lengths, strata, k_far, k_near)
        for i, (u, (far, near)) in enumerate(zip(members, fences)):
            length = gene_lengths.get(u.gene) if far is not None else None
            st = detail.get(i, {})
            out.append({
                "labeling": labeling, "Gene": u.gene,
                **{key: st.get(key) for key in ("stratum", "length_min", "length_max", "q1", "median", "q3", "iqr", "pool_iqr",
                                                 "scale", "scale_source")},
                "stratum_genes": len(st["members"]) if st else None,
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
    remove_extreme: bool = True,
    remove_dubious: bool = True,
    strata: int = 10,
    k_far: float = 3.0,
    k_near: float = 1.5,
) -> Dict[Unit, str]:
    """Units to remove, mapped to their category (``Dubious``, ``Extreme`` or ``Both``).

    A unit is Extreme when its count exceeds the far fence of its length stratum, Dubious when it exceeds the near
    fence and carries a train position, Both when both hold (see :func:`stratified_fences`). Each rule has its own
    switch: a unit is removed when it is Extreme and ``remove_extreme`` is on, or Dubious and ``remove_dubious`` is on;
    its category names the marks it carries, whichever switch removed it. Both switches off removes nothing.
    """
    if not (remove_extreme or remove_dubious):
        return {}
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
            if (ext and remove_extreme) or (dub and remove_dubious):
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
    With ``genes``, return only those literal IDs and resolve annotation producer
    suffixes consistently with the observed/null annotation whitelist.
    """
    annotation: Dict[str, float] = {}
    annotation_sources = {}
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
                    annotation_sources[gene] = source
            if known_genes is None:
                known_genes = seen
    from ..data.loaders import annotation_keys
    matched = annotation_keys(known_genes, annotation) if genes is not None else {g: g for g in annotation}
    out = {gene: annotation[key] for gene, key in matched.items()}
    if sources is not None:
        sources.update({gene: annotation_sources[key] for gene, key in matched.items()})
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


def read_gene_lengths_table(path: str) -> Tuple[Dict[str, float], Dict[str, str], int]:
    """gene_lengths.tsv (columns gene, length and optionally length_source) to ``(lengths, sources, n_genes)``.

    The table is the one source of gene lengths of a run: the observed filter and every null block read it, so the
    rule sees the same lengths on both sides. A gene whose length is empty or NA is listed in ``n_genes`` but has no
    entry in ``lengths``, and the rule never flags it.
    """
    lengths: Dict[str, float] = {}
    sources: Dict[str, str] = {}
    n_genes = 0
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        names = {key.strip().lower(): key for key in (reader.fieldnames or [])}
        if not {'gene', 'length'} <= names.keys():
            raise ValueError('Gene lengths table needs gene and length columns')
        seen = set()
        for row in reader:
            gene = (row[names['gene']] or '').strip()
            if not gene or gene in seen:
                raise ValueError(f'Invalid or duplicate gene in the lengths table: {gene}')
            seen.add(gene)
            n_genes += 1
            text = (row[names['length']] or '').strip()
            if text in ('', 'NA'):
                continue
            length = float(text)
            if not math.isfinite(length) or length <= 0:
                raise ValueError(f'Invalid length for {gene}')
            lengths[gene] = length
            if 'length_source' in names:
                sources[gene] = (row[names['length_source']] or '').strip() or 'table'
    return lengths, sources, n_genes
