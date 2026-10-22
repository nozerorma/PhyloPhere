# null_merge.py — Merge of the tile tables of one block of cycles into the gene x cycle statistics of the null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Block merge of the permutation null.

Imported by: merge_null_block.py (CAAS_CORE_MERGE, one task per block of cycles)
Inputs: the ``units`` and ``scores`` tables of every tile of one block (src/utils/null_tables.py), the cycles of the
        block, the gene lengths and the gene-removal settings
Outputs: per block, the genes x cycles matrices of the gene statistics, the removed units, the roster, the per-cycle
         sample and quantiles of the scores, the per-cycle summary, the table of scored positions and the number of
         cycles that detect each position (``detections``)

Why a block and not a gene batch: gene removal and the reference pools are decided within a cycle, over every gene of
that cycle (core.postproc.gene_removal calibrates each cycle against its own units, and the gene score is calibrated
against the pool of the positions of its cycle), so a cycle must be merged with all its genes together. A cycle belongs
to one block, so blocks are independent and the memory of a merge follows the block, not the number of cycles of the run.

Gene statistic, per (gene, cycle, direction) with n scored positions and largest score ``max``:
    max           the largest position score of the gene in the direction
    fmax          F(max), F the ECDF of the position scores of the cycle in the direction: the count of pool values up to
                  ``max + TIE_TOL`` over the pool size (core.scores.TIE_TOL)
    size_adj_max  F(max) ** n
The pool of a cycle holds the scored positions of the units that gene removal keeps; the ``all`` pool holds one entry per
position (its best side). A unit that is removed, or has no scored position in a direction, has no value there (0 in the
matrices). F(max) ** n is computed with Python floats, the same expression as core.scores.size_adj_max, so its bits do
not depend on a vectorised power.

The sample and the quantiles describe the rows of the null (scored or not) that remain after removal: a uniform sample
without replacement of ``sample_per_cycle`` rows per cycle, drawn with a generator seeded by the pipeline seed and the
cycle. The rows without a score are not stored in the tables, so a sampled row of that kind carries its gene and a
score of 0 and no position (``NA``).
"""

from __future__ import annotations

import csv
import json
import logging
import random
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc

from src.core.postproc import GeneUnit, gene_removal
from src.core.scores import SCORE_RULE, TIE_TOL, cycle_order_key
from src.utils.null_tables import SCORES_SCHEMA, UNITS_SCHEMA, read_scores, read_units

logger = logging.getLogger(__name__)

DIRECTIONS = {"global": "all", "top": "top", "bottom": "bottom"}   # matrix name -> direction of core.scores
GENE_STATS = ("size_adj_max", "fmax", "max")
QUANT_LEVELS = [5, 10, 25, 50, 75, 90, 95]
HIST_BINS = 101   # score bins of width 0.01: floor(100 * score), the score 1.0 in the last one


def _vocabulary_codes(column: pa.ChunkedArray, vocabulary: Sequence[str]) -> np.ndarray:
    """Index of each value of a string column in `vocabulary` (which must contain them all)."""
    codes = pc.index_in(column.combine_chunks() if column.num_chunks else pa.array([], pa.string()),
                        value_set=pa.array(list(vocabulary), pa.string()))
    arr = codes.to_numpy(zero_copy_only=False)
    if np.isnan(arr.astype(float)).any():
        raise ValueError("a value is missing from its vocabulary")
    return arr.astype(np.int64)


class BlockTables:
    """Units and scores of one block, as numpy arrays indexed by gene and cycle codes."""

    def __init__(self, tile_dirs: Sequence[Path], cycles: Sequence[str], universe: Optional[Iterable[str]] = None):
        units = [read_units(d) for d in tile_dirs]
        scores = [read_scores(d) for d in tile_dirs]
        units = pa.concat_tables(units) if units else UNITS_SCHEMA.empty_table()
        scores = pa.concat_tables(scores) if scores else SCORES_SCHEMA.empty_table()
        self.cycles = list(cycles)
        stray = sorted(set(pc.unique(units["cycle"]).to_pylist()) - set(self.cycles))
        if stray:
            raise ValueError(f"units of {len(stray)} cycle(s) that are not in the roster of the block: {', '.join(stray[:5])}")
        block_genes = set(pc.unique(units["gene"]).to_pylist())
        self.genes = sorted(block_genes | set(universe or ()))
        self.block_genes = block_genes
        self.u_gene = _vocabulary_codes(units["gene"], self.genes)
        self.u_cycle = _vocabulary_codes(units["cycle"], self.cycles)
        self.units = {name: units[name].to_numpy(zero_copy_only=False) for name in units.column_names
                      if name not in ("gene", "cycle")}
        self.s_gene = _vocabulary_codes(scores["gene"], self.genes)
        self.s_cycle = _vocabulary_codes(scores["cycle"], self.cycles)
        self.s_position = scores["position"].to_numpy(zero_copy_only=False).astype(np.int64)
        self.s_top = scores["top"].to_numpy(zero_copy_only=False).astype(np.float64)
        self.s_bottom = scores["bottom"].to_numpy(zero_copy_only=False).astype(np.float64)


def removal_units(t: BlockTables, lengths: Dict[str, float], mode: str, iqr_multiplier: float,
                  extreme_percentile: float) -> np.ndarray:
    """Boolean per unit: removed by gene removal (core.postproc.gene_removal, each cycle against its own units)."""
    removed = np.zeros(len(t.u_gene), dtype=bool)
    if (mode or "none").lower() == "none" or not len(t.u_gene):
        return removed
    n_pos, clust = t.units["n_pos"], t.units["clust"]
    stream = (GeneUnit(t.cycles[c], t.genes[g], int(n), bool(k))
              for c, g, n, k in zip(t.u_cycle.tolist(), t.u_gene.tolist(), n_pos.tolist(), clust.tolist()))
    flagged = gene_removal(stream, lengths, mode=mode, iqr_multiplier=iqr_multiplier,
                           extreme_percentile=extreme_percentile)
    if flagged:
        cyc_index = {c: i for i, c in enumerate(t.cycles)}
        gene_index = {g: i for i, g in enumerate(t.genes)}
        keys = {cyc_index[c] * len(t.genes) + gene_index[g] for c, g in flagged}
        removed = np.isin(t.u_cycle * len(t.genes) + t.u_gene, np.fromiter(keys, dtype=np.int64, count=len(keys)))
    return removed


def _flagged_pairs(t: BlockTables, removed: np.ndarray) -> List[Tuple[str, str]]:
    return sorted((t.cycles[c], t.genes[g]) for c, g in zip(t.u_cycle[removed].tolist(), t.u_gene[removed].tolist()))


def kept_scores(t: BlockTables, removed: np.ndarray) -> np.ndarray:
    """Boolean per score row: its unit is not removed."""
    if not removed.any():
        return np.ones(len(t.s_gene), dtype=bool)
    keys = (t.u_cycle * len(t.genes) + t.u_gene)[removed]
    return ~np.isin(t.s_cycle * len(t.genes) + t.s_gene, keys)


def reference_pools(t: BlockTables, keep: np.ndarray) -> Dict[int, Dict[str, np.ndarray]]:
    """Per cycle code and direction, the ascending position scores of the kept units (``all``: best side per position)."""
    cyc, top, bottom = t.s_cycle[keep], t.s_top[keep], t.s_bottom[keep]
    best = np.fmax(top, bottom)
    order = np.argsort(cyc, kind="stable")
    cyc, top, bottom, best = cyc[order], top[order], bottom[order], best[order]
    bounds = np.searchsorted(cyc, np.arange(len(t.cycles) + 1), side="left")
    pools: Dict[int, Dict[str, np.ndarray]] = {}
    for c in range(len(t.cycles)):
        lo, hi = bounds[c], bounds[c + 1]
        if lo == hi:
            continue
        tp, bt = top[lo:hi], bottom[lo:hi]
        pools[c] = {"all": np.sort(best[lo:hi]), "top": np.sort(tp[~np.isnan(tp)]), "bottom": np.sort(bt[~np.isnan(bt)])}
    return pools


def gene_statistics(t: BlockTables, removed: np.ndarray, pools: Dict[int, Dict[str, np.ndarray]]):
    """The matrices of one block: {part: {matrix name: genes x cycles array}} for part in max, fmax, size_adj_max, n.

    A cell is 0 where the unit has no value (removed, absent, or no scored position in the direction).
    """
    n_genes, n_cycles = len(t.genes), len(t.cycles)
    out: Dict[str, Dict[str, np.ndarray]] = {p: {} for p in ("max", "fmax", "size_adj_max", "n")}
    keep_unit = ~removed
    for name, direction in DIRECTIONS.items():
        n_d, max_d = t.units[f"n_{direction}"], t.units[f"max_{direction}"]
        m_max = np.zeros((n_genes, n_cycles))
        m_fmax = np.zeros((n_genes, n_cycles))
        m_sam = np.zeros((n_genes, n_cycles))
        m_n = np.zeros((n_genes, n_cycles), dtype=np.int32)
        sel = np.flatnonzero(keep_unit & (n_d > 0))
        for c in np.unique(t.u_cycle[sel]).tolist():
            idx = sel[t.u_cycle[sel] == c]
            pool = pools[c][direction]
            mx = max_d[idx].astype(np.float64)
            f = np.searchsorted(pool, mx + TIE_TOL, side="right") / len(pool)
            n = n_d[idx].astype(np.int64)
            g = t.u_gene[idx]
            m_max[g, c] = mx
            m_fmax[g, c] = f
            m_n[g, c] = n
            m_sam[g, c] = [x ** k for x, k in zip(f.tolist(), n.tolist())]
        out["max"][name], out["fmax"][name], out["size_adj_max"][name], out["n"][name] = m_max, m_fmax, m_sam, m_n
    return out


def write_matrices(directory: Path, genes: Sequence[str], cycles: Sequence[str], parts, gene_stat: str) -> None:
    """The matrices as raw little-endian column-major arrays (``<part>.<name>.f8`` or ``.i4``) and ``manifest.json``.

    A matrix has the genes in rows and the cycles in columns; the manifest lists them and names the matrices that
    form the statistic chosen by ``gene_stat`` (``stat``). assemble_block_matrices.R turns the directory into an RDS.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for part, mats in parts.items():
        for name, mat in mats.items():
            ext = "i4" if part == "n" else "f8"
            fname = f"{part}.{name}.{ext}"
            (directory / fname).write_bytes(mat.astype("<i4" if part == "n" else "<f8").tobytes(order="F"))
            files[f"{part}.{name}"] = fname
    (directory / "manifest.json").write_text(json.dumps({
        "genes": list(genes), "cycles": list(cycles), "gene_stat": gene_stat, "files": files}))


def position_table(t: BlockTables, keep: np.ndarray) -> pa.Table:
    """The scored positions of the kept units: gene, cycle, position and the best-side score ``all``."""
    best = np.fmax(t.s_top[keep], t.s_bottom[keep])
    genes = pa.array(np.asarray(t.genes, dtype=object)[t.s_gene[keep]] if len(t.genes) else [], pa.string())
    cycles = pa.array(np.asarray(t.cycles, dtype=object)[t.s_cycle[keep]] if len(t.cycles) else [], pa.string())
    return pa.table({"gene": genes, "cycle": cycles, "position": pa.array(t.s_position[keep].astype(np.int32)),
                     "all": pa.array(best)})


def position_detections(t: BlockTables, keep: np.ndarray) -> pa.Table:
    """Per scored position of the kept units: gene, position, nd (cycles of the block that score it on either side) and
    ``first`` (the earliest of those cycles in replay order)."""
    n_genes = len(t.genes)
    key = t.s_gene[keep] * (1 << 32) + t.s_position[keep]
    rank_of = {c: r for r, c in enumerate(sorted(t.cycles, key=cycle_order_key))}
    ranks = np.fromiter((rank_of[c] for c in t.cycles), dtype=np.int64, count=len(t.cycles))[t.s_cycle[keep]]
    uniq, inverse = np.unique(key, return_inverse=True)
    nd = np.bincount(inverse, minlength=len(uniq))
    first = np.full(len(uniq), np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(first, inverse, ranks)
    by_rank = sorted(t.cycles, key=cycle_order_key)
    genes = np.asarray(t.genes, dtype=object)[uniq >> 32] if n_genes and len(uniq) else []
    return pa.table({"gene": pa.array(genes, pa.string()), "position": pa.array((uniq & 0xFFFFFFFF).astype(np.int32)),
                     "nd": pa.array(nd.astype(np.int32)),
                     "first": pa.array([by_rank[r] for r in first.tolist()], pa.string())})


def cycle_summaries(t: BlockTables, removed: np.ndarray, keep: np.ndarray):
    """Per cycle and direction: scored entries, genes with an entry and the score histogram (bins of 0.01).

    Returns (summary rows [(cycle, direction, n, n_genes)], histogram rows [(cycle, direction, bin, n)]).
    """
    n_cycles = len(t.cycles)
    keep_unit = ~removed
    summary, hist = [], []
    best = np.fmax(t.s_top[keep], t.s_bottom[keep])
    s_cyc = t.s_cycle[keep]
    for name, direction in DIRECTIONS.items():
        n_d = t.units[f"n_{direction}"]
        n_cycle = np.bincount(t.u_cycle[keep_unit], weights=n_d[keep_unit], minlength=n_cycles)
        genes_cycle = np.bincount(t.u_cycle[keep_unit & (n_d > 0)], minlength=n_cycles)
        values = {"all": best, "top": t.s_top[keep], "bottom": t.s_bottom[keep]}[direction]
        ok = ~np.isnan(values)
        bins = np.minimum(np.floor(values[ok] * 100.0).astype(np.int64), HIST_BINS - 1)
        counts = np.bincount(s_cyc[ok] * HIST_BINS + bins, minlength=n_cycles * HIST_BINS).reshape(n_cycles, HIST_BINS)
        for c, cyc in enumerate(t.cycles):
            summary.append((cyc, name, int(n_cycle[c]), int(genes_cycle[c])))
            for b in np.flatnonzero(counts[c]).tolist():
                hist.append((cyc, name, b, int(counts[c, b])))
    return summary, hist


def sample_and_quantiles(t: BlockTables, removed: np.ndarray, keep: np.ndarray, seed: int, per_cycle: int):
    """Per-cycle uniform sample of the rows of the null and its quantiles.

    The rows of a cycle are ordered by gene and, inside a unit, scored rows first (position order, top before bottom)
    and unscored rows after. Returns (sample rows [(gene, position or 'NA', cycle, score)], quantile records).
    """
    n_genes = len(t.genes)
    keep_unit = np.flatnonzero(~removed)
    unit_cycle, unit_gene = t.u_cycle[keep_unit], t.u_gene[keep_unit]
    order = np.lexsort((unit_gene, unit_cycle))
    keep_unit, unit_cycle, unit_gene = keep_unit[order], unit_cycle[order], unit_gene[order]
    n_rows = t.units["n_rows"][keep_unit].astype(np.int64)
    # scores of the kept units, ordered by (cycle, gene, position)
    s_idx = np.flatnonzero(keep)
    s_order = np.lexsort((t.s_position[s_idx], t.s_gene[s_idx], t.s_cycle[s_idx]))
    s_idx = s_idx[s_order]
    s_key = t.s_cycle[s_idx] * n_genes + t.s_gene[s_idx]
    cycle_bounds = np.searchsorted(unit_cycle, np.arange(len(t.cycles) + 1), side="left")
    sample_rows, quantiles = [], []
    for c, cyc in enumerate(t.cycles):
        lo, hi = cycle_bounds[c], cycle_bounds[c + 1]
        if lo == hi:
            continue
        rows = n_rows[lo:hi]
        total = int(rows.sum())
        if total == 0:
            continue
        cum = np.cumsum(rows)
        rng = random.Random(f"{seed}:{cyc}")
        picks = list(range(total)) if total <= per_cycle else sorted(rng.sample(range(total), per_cycle))
        values = []
        for i in picks:
            u = int(np.searchsorted(cum, i, side="right"))
            j = i - (int(cum[u - 1]) if u else 0)
            unit = keep_unit[lo + u]
            gene = t.genes[t.u_gene[unit]]
            n_scored = int(t.units["n_top"][unit]) + int(t.units["n_bottom"][unit])
            if j < n_scored:
                start = int(np.searchsorted(s_key, c * n_genes + int(t.u_gene[unit]), side="left"))
                stop = start + int(t.units["n_all"][unit])
                top, bottom = t.s_top[s_idx[start:stop]], t.s_bottom[s_idx[start:stop]]
                per_position = (~np.isnan(top)).astype(int) + (~np.isnan(bottom)).astype(int)
                p = int(np.searchsorted(np.cumsum(per_position), j, side="right"))
                before = int(per_position[:p].sum())
                side_values = [v for v in (top[p], bottom[p]) if not np.isnan(v)]
                score, position = float(side_values[j - before]), int(t.s_position[s_idx[start + p]])
            else:
                score, position = 0.0, "NA"
            values.append(score)
            sample_rows.append((gene, position, cyc, score))
        arr = np.asarray(values, dtype=float)
        rec = {"cycle": cyc, "metric": "asr_path_score", "n_sampled": len(arr), "n_total": total, "mean": float(arr.mean())}
        for q, v in zip(QUANT_LEVELS, np.percentile(arr, QUANT_LEVELS)):
            rec[f"q{q}"] = float(v)
        quantiles.append(rec)
    return sample_rows, quantiles


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def merge_block(tile_dirs: Sequence[Path], out_dir: Path, cycles: Sequence[str], *, universe: Optional[Iterable[str]] = None,
                lengths: Optional[Dict[str, float]] = None, gene_filter_mode: str = "none", iqr_multiplier: float = 3.0,
                extreme_percentile: float = 0.99, gene_stat: str = "size_adj_max", seed: int = 1998,
                sample_per_cycle: int = 200) -> dict:
    """Merge the tiles of one block into `out_dir`; returns counts (units, removed units, scored positions, cycles present)."""
    if gene_stat not in GENE_STATS:
        raise ValueError(f"gene_stat must be one of {GENE_STATS}, not {gene_stat!r}")
    import pyarrow.parquet as pq

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t = BlockTables(tile_dirs, cycles, universe)
    removed = removal_units(t, lengths or {}, gene_filter_mode, iqr_multiplier, extreme_percentile)
    keep = kept_scores(t, removed)
    pools = reference_pools(t, keep)
    parts = gene_statistics(t, removed, pools)

    write_matrices(out_dir / "matrices", t.genes, t.cycles, parts, gene_stat)
    # the matrices of the chosen statistic, named for the consumers
    manifest = json.loads((out_dir / "matrices" / "manifest.json").read_text())
    manifest["stat"] = {name: f"{gene_stat}.{name}" for name in DIRECTIONS}
    manifest["n_units"] = int(len(t.u_gene))
    (out_dir / "matrices" / "manifest.json").write_text(json.dumps(manifest))

    pq.write_table(position_table(t, keep), out_dir / "positions.parquet", compression="zstd")
    pq.write_table(position_detections(t, keep), out_dir / "detections.parquet", compression="zstd")
    if (gene_filter_mode or "none").lower() != "none":
        _write_tsv(out_dir / "removed_units.tsv", ("cycle", "Gene"), _flagged_pairs(t, removed))
    summary, hist = cycle_summaries(t, removed, keep)
    _write_tsv(out_dir / "cycle_summary.tsv", ("cycle", "direction", "n", "n_genes"), summary)
    _write_tsv(out_dir / "cycle_hist.tsv", ("cycle", "direction", "bin", "n"), hist)
    sample, quantiles = sample_and_quantiles(t, removed, keep, seed, sample_per_cycle)
    _write_tsv(out_dir / "perm_pos_sample.tsv", ("Gene", "Position", "cycle", "asr_path_score"), sample)
    quant_fields = ["cycle", "metric", "n_sampled", "n_total", "mean"] + [f"q{q}" for q in QUANT_LEVELS]
    _write_tsv(out_dir / "perm_pos_quantiles.tsv", quant_fields, ([rec[f] for f in quant_fields] for rec in quantiles))
    present = sorted({t.cycles[c] for c in np.unique(t.u_cycle[~removed & (t.units["n_rows"] > 0)]).tolist()})
    (out_dir / "cycles.txt").write_text("".join(f"{c}\n" for c in t.cycles))
    counts = {"units": int(len(t.u_gene)), "removed_units": int(removed.sum()), "scored_positions": int(keep.sum()),
              "cycles": len(t.cycles), "cycles_present": present, "score_rule": SCORE_RULE}
    (out_dir / "counts.json").write_text(json.dumps(counts))
    return counts
