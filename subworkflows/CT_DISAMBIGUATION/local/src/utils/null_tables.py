# null_tables.py — Tile tables of the permutation null: one (gene, cycle) unit row and one score row per scored position.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Compact tables of the permutation null, written by the replay of a tile (a batch of genes x one block of cycles) and read
by the block merge (null_merge.py).

Imported by: src/utils/gene_wrapper.py (replay of a gene), src/utils/null_merge.py (readers)
Inputs: the pooled per-position records of one cycle of one gene
Outputs: ``units``, ``scores`` and ``scores_diag`` tables (Parquet) and the readers that give them back as numpy arrays

Two levels of the same unit are kept. The position level is limited to the positional universe U of the gene (None = no
restriction): it feeds the position p-values and everything that conditions on U. The gene level covers every record
that survives the trains, whatever its position: it feeds the gene statistic of each cycle, so that the null and the
observed labeling are summarized over the same kind of positions.

``units``: one row per (gene, cycle) with a detection
    gene, cycle
    n_pos       distinct detected positions of the unit, every row counted (scored or not, trains included): the
                count of core.gene_rules.gene_removal
    clust       some detected position lies in a cluster train (the dubious-gene test of gene_removal)
    n_rows      records kept once the trains are removed (when they are removed) and inside U: the rows of the null the
                position sample draws from
    n_top, n_bottom, n_all
                positions inside U with a score > 0 on that side, or on either side (``all``), among those rows
    max_top, max_bottom, max_all
                largest score of the direction among them (NaN when n = 0)
    gn_rows, gn_top, gn_bottom, gn_all, gmax_top, gmax_bottom, gmax_all
                the same quantities over every record kept once the trains are removed, inside U or not
``scores``: one row per (gene, cycle, position inside U) with a score > 0 on some side, among the kept rows
    gene, cycle, position, top, bottom   (NaN where the side has no score > 0)
``scores_diag``: one row per (gene, cycle, position outside U) with a score > 0 on some side, for the leading cycles
    only (``diag_cycles``); the score is reduced to its bin of width 0.01 (floor(100 * score), at most 100)
    gene, cycle, top, bottom   (uint8, 255 where the side has no score > 0)

The records without a score (about half of the rows of the null) are not stored; inside U they are the difference
``n_rows - n_top - n_bottom`` of their unit. Gene removal is a decision on a whole (gene, cycle) unit and the trains are
decided per cycle, so every value of a unit is final when the tile is written; the only inputs the merge adds are the
removal thresholds (which need every gene of the cycle) and the reference pools of the cycle.

The strict pooling of the same replay is a second pair of tables of the same layout.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Tuple

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from src.core.trains import ctrain
from src.core.scores import collapse_sides, cycle_order_key, position_score

__all__ = ["UNITS_SCHEMA", "SCORES_SCHEMA", "DIAG_SCHEMA", "CycleReducer", "CycleResult", "gene_tables", "concat_tables",
           "write_tables", "read_units", "read_scores", "read_diag", "TABLE_NAMES", "TILE_LAYOUT", "DIAG_NA", "score_bin",
           "diag_cycle"]

NAN = float("nan")

# Layout of the tile tables, stamped in meta.json: a tile of another layout is refused by the readers.
TILE_LAYOUT = 2

UNITS_SCHEMA = pa.schema([
    ("gene", pa.string()), ("cycle", pa.string()), ("n_pos", pa.int32()), ("clust", pa.bool_()),
    ("n_rows", pa.int32()), ("n_top", pa.int32()), ("n_bottom", pa.int32()), ("n_all", pa.int32()),
    ("max_top", pa.float64()), ("max_bottom", pa.float64()), ("max_all", pa.float64()),
    ("gn_rows", pa.int32()), ("gn_top", pa.int32()), ("gn_bottom", pa.int32()), ("gn_all", pa.int32()),
    ("gmax_top", pa.float64()), ("gmax_bottom", pa.float64()), ("gmax_all", pa.float64()),
])
SCORES_SCHEMA = pa.schema([
    ("gene", pa.string()), ("cycle", pa.string()), ("position", pa.int32()),
    ("top", pa.float64()), ("bottom", pa.float64()),
])
# Bin of a score outside U: floor(100 * score), at most 100; DIAG_NA where the side has no score.
DIAG_NA = 255
DIAG_SCHEMA = pa.schema([
    ("gene", pa.string()), ("cycle", pa.string()), ("top", pa.uint8()), ("bottom", pa.uint8()),
])
# File names of the tables of a tile (the strict pooling lives in a sibling directory with the same names).
TABLE_NAMES = {"units": "units.parquet", "scores": "scores.parquet", "diag": "scores_diag.parquet"}


def score_bin(value: float) -> int:
    """Bin of width 0.01 of a score (the histogram of the null report): floor(100 * score), the score 1.0 in bin 100."""
    return DIAG_NA if math.isnan(value) else min(int(math.floor(value * 100.0)), 100)


def diag_cycle(cycle: str, diag_cycles: int) -> bool:
    """Whether the cycle keeps the bins of its scores outside U: b_0 and b_1..b_<diag_cycles> (0 = none)."""
    if diag_cycles <= 0:
        return False
    key = cycle_order_key(cycle.split("~", 1)[0])
    return key[0] == 0 and key[1] <= diag_cycles


class CycleResult(NamedTuple):
    """One (gene, cycle) unit and its scored positions, reduced from the pooled records of the cycle."""
    cycle: str
    n_pos: int
    clust: bool
    n_rows: int
    n_top: int
    n_bottom: int
    n_all: int
    max_top: float
    max_bottom: float
    max_all: float
    positions: List[int]   # scored positions inside U, ascending
    top: List[float]       # per scored position, NaN where the side has no score
    bottom: List[float]
    gn_rows: int           # gene level: every record kept once the trains are removed, inside U or not
    gn_top: int
    gn_bottom: int
    gn_all: int
    gmax_top: float
    gmax_bottom: float
    gmax_all: float
    diag_top: List[int]    # bins of the positions scored outside U (empty when the cycle keeps none)
    diag_bottom: List[int]


class CycleReducer:
    """Reduces the pooled records of one cycle of one gene to a CycleResult.

    The cluster trains are computed here over the positions of every record of the cycle (core.trains.ctrain, with
    the gene's position -> column map when there is one); they are flagged only with ``postproc_filter`` and their
    rows are dropped only with ``remove_clusters``, but they count for ``n_pos`` and ``clust`` either way.

    Two levels are summarized from the records that survive the trains. The gene level (``gn_*``, ``gmax_*``) takes all
    of them: the gene maximum of a labeling is the largest score among the positions it detects, wherever they lie.
    The position level (``n_*``, ``max_*``, the scored positions) takes only the gene's positions of the positional
    universe (``universe``, None = no restriction): it is the set the position p-values condition on. With
    ``diag_cycles`` > 0 the cycles b_0 and b_1..b_<diag_cycles> also keep the bins of the scores outside the universe.
    The detections count for ``n_pos`` and the trains whatever the universe: the dubious-gene test and the cluster filter
    judge the labeling's detections, as they judge the observed ones.
    """

    def __init__(self, postproc_filter: bool, remove_clusters: bool, clust_minlen: int = 3, clust_maxcaas: float = 0.7,
                 train_columns: Optional[Dict[int, int]] = None, clust_min_positions: Optional[int] = None,
                 universe: Optional[Iterable[int]] = None, diag_cycles: int = 0):
        self.universe = None if universe is None else frozenset(int(p) for p in universe)
        self.postproc_filter = postproc_filter
        self.remove_clusters = remove_clusters
        self.clust_minlen = clust_minlen
        self.clust_min_positions = clust_min_positions
        self.clust_maxcaas = clust_maxcaas
        self.train_columns = train_columns
        self.diag_cycles = int(diag_cycles)

    def _inside(self, position: int) -> bool:
        return self.universe is None or position in self.universe

    def __call__(self, cycle: str, records: Iterable[Any]) -> Optional[CycleResult]:
        """The unit of the cycle, or None when the cycle has no record with a position (it leaves no unit)."""
        rows: List[Tuple[int, float, str]] = []
        for r in records:
            pos = getattr(r, "position", None)
            if pos is None:
                continue
            score = getattr(r, "asr_path_score", 0.0)
            rows.append((int(pos), 0.0 if score is None else float(score), getattr(r, "side", "none") or "none"))
        if not rows:
            return None
        positions = {p for p, _, _ in rows}
        train = set()
        if self.postproc_filter:
            train = set(ctrain(list(positions), self.clust_maxcaas, self.clust_minlen, self.train_columns, self.clust_min_positions))
        kept = [row for row in rows if not (self.remove_clusters and row[0] in train)]

        # score of each (position, side) (the last record wins, as a record key is unique by construction), then per
        # position the sides that score
        by_key: Dict[Tuple[int, str], float] = {}
        for pos, score, side in kept:
            by_key[(pos, side)] = score
        sides: Dict[int, Dict[str, float]] = {}
        for (pos, side), score in by_key.items():
            value = position_score(score)
            if value is not None:
                sides.setdefault(pos, {})[side] = value
        # (top, bottom, best side) of every scored position, whatever its place
        values: Dict[int, Tuple[float, float, float]] = {}
        for pos in sorted(sides):
            if set(sides[pos]) - {"top", "bottom"}:
                raise ValueError(f"cycle {cycle}: position {pos} has a score on side(s) {sorted(sides[pos])}")
            by_dir = collapse_sides(sides[pos])
            values[pos] = (by_dir.get("top", NAN), by_dir.get("bottom", NAN), by_dir["all"])

        def summary(selected: List[int]):
            t = [values[p][0] for p in selected if not math.isnan(values[p][0])]
            b = [values[p][1] for p in selected if not math.isnan(values[p][1])]
            a = [values[p][2] for p in selected]
            return (len(t), len(b), len(a), max(t) if t else NAN, max(b) if b else NAN, max(a) if a else NAN)

        inside = [p for p in values if self._inside(p)]
        n_t, n_b, n_a, m_t, m_b, m_a = summary(inside)
        g_t, g_b, g_a, gm_t, gm_b, gm_a = summary(list(values))
        diag_top: List[int] = []
        diag_bottom: List[int] = []
        if self.universe is not None and diag_cycle(cycle, self.diag_cycles):
            for pos in values:
                if not self._inside(pos):
                    diag_top.append(score_bin(values[pos][0]))
                    diag_bottom.append(score_bin(values[pos][1]))
        return CycleResult(
            cycle, len(positions), bool(train), sum(1 for row in kept if self._inside(row[0])), n_t, n_b, n_a, m_t, m_b, m_a,
            inside, [values[p][0] for p in inside], [values[p][1] for p in inside],
            len(kept), g_t, g_b, g_a, gm_t, gm_b, gm_a, diag_top, diag_bottom)


def gene_tables(gene: str, results: List[CycleResult]) -> Tuple[pa.Table, pa.Table, pa.Table]:
    """The units, scores and scores_diag tables of one gene from its per-cycle results (cycles with a detection only)."""
    n = len(results)
    units = pa.table({
        "gene": pa.array([gene] * n, pa.string()),
        "cycle": pa.array([r.cycle for r in results], pa.string()),
        "n_pos": pa.array([r.n_pos for r in results], pa.int32()),
        "clust": pa.array([r.clust for r in results], pa.bool_()),
        "n_rows": pa.array([r.n_rows for r in results], pa.int32()),
        "n_top": pa.array([r.n_top for r in results], pa.int32()),
        "n_bottom": pa.array([r.n_bottom for r in results], pa.int32()),
        "n_all": pa.array([r.n_all for r in results], pa.int32()),
        "max_top": pa.array([r.max_top for r in results], pa.float64()),
        "max_bottom": pa.array([r.max_bottom for r in results], pa.float64()),
        "max_all": pa.array([r.max_all for r in results], pa.float64()),
        "gn_rows": pa.array([r.gn_rows for r in results], pa.int32()),
        "gn_top": pa.array([r.gn_top for r in results], pa.int32()),
        "gn_bottom": pa.array([r.gn_bottom for r in results], pa.int32()),
        "gn_all": pa.array([r.gn_all for r in results], pa.int32()),
        "gmax_top": pa.array([r.gmax_top for r in results], pa.float64()),
        "gmax_bottom": pa.array([r.gmax_bottom for r in results], pa.float64()),
        "gmax_all": pa.array([r.gmax_all for r in results], pa.float64()),
    }, schema=UNITS_SCHEMA)
    counts = [len(r.positions) for r in results]
    total = sum(counts)
    scores = pa.table({
        "gene": pa.array([gene] * total, pa.string()),
        "cycle": pa.array([r.cycle for r, c in zip(results, counts) for _ in range(c)], pa.string()),
        "position": pa.array([p for r in results for p in r.positions], pa.int32()),
        "top": pa.array([v for r in results for v in r.top], pa.float64()),
        "bottom": pa.array([v for r in results for v in r.bottom], pa.float64()),
    }, schema=SCORES_SCHEMA)
    outside = [len(r.diag_top) for r in results]
    diag = pa.table({
        "gene": pa.array([gene] * sum(outside), pa.string()),
        "cycle": pa.array([r.cycle for r, c in zip(results, outside) for _ in range(c)], pa.string()),
        "top": pa.array([v for r in results for v in r.diag_top], pa.uint8()),
        "bottom": pa.array([v for r in results for v in r.diag_bottom], pa.uint8()),
    }, schema=DIAG_SCHEMA)
    return units, scores, diag


def concat_tables(parts: List[Tuple[pa.Table, pa.Table, pa.Table]]) -> Tuple[pa.Table, pa.Table, pa.Table]:
    """The three tables of several genes, ordered by (gene, cycle) so the files do not depend on the order the genes finished in."""
    if not parts:
        return UNITS_SCHEMA.empty_table(), SCORES_SCHEMA.empty_table(), DIAG_SCHEMA.empty_table()
    keys = [("gene", "ascending"), ("cycle", "ascending")]
    units = pa.concat_tables([part[0] for part in parts]).sort_by(keys)
    scores = pa.concat_tables([part[1] for part in parts]).sort_by(keys + [("position", "ascending")])
    diag = pa.concat_tables([part[2] for part in parts]).sort_by(keys)
    return units, scores, diag


def write_tables(directory: Path, units: pa.Table, scores: pa.Table, diag: Optional[pa.Table] = None,
                 meta: Optional[Dict[str, Any]] = None) -> None:
    """Write ``units.parquet``, ``scores.parquet`` and ``scores_diag.parquet`` (zstd) under `directory`, with `meta` as ``meta.json``."""
    from src.core.scores import SCORE_RULE
    record = dict(meta or {})
    if record.get('score_rule', SCORE_RULE) != SCORE_RULE:
        raise ValueError('Cannot write a tile with a retired score rule')
    record['score_rule'] = SCORE_RULE
    record['layout'] = TILE_LAYOUT
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    pq.write_table(units, directory / TABLE_NAMES["units"], compression="zstd")
    pq.write_table(scores, directory / TABLE_NAMES["scores"], compression="zstd")
    pq.write_table(diag if diag is not None else DIAG_SCHEMA.empty_table(), directory / TABLE_NAMES["diag"], compression="zstd")
    (directory / "meta.json").write_text(json.dumps(record, sort_keys=True))


def validate_tile(path):
    from src.core.scores import SCORE_RULE
    directory = path if path.is_dir() else path.parent
    stamp = directory / 'meta.json'
    meta = json.loads(stamp.read_text()) if stamp.is_file() else {}
    if meta.get('score_rule') != SCORE_RULE:
        raise ValueError(f'Unverifiable/retired CAAS tile score rule: {directory}; regenerate observed and null artifacts together')
    if meta.get('layout') != TILE_LAYOUT:
        raise ValueError(f'CAAS tile of another layout (found {meta.get("layout")}, need {TILE_LAYOUT}): {directory}; '
                         'regenerate observed and null artifacts together')


def read_units(path: Path) -> pa.Table:
    """The units table of a tile directory (or of a single parquet file)."""
    path = Path(path)
    validate_tile(path)
    return pq.read_table(path / TABLE_NAMES["units"] if path.is_dir() else path, schema=UNITS_SCHEMA)


def read_scores(path: Path) -> pa.Table:
    """The scores table of a tile directory (or of a single parquet file)."""
    path = Path(path)
    validate_tile(path)
    return pq.read_table(path / TABLE_NAMES["scores"] if path.is_dir() else path, schema=SCORES_SCHEMA)


def read_diag(path: Path) -> pa.Table:
    """The scores_diag table of a tile directory (or of a single parquet file)."""
    path = Path(path)
    validate_tile(path)
    return pq.read_table(path / TABLE_NAMES["diag"] if path.is_dir() else path, schema=DIAG_SCHEMA)
