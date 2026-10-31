# null_tables.py — Tile tables of the permutation null: one (gene, cycle) unit row and one score row per scored position.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Compact tables of the permutation null, written by the replay of a tile (a batch of genes x one block of cycles) and read
by the block merge (null_merge.py).

Imported by: src/utils/gene_wrapper.py (replay of a gene), src/utils/null_merge.py (readers)
Inputs: the pooled per-position records of one cycle of one gene
Outputs: ``units`` and ``scores`` tables (Parquet, float64 values) and the readers that give them back as numpy arrays

``units``: one row per (gene, cycle) with a detection
    gene, cycle
    n_pos       distinct detected positions of the unit, every row counted (scored or not, trains included): the
                count of core.postproc.gene_removal
    clust       some detected position lies in a cluster train (the dubious-gene test of gene_removal)
    n_rows      records kept once the trains are removed (when they are removed): the rows of the null the position
                sample draws from
    n_top, n_bottom, n_all
                positions with a score > 0 on that side, or on either side (``all``), among the kept rows
    max_top, max_bottom, max_all
                largest score of the direction (NaN when n = 0)
``scores``: one row per (gene, cycle, position) with a score > 0 on some side, among the kept rows
    gene, cycle, position, top, bottom   (NaN where the side has no score > 0)

The records without a score (about half of the rows of the null) are not stored; they are the difference
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

from src.core.postproc import ctrain
from src.core.scores import collapse_sides, position_score

__all__ = ["UNITS_SCHEMA", "SCORES_SCHEMA", "CycleReducer", "CycleResult", "gene_tables", "concat_tables",
           "write_tables", "read_units", "read_scores", "TABLE_NAMES"]

NAN = float("nan")

UNITS_SCHEMA = pa.schema([
    ("gene", pa.string()), ("cycle", pa.string()), ("n_pos", pa.int32()), ("clust", pa.bool_()),
    ("n_rows", pa.int32()), ("n_top", pa.int32()), ("n_bottom", pa.int32()), ("n_all", pa.int32()),
    ("max_top", pa.float64()), ("max_bottom", pa.float64()), ("max_all", pa.float64()),
])
SCORES_SCHEMA = pa.schema([
    ("gene", pa.string()), ("cycle", pa.string()), ("position", pa.int32()),
    ("top", pa.float64()), ("bottom", pa.float64()),
])
# File names of the two tables of a tile (the strict pooling lives in a sibling directory with the same names).
TABLE_NAMES = {"units": "units.parquet", "scores": "scores.parquet"}


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
    positions: List[int]   # scored positions, ascending
    top: List[float]       # per scored position, NaN where the side has no score
    bottom: List[float]


class CycleReducer:
    """Reduces the pooled records of one cycle of one gene to a CycleResult.

    The cluster trains are computed here over the positions of every record of the cycle (core.postproc.ctrain, with
    the gene's position -> column map when there is one); they are flagged only with ``postproc_filter`` and their
    rows are dropped only with ``remove_clusters``, but they count for ``n_pos`` and ``clust`` either way.

    ``universe`` (the gene's positions of the positional universe, None = no restriction) limits what is scored: rows
    outside it are dropped before the scores, the counts and the maxima of the unit are computed, so the gene maximum of
    a cycle is the maximum over universe positions. The detections themselves still count for ``n_pos`` and the trains:
    the dubious-gene test and the cluster filter judge the labeling's detections, as they judge the observed ones.
    """

    def __init__(self, postproc_filter: bool, remove_clusters: bool, clust_minlen: int = 3, clust_maxcaas: float = 0.7,
                 train_columns: Optional[Dict[int, int]] = None, clust_min_positions: Optional[int] = None,
                 universe: Optional[Iterable[int]] = None):
        self.universe = None if universe is None else frozenset(int(p) for p in universe)
        self.postproc_filter = postproc_filter
        self.remove_clusters = remove_clusters
        self.clust_minlen = clust_minlen
        self.clust_min_positions = clust_min_positions
        self.clust_maxcaas = clust_maxcaas
        self.train_columns = train_columns

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
        kept = [row for row in rows if not (self.remove_clusters and row[0] in train)
                and (self.universe is None or row[0] in self.universe)]

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
        scored = sorted(sides)
        top, bottom, allv = [], [], []
        for pos in scored:
            if set(sides[pos]) - {"top", "bottom"}:
                raise ValueError(f"cycle {cycle}: position {pos} has a score on side(s) {sorted(sides[pos])}")
            by_dir = collapse_sides(sides[pos])
            top.append(by_dir.get("top", NAN))
            bottom.append(by_dir.get("bottom", NAN))
            allv.append(by_dir["all"])
        t = [v for v in top if not math.isnan(v)]
        b = [v for v in bottom if not math.isnan(v)]
        return CycleResult(
            cycle, len(positions), bool(train), len(kept), len(t), len(b), len(allv),
            max(t) if t else NAN, max(b) if b else NAN, max(allv) if allv else NAN, scored, top, bottom)


def gene_tables(gene: str, results: List[CycleResult]) -> Tuple[pa.Table, pa.Table]:
    """The units and scores tables of one gene from its per-cycle results (cycles with a detection only)."""
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
    return units, scores


def concat_tables(parts: List[Tuple[pa.Table, pa.Table]]) -> Tuple[pa.Table, pa.Table]:
    """Both tables of several genes, ordered by (gene, cycle) so the files do not depend on the order the genes finished in."""
    if not parts:
        return UNITS_SCHEMA.empty_table(), SCORES_SCHEMA.empty_table()
    keys = [("gene", "ascending"), ("cycle", "ascending")]
    units = pa.concat_tables([u for u, _ in parts]).sort_by(keys)
    scores = pa.concat_tables([s for _, s in parts]).sort_by(keys + [("position", "ascending")])
    return units, scores


def write_tables(directory: Path, units: pa.Table, scores: pa.Table, meta: Optional[Dict[str, Any]] = None) -> None:
    """Write ``units.parquet`` and ``scores.parquet`` (zstd) under `directory`, with `meta` as ``meta.json``."""
    from src.core.scores import SCORE_RULE
    record = dict(meta or {})
    if record.get('score_rule', SCORE_RULE) != SCORE_RULE:
        raise ValueError('Cannot write a tile with a retired score rule')
    record['score_rule'] = SCORE_RULE
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    pq.write_table(units, directory / TABLE_NAMES["units"], compression="zstd")
    pq.write_table(scores, directory / TABLE_NAMES["scores"], compression="zstd")
    (directory / "meta.json").write_text(json.dumps(record, sort_keys=True))


def validate_tile(path):
    from src.core.scores import SCORE_RULE
    directory = path if path.is_dir() else path.parent
    stamp = directory / 'meta.json'
    if not stamp.is_file() or json.loads(stamp.read_text()).get('score_rule') != SCORE_RULE:
        raise ValueError(f'Unverifiable/retired CAAS tile score rule: {directory}; regenerate observed and null artifacts together')


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
