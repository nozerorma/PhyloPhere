"""The scripts that turn the tiles of a block into the null: merge_null_block.py, assemble_block_matrices.R,
reduce_null_blocks.py and scoring_caas_perms.R, end to end on synthetic tiles.

The roster is the one of the labelings file (every replayed cycle, the ones with no detection included, the real labeling
left out), N is the roster's size, and caas_perms.rds holds one column per roster cycle and the stamp of the gene statistic.
"""
import json
import shutil
import subprocess
import sys

import pyarrow.parquet as pq
import pytest

from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1]

LOCAL = PIPELINE / "subworkflows/CT_DISAMBIGUATION/local"
SCORING = PIPELINE / "subworkflows/SCORING/local/src"
sys.path.insert(0, str(LOCAL))

from src.utils.null_tables import CycleReducer, concat_tables, gene_tables, write_tables  # noqa: E402
from src.utils.perm_stream import SlimRecord  # noqa: E402

RSCRIPT = shutil.which("Rscript")
pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="Rscript not available")


def _tile(directory, units, by_labeling=False):
    """units: {(gene, cycle): [(position, score, side)]}."""
    parts = {}
    for (gene, cycle), recs in units.items():
        parts.setdefault(gene, []).append(CycleReducer(False, True, 3, 0.7, None)(cycle, [SlimRecord(p, s, side, 1, 1) for p, s, side in recs]))
    write_tables(directory, *concat_tables([gene_tables(g, r) for g, r in parts.items()]), meta={"by_labeling": by_labeling})
    return directory


def _labelings(path, tags):
    path.write_text("".join(f"{t}\tsp1,sp2\tsp3,sp4\n" for t in tags))
    return path


def _merge(tmp_path, tiles, labelings, name, *extra):
    out = tmp_path / name
    r = subprocess.run([sys.executable, str(LOCAL / "merge_null_block.py"), "--tiles", *map(str, tiles), "--labelings", str(labelings),
                        "--output-dir", str(out), *extra], cwd=LOCAL, capture_output=True, text=True)
    return r, out


def _rds(out):
    r = subprocess.run([RSCRIPT, str(SCORING / "assemble_block_matrices.R"), "--dir", str(out / "matrices"), "--output", str(out / "matrices.rds")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    shutil.rmtree(out / "matrices")
    return out


def _caas_perms(tmp_path, blocks, *extra):
    r = subprocess.run([RSCRIPT, str(SCORING / "scoring_caas_perms.R"), "--blocks", *[str(b / "matrices.rds") for b in blocks], *extra,
                        "--output", str(tmp_path / "caas_perms.rds")], capture_output=True, text=True)
    return r, tmp_path / "caas_perms.rds"


def _r(expr):
    r = subprocess.run([RSCRIPT, "-e", expr], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    return r.stdout


def test_the_roster_is_the_labelings_without_the_real_one_and_keeps_empty_cycles(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_1"): [(3, 0.5, "top")]})
    lab = _labelings(tmp_path / "lab.tab", ["b_1", "b_2", "b_3"])
    r, out = _merge(tmp_path, [tile], lab, "blk")
    assert r.returncode == 0, r.stderr[-800:]
    assert (out / "cycles.txt").read_text().split() == ["b_1", "b_2", "b_3"]
    assert json.loads((out / "counts.json").read_text())["cycles_present"] == ["b_1"]
    assert "2 of 3 cycles of the block left no row" in r.stderr


def test_the_real_labeling_is_a_block_of_its_own(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_0"): [(3, 0.5, "top")]})
    lab = _labelings(tmp_path / "lab.tab", ["b_0"])
    r, out = _merge(tmp_path, [tile], lab, "b0", "--include-b0")
    assert r.returncode == 0, r.stderr[-800:]
    assert (out / "cycles.txt").read_text().split() == ["b_0"]
    assert _merge(tmp_path, [tile], lab, "bad")[0].returncode != 0           # b_0 is no cycle of a block of permuted cycles


def test_the_roster_keeps_the_grain_of_the_tiles_when_hypotheses_were_not_pooled(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_1~H1"): [(3, 0.5, "top")]}, by_labeling=True)
    lab = _labelings(tmp_path / "lab.tab", ["b_1~H1", "b_1~H2"])
    r, out = _merge(tmp_path, [tile], lab, "blk")
    assert r.returncode == 0, r.stderr[-800:]
    assert (out / "cycles.txt").read_text().split() == ["b_1~H1", "b_1~H2"]


def test_a_cycle_with_units_that_is_not_in_the_labelings_is_an_error(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_9"): [(3, 0.5, "top")]})
    r, _ = _merge(tmp_path, [tile], _labelings(tmp_path / "lab.tab", ["b_1", "b_2"]), "blk")
    assert r.returncode != 0 and "not in the roster" in r.stderr


def test_the_rds_has_a_column_for_every_cycle_of_the_roster_and_the_gene_statistic_stamp(tmp_path):
    t1 = _tile(tmp_path / "t1", {("G1", "b_1"): [(3, 0.5, "top"), (4, 0.7, "bottom")], ("G2", "b_1"): [(3, 0.4, "top")]})
    t2 = _tile(tmp_path / "t2", {("G1", "b_3"): [(3, 0.6, "top")]})
    r1, blk1 = _merge(tmp_path, [t1], _labelings(tmp_path / "l1.tab", ["b_1", "b_2"]), "blk1")
    r2, blk2 = _merge(tmp_path, [t2], _labelings(tmp_path / "l2.tab", ["b_3", "b_4"]), "blk2")
    assert r1.returncode == 0 and r2.returncode == 0, r1.stderr[-400:] + r2.stderr[-400:]
    universe = tmp_path / "universe.txt"
    universe.write_text("Gene\nG1\nG2\nG3\n")
    roster = tmp_path / "roster.txt"
    roster.write_text("b_1\nb_2\nb_3\nb_4\n")
    r, rds = _caas_perms(tmp_path, [_rds(blk1), _rds(blk2)], "--cycles", str(roster), "--universe", str(universe))
    assert r.returncode == 0, r.stderr[-800:]
    out = _r(f'x <- readRDS("{rds}"); cat(x$gene_stat, names(x$caas_corStat_byrank), colnames(x$caas_corStat_byrank$global), rownames(x$caas_corStat_byrank$global), sum(x$caas_corStat_byrank$global > 0))')
    assert out.split() == ["max", "global", "top", "bottom", "b_1", "b_2", "b_3", "b_4", "G1", "G2", "G3", "3"]


def test_a_roster_that_misses_a_cycle_with_units_stops_the_build(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_1"): [(3, 0.5, "top")]})
    r, blk = _merge(tmp_path, [tile], _labelings(tmp_path / "l.tab", ["b_1"]), "blk")
    assert r.returncode == 0
    roster = tmp_path / "roster.txt"
    roster.write_text("b_2\n")
    r, _ = _caas_perms(tmp_path, [_rds(blk)], "--cycles", str(roster))
    assert r.returncode != 0 and "not in" in r.stderr


@pytest.mark.parametrize("stat", ["size_adj_max", "fmax"])
def test_retired_statistics_are_rejected_before_merging(tmp_path, stat):
    tile = _tile(tmp_path / "t", {("G1", "b_1"): [(3, 0.5, "top")]})
    result, _ = _merge(tmp_path, [tile], _labelings(tmp_path / "l.tab", ["b_1"]), "blk", "--gene-stat", stat)
    assert result.returncode != 0 and "invalid choice" in result.stderr


def test_stat_is_raw_max_and_transformed_parts_are_not_written(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_1"): [(3, 0.5, "top"), (4, 0.25, "top")], ("G2", "b_1"): [(3, 0.5, "top")]})
    result, blk = _merge(tmp_path, [tile], _labelings(tmp_path / "l.tab", ["b_1"]), "blk")
    assert result.returncode == 0, result.stderr[-800:]
    rds = _rds(blk) / "matrices.rds"
    out = _r(f'x <- readRDS("{rds}"); cat(x$stat$top["G1", "b_1"], x$parts$n$top["G1", "b_1"], names(x$parts))')
    assert out.split() == ["0.5", "2", "max", "n"]


def test_an_null_with_no_unit_gives_an_empty_rds(tmp_path):
    tile = _tile(tmp_path / "t", {})
    r, blk = _merge(tmp_path, [tile], _labelings(tmp_path / "l.tab", ["b_1", "b_2"]), "blk")
    assert r.returncode == 0, r.stderr[-800:]
    r, rds = _caas_perms(tmp_path, [_rds(blk)])
    assert r.returncode == 0 and "no scored rows" in (r.stdout + r.stderr)
    assert _r(f'x <- readRDS("{rds}"); cat(length(x$caas_corStat_byrank), x$gene_stat)').split() == ["0", "max"]


def _reduce(tmp_path, blocks):
    out = tmp_path / "reduced"
    r = subprocess.run([sys.executable, str(LOCAL / "reduce_null_blocks.py"), "--blocks", *map(str, blocks), "--output-dir", str(out)],
                       capture_output=True, text=True)
    return r, out


def test_the_reduce_joins_the_roster_and_the_tables_of_the_blocks(tmp_path):
    t1 = _tile(tmp_path / "t1", {("G1", "b_1"): [(3, 0.5, "top")], ("G2", "b_1"): [(5, 0.4, "bottom")]})
    t2 = _tile(tmp_path / "t2", {("G1", "b_10"): [(3, 0.6, "top")]})
    _, a = _merge(tmp_path, [t1], _labelings(tmp_path / "l1.tab", ["b_1", "b_2"]), "a")
    _, b = _merge(tmp_path, [t2], _labelings(tmp_path / "l2.tab", ["b_9", "b_10"]), "b")
    r, out = _reduce(tmp_path, [a, b])
    assert r.returncode == 0, r.stderr[-800:]
    assert (out / "cycle_roster.txt").read_text().split() == ["b_1", "b_10", "b_2", "b_9"]            # text order, as the roster of the null always was
    pos = {(g, p): (nd, first) for g, p, nd, first in
           ((r_["gene"], r_["position"], r_["nd"], r_["first_rank"]) for r_ in pq.read_table(out / "null_positions.parquet").to_pylist())}
    # replay order b_1, b_2, b_9, b_10: G1 position 3 is scored in b_1 and b_10, first in rank 0
    assert pos == {("G1", 3): (2, 0), ("G2", 5): (1, 0)}
    assert json.loads((out / "counts.json").read_text())["cycles"] == 4
    assert (out / "cycle_summary.tsv").read_text().count("\n") == 1 + 4 * 3


def test_the_reduce_without_the_tables_of_a_block_fails(tmp_path):
    tile = _tile(tmp_path / "t", {("G1", "b_1"): [(3, 0.5, "top")]})
    _, a = _merge(tmp_path, [tile], _labelings(tmp_path / "l.tab", ["b_1"]), "a")
    (a / "perm_pos_sample.tsv").unlink()
    r, _ = _reduce(tmp_path, [a])
    assert r.returncode != 0


@pytest.mark.parametrize('fault', ['same_block', 'overlap_empty_cycle', 'repeated_cycle', 'wrong_count'])
def test_reduce_refuses_inconsistent_cycle_ownership(tmp_path, fault):
    tile = _tile(tmp_path / 't', {('G1', 'b_1'): [(3, 0.5, 'top')]})
    r, a = _merge(tmp_path, [tile], _labelings(tmp_path / 'la.tab', ['b_1', 'b_2']), 'a')
    assert r.returncode == 0, r.stderr
    blocks = [a]
    if fault == 'same_block':
        blocks.append(a)
    elif fault == 'overlap_empty_cycle':
        tile_b = _tile(tmp_path / 'tb', {('G1', 'b_3'): [(3, 0.6, 'top')]})
        r, b = _merge(tmp_path, [tile_b], _labelings(tmp_path / 'lb.tab', ['b_2', 'b_3']), 'b')
        assert r.returncode == 0, r.stderr
        blocks.append(b)
    elif fault == 'repeated_cycle':
        with (a / 'cycles.txt').open('a') as f:
            f.write('b_2\n')
    else:
        meta = json.loads((a / 'counts.json').read_text())
        meta['cycles'] += 1
        (a / 'counts.json').write_text(json.dumps(meta))
    r, out = _reduce(tmp_path, blocks)
    assert r.returncode != 0, 'inconsistent cycle ownership was accepted'
    assert 'cycle' in r.stderr.lower()
    assert not (out / 'cycle_roster.txt').exists(), 'an invalid roster was published'


@pytest.mark.parametrize('fault', ['same_block', 'overlap_empty_cycle'])
def test_rds_assembly_refuses_overlapping_cycles(tmp_path, fault):
    tile = _tile(tmp_path / 't', {('G1', 'b_1'): [(3, 0.5, 'top')]})
    r, a = _merge(tmp_path, [tile], _labelings(tmp_path / 'la.tab', ['b_1', 'b_2']), 'a')
    assert r.returncode == 0, r.stderr
    a = _rds(a)
    b = a
    if fault == 'overlap_empty_cycle':
        tile_b = _tile(tmp_path / 'tb', {('G1', 'b_3'): [(3, 0.6, 'top')]})
        r, b = _merge(tmp_path, [tile_b], _labelings(tmp_path / 'lb.tab', ['b_2', 'b_3']), 'b')
        assert r.returncode == 0, r.stderr
        b = _rds(b)
    r, out = _caas_perms(tmp_path, [a, b])
    assert r.returncode != 0, 'overlapping matrix columns were silently overwritten'
    assert 'cycle' in r.stderr.lower()
    assert not out.exists()
