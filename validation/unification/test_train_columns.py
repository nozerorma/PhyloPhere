"""Cluster trains measured in untrimmed alignment columns: core, observed filter and null agree, and off changes nothing."""
import random
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.columns import column_map, find_file, gene_columns, index_files, read_map  # noqa: E402
from src.core.postproc import ctrain, train_flags  # noqa: E402
from src.utils.gene_wrapper import _gene_train_columns, _perms_worker_finalize  # noqa: E402

FILTER = ROOT / "subworkflows/CT_POSTPROC/local/src/filter_caas_clusters-param.py"


def _legacy_ctrain(positions, maxcaas=0.7, minlen=3):
    """ctrain as it was before the optional column map."""
    uniq = sorted({int(p) for p in positions})
    n = len(uniq)
    if n < minlen:
        return []
    bad = set()
    for r in range(n):
        for l in range(r + 1):
            span = uniq[r] - uniq[l] + 1
            if span >= minlen and (r - l + 1) / span >= maxcaas:
                bad.update(uniq[l:r + 1])
    return sorted(bad)


def _write_map(path, selected_ori, n_ori):
    """MAP table whose selected untrimmed columns are `selected_ori` (increasing) out of 1..n_ori."""
    sel = {c: i + 1 for i, c in enumerate(selected_ori)}
    with open(path, "w") as fh:
        fh.write("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n")
        for c in range(1, n_ori + 1):
            fh.write(f"{c}\t{'selected' if c in sel else 'removed'}\t{sel.get(c, 'NA')}\t{sel.get(c, 'NA')}\n")


# trimmed positions 9, 10, 11 (prot columns 10..12) sit at untrimmed columns 10, 14, 20
SPREAD = list(range(1, 10)) + [10, 14, 20] + list(range(21, 60))


# ── core.columns ─────────────────────────────────────────────────────────────

def test_column_map_gives_the_untrimmed_column_of_each_zero_based_position(tmp_path):
    _write_map(tmp_path / "G.map.tsv", [2, 3, 6], 7)
    removed, ori = read_map(tmp_path / "G.map.tsv")
    assert removed == [True, False, False, True, True, False, True] and ori == {1: 2, 2: 3, 3: 6}
    assert column_map(tmp_path / "G.map.tsv") == {0: 2, 1: 3, 2: 6}


def test_read_map_rejects_unknown_status_and_broken_numbering(tmp_path):
    (tmp_path / "a.tsv").write_text("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n1\tkept\t1\t1\n")
    with pytest.raises(ValueError, match="status"):
        read_map(tmp_path / "a.tsv")
    (tmp_path / "b.tsv").write_text("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n1\tselected\t1\t2\n")
    with pytest.raises(ValueError, match="prot_ali_col"):
        read_map(tmp_path / "b.tsv")
    (tmp_path / "c.tsv").write_text("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n2\tselected\t1\t1\n")
    with pytest.raises(ValueError, match="ori_codon_col"):
        read_map(tmp_path / "c.tsv")


def test_a_gene_is_found_whatever_its_species_and_a_missing_or_ambiguous_one_is_reported(tmp_path):
    for n in ("A.Homo_sapiens.map.tsv", "B.Lemur_catta.map.tsv", "C.1.Homo_sapiens.map.tsv", "D.Homo_sapiens.map.tsv",
              "D.Papio_anubis.map.tsv"):
        _write_map(tmp_path / n, [1, 2], 2)
    idx = index_files(tmp_path, ".map.tsv")
    assert set(idx) == {"A", "B", "C", "D"} and idx["D"] is None
    assert gene_columns(idx, "B", ".map.tsv") == {0: 1, 1: 2}
    assert gene_columns(idx, "Z", ".map.tsv") is None
    with pytest.raises(ValueError, match="several files"):
        gene_columns(idx, "D", ".map.tsv")
    with pytest.raises(FileNotFoundError):
        find_file(idx, "Z", ".map.tsv")


def test_an_inconsistent_map_file_raises_instead_of_being_used(tmp_path):
    (tmp_path / "G.map.tsv").write_text("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n1\tmaybe\t1\t1\n")
    with pytest.raises(ValueError):
        gene_columns(index_files(tmp_path, ".map.tsv"), "G", ".map.tsv")


# ── ctrain ───────────────────────────────────────────────────────────────────

def test_without_columns_ctrain_is_what_it_was():
    rng = random.Random(3)
    for _ in range(300):
        pos = rng.sample(range(80), rng.randint(0, 30))
        maxcaas, minlen = rng.choice([0.5, 0.6, 0.7, 0.9]), rng.choice([1, 2, 3, 4, 10])
        assert ctrain(pos, maxcaas, minlen) == _legacy_ctrain(pos, maxcaas, minlen)


def test_a_constant_shift_of_the_columns_changes_nothing():
    rng = random.Random(5)
    for _ in range(200):
        pos = rng.sample(range(80), rng.randint(0, 30))
        assert ctrain(pos, 0.7, 3, {p: p + 1 for p in range(80)}) == ctrain(pos, 0.7, 3)


def test_columns_between_two_positions_count_towards_the_span():
    cols = {p: c - 1 for p, c in zip(range(len(SPREAD)), SPREAD)}          # position p -> untrimmed column
    assert ctrain([9, 10, 11], columns=cols) == [] and ctrain([9, 10, 11]) == [9, 10, 11]
    assert ctrain([0, 1, 2], columns={0: 1, 1: 2, 2: 4}) == [0, 1, 2]          # 3 of 4 = 0.75 >= 0.7
    assert ctrain([0, 1, 2], columns={0: 1, 1: 3, 2: 5}) == []                 # 3 of 5 = 0.6 < 0.7


def test_flags_in_untrimmed_columns_are_a_subset_of_the_trimmed_flags_with_the_default_parameters():
    rng = random.Random(9)
    for _ in range(300):
        ori = sorted(rng.sample(range(1, 121), 70))
        cols = {p: c for p, c in enumerate(ori)}
        pos = rng.sample(range(70), rng.randint(0, 30))
        assert set(ctrain(pos, 0.7, 3, cols)) <= set(ctrain(pos, 0.7, 3))


def test_a_position_without_a_column_or_two_sharing_one_raise():
    with pytest.raises(ValueError, match="no column"):
        ctrain([0, 1, 2], columns={0: 1, 1: 2})
    with pytest.raises(ValueError, match="share a column"):
        ctrain([0, 1, 2], columns={0: 1, 1: 1, 2: 2})


def test_train_flags_passes_the_columns_to_every_key():
    cols = {p: c for p, c in enumerate([1, 2, 4, 50, 51, 52])}
    f = train_flags({("b_0", "US"): [0, 1, 2], ("b_1", "US"): [3, 4, 5]}, 0.7, 3, cols)
    assert f == {("b_0", "US"): {0, 1, 2}, ("b_1", "US"): {3, 4, 5}}
    assert train_flags({"k": [0, 1, 2]}, 0.7, 3, {0: 1, 1: 3, 2: 5}) == {"k": set()}


# ── observed filter and null finalize ────────────────────────────────────────

def _observed(tmp_path, rows, *extra):
    disc = tmp_path / "disc.tsv"
    pd.DataFrame(rows, columns=["Gene", "Position", "caap_group"]).to_csv(disc, sep="\t", index=False)
    r = subprocess.run([sys.executable, str(FILTER), "-i", str(disc), *extra], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:] + r.stdout[-400:]
    out = pd.read_csv(tmp_path / "disc.filtered.minlen3.maxcaas70.tsv", sep="\t")
    return {(g, p) for g, p, f in zip(out["Gene"], out["Position"], out["clustering_flag"]) if f == "Discarded"}


ROWS = ([("A", p, "US") for p in (9, 10, 11, 40)]                 # spread by removed columns in untrimmed coordinates
        + [("C", p, "US") for p in (3, 4, 5, 40)]                 # contiguous in both coordinates
        + [("N", p, "US") for p in (9, 10, 11, 40)])              # no MAP file


def _maps(tmp_path):
    d = tmp_path / "maps"
    d.mkdir()
    _write_map(d / "A.Homo_sapiens.map.tsv", SPREAD, 59 + 0)
    _write_map(d / "C.Lemur_catta.map.tsv", list(range(1, 60)), 59)
    return d


def test_observed_filter_without_a_map_directory_is_what_it_was(tmp_path):
    got = _observed(tmp_path, ROWS)
    assert got == {("A", 9), ("A", 10), ("A", 11), ("C", 3), ("C", 4), ("C", 5), ("N", 9), ("N", 10), ("N", 11)}


def test_observed_filter_with_a_map_directory_unflags_spread_positions_and_keeps_genes_without_a_map(tmp_path):
    got = _observed(tmp_path, ROWS, "--map-dir", str(_maps(tmp_path)))
    assert got == {("C", 3), ("C", 4), ("C", 5), ("N", 9), ("N", 10), ("N", 11)}   # A is spread out, N has no MAP


Rec = namedtuple("Rec", "position caap_group asr_path_score side")


def _null_flags(gene, positions, columns):
    recs = [Rec(p, "US", 0.5, "top") for p in positions]
    _, rows = _perms_worker_finalize(gene, [("c1", recs)], 1, True, 3, 0.7, columns)
    return {r["Position"] for r in rows if r["clust"] == 1}


def test_null_clust_flags_follow_the_same_columns(tmp_path):
    idx = index_files(_maps(tmp_path), ".map.tsv")
    assert _null_flags("A", [9, 10, 11, 40], None) == {9, 10, 11}
    assert _null_flags("A", [9, 10, 11, 40], gene_columns(idx, "A", ".map.tsv")) == set()
    assert _null_flags("C", [3, 4, 5, 40], gene_columns(idx, "C", ".map.tsv")) == {3, 4, 5}


def test_observed_and_null_flag_the_same_positions_with_the_same_map_directory(tmp_path):
    maps = _maps(tmp_path)
    idx = index_files(maps, ".map.tsv")
    observed = _observed(tmp_path, ROWS, "--map-dir", str(maps))
    null = set()
    for gene, pos in (("A", [9, 10, 11, 40]), ("C", [3, 4, 5, 40]), ("N", [9, 10, 11, 40])):
        missing = []
        cols = _gene_train_columns(idx, gene, ".map.tsv", missing)
        assert missing == (["N"] if gene == "N" else [])
        null |= {(gene, p) for p in _null_flags(gene, pos, cols)}
    assert null == observed


def test_the_helper_returns_none_without_an_index():
    missing = []
    assert _gene_train_columns(None, "A", ".map.tsv", missing) is None and missing == []
