"""trains_quality.py: alignment-quality measures for union-only trains, with every source optional."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "bin"))
import trains_quality as tq  # noqa: E402

GOLDEN = HERE / "golden/pepc_c4_complete"

# Untrimmed alignment of four sequences and six codon columns; columns 2 and 5 are removed by the trimmer.
RAW = ["ATG" "---" "ATG" "ATG" "---" "ATG",
       "ATG" "---" "ATG" "ATG" "---" "ATG",
       "ATG" "ATG" "ATG" "ATG" "---" "---",
       "ATG" "ATG" "NNN" "ATG" "ATG" "ATG"]
GAP_ORI = [0.0, 0.5, 0.25, 0.0, 0.75, 0.25]
SELECTED = [1, 3, 4, 6]


def _write_gene(d, gene="G", g=None, status=None):
    (d / "raw").mkdir(exist_ok=True)
    (d / "raw" / f"{gene}.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(RAW)))
    status = status or ["removed" if c not in SELECTED else "selected" for c in range(1, 7)]
    prot = iter(range(1, 5))
    with open(d / f"{gene}.map.tsv", "w") as fh:
        fh.write("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n")
        for c, s in enumerate(status, 1):
            fh.write(f"{c}\t{s}\t{c}\t{next(prot) if s == 'selected' else 'NA'}\n")
    g = [GAP_ORI[c - 1] for c in SELECTED] if g is None else g
    with open(d / f"{gene}.entropy.tsv", "w") as fh:
        fh.write("gene\tposition\tt\tr\tg\tC_trident\tvariability\tn_seqs\n")
        for p, gp in enumerate(g, 1):
            fh.write(f"{gene}\t{p}\t0\t0\t{gp:.6f}\t0\t{0.1 * p:.6f}\t4\n")


def _sources(d):
    return dict(entropy_dir=d, map_dir=d, raw_dir=d / "raw", entropy_suffix=".entropy.tsv",
                map_suffix=".map.tsv", raw_suffix=".fa")


def test_window_mean_clips_at_the_edges_and_ignores_nan():
    assert tq.window_mean([1, 2, 3, 4], 0, 1) == 1.5
    assert tq.window_mean([1, 2, 3, 4], 3, 2) == 3.0
    assert tq.window_mean([np.nan, 2.0, np.nan], 0, 1) == 2.0
    assert np.isnan(tq.window_mean([np.nan, np.nan], 0, 1))


def test_codon_gap_fraction_counts_any_non_acgt_character():
    assert tq.codon_gap_fraction(RAW).tolist() == GAP_ORI


def test_annotate_gene_reads_each_measure_from_its_own_source(tmp_path):
    _write_gene(tmp_path)
    ent = tq.load_entropy(tmp_path / "G.entropy.tsv")
    removed, ori = tq.load_map(tmp_path / "G.map.tsv")
    gap = tq.codon_gap_fraction(RAW)
    q = tq.annotate_gene([0, 1, 2, 3], ent, removed, ori, gap, window=1).set_index("position")
    # position 1 is trimmed column 2 = untrimmed column 3: neighbours 2..4, column 2 removed
    assert q.loc[1, "n_removed_flank"] == 1 and q.loc[1, "gap_pre"] == 0.25
    assert q.loc[1, "gap_pre_win"] == pytest.approx((0.5 + 0.25 + 0.0) / 3)
    assert q.loc[1, "g_win"] == pytest.approx((0.0 + 0.25 + 0.0) / 3)
    # position 2 is untrimmed column 4: neighbours 3..5, column 5 removed, and its gap fraction enters the window
    assert q.loc[2, "n_removed_flank"] == 1 and q.loc[2, "gap_pre_win"] == pytest.approx((0.25 + 0 + 0.75) / 3)
    # edges: the window is clipped
    assert q.loc[0, "gap_pre_win"] == pytest.approx((0.0 + 0.5) / 2)
    assert q.loc[3, "gap_pre_win"] == pytest.approx((0.75 + 0.25) / 2)
    assert q.loc[1, "variability"] == pytest.approx(0.2)


def test_each_source_can_be_left_out(tmp_path):
    _write_gene(tmp_path)
    ent = tq.load_entropy(tmp_path / "G.entropy.tsv")
    removed, ori = tq.load_map(tmp_path / "G.map.tsv")
    assert set(tq.annotate_gene([1], ent)) == {"position", "g", "g_win", "variability"}
    assert set(tq.annotate_gene([1], None, removed, ori)) == {"position", "n_removed_flank"}
    both = tq.annotate_gene([1], None, removed, ori, tq.codon_gap_fraction(RAW))
    assert set(both) == {"position", "n_removed_flank", "gap_pre", "gap_pre_win"}


def test_inconsistent_inputs_raise(tmp_path):
    _write_gene(tmp_path)
    ent = tq.load_entropy(tmp_path / "G.entropy.tsv")
    removed, ori = tq.load_map(tmp_path / "G.map.tsv")
    gap = tq.codon_gap_fraction(RAW)
    with pytest.raises(ValueError, match="differs from the raw gap"):
        tq.annotate_gene([1], ent.assign(g=ent["g"] + 0.1), removed, ori, gap)
    with pytest.raises(ValueError, match="columns, the map selects"):
        tq.annotate_gene([1], ent.iloc[:3], removed, ori)
    with pytest.raises(ValueError, match="beyond the trimmed"):
        tq.annotate_gene([4], ent, removed, ori, gap)
    with pytest.raises(ValueError, match="untrimmed columns"):
        tq.annotate_gene([1], ent, removed, ori, np.r_[gap, 0.0])


def test_map_with_unknown_status_or_broken_numbering_is_rejected(tmp_path):
    _write_gene(tmp_path, status=["removed", "kept", "selected", "selected", "removed", "selected"])
    with pytest.raises(ValueError, match="status"):
        tq.load_map(tmp_path / "G.map.tsv")
    p = tmp_path / "bad.map.tsv"
    p.write_text("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n1\tselected\t1\t2\n2\tselected\t2\t1\n")
    with pytest.raises(ValueError, match="prot_ali_col"):
        tq.load_map(p)


def test_classes_separate_union_only_both_and_unflagged_positions():
    rows = [("g", "US", "H1", p) for p in (10, 12, 100)] + [("g", "US", "H2", p) for p in (11, 13)]
    rows += [("g", "GS1", "H1", p) for p in (5, 6, 7)]
    c = tq.classify(pd.DataFrame(rows, columns=["gene", "caap_group", "trait", "position"]))
    got = {(r.caap_group, r.position): r.cls for r in c.itertuples()}
    assert {got[("US", p)] for p in (10, 11, 12, 13)} == {"only_union"}
    assert got[("US", 100)] == "none"
    assert {got[("GS1", p)] for p in (5, 6, 7)} == {"both"}


def test_collect_annotates_and_reports_genes_it_cannot_annotate(tmp_path):
    _write_gene(tmp_path)
    d = pd.DataFrame([("G", "US", "H1", p) for p in (0, 1, 2)] + [("MISSING", "US", "H1", 3)],
                     columns=["gene", "caap_group", "trait", "position"])
    rows, skipped = tq.collect(d, **_sources(tmp_path), window=1)
    assert set(rows["gene"]) == {"G"} and list(skipped) == ["MISSING"]
    assert rows["cls"].eq("both").all() and rows["gap_pre"].notna().all()


def test_summary_pairs_genes_and_handles_ties():
    recs = []
    for i in range(6):
        recs += [(f"g{i}", "US", 1, "only_union", 0.5), (f"g{i}", "US", 2, "none", 0.1),
                 (f"g{i}", "US", 3, "both", 0.5)]
    s = tq.summarize(pd.DataFrame(recs, columns=["gene", "caap_group", "position", "cls", "g"])).set_index(["a", "b"])
    r = s.loc[("only_union", "none")]
    assert r["n_genes"] == 6 and r["median_a_minus_b"] == pytest.approx(0.4)
    assert r["wilcoxon_p"] == pytest.approx(2 / 64)  # exact two-sided p, six positive differences
    assert s.loc[("only_union", "both"), "median_a_minus_b"] == 0 and np.isnan(s.loc[("only_union", "both"), "wilcoxon_p"])


def test_summary_averages_within_gene_before_comparing():
    # gene a has many only_union positions with high g, gene b one with low g: the gene means decide
    recs = [("a", "US", p, "only_union", 0.9) for p in range(10)] + [("a", "US", 99, "none", 0.1)]
    recs += [("b", "US", 1, "only_union", 0.2), ("b", "US", 2, "none", 0.1)]
    s = tq.summarize(pd.DataFrame(recs, columns=["gene", "caap_group", "position", "cls", "g"]))
    assert s.loc[0, "mean_a_minus_b"] == pytest.approx((0.8 + 0.1) / 2)


def _pepc_entropy(dst):
    """Entropy table of PEPC.fasta from the pipeline's own functions, in its file format."""
    import compute_variability as cv
    seqs = cv.read_fasta(GOLDEN / "PEPC.fasta")
    per_col, _ = cv.compute_per_column(seqs, cv.henikoff_weights([s for _, s in seqs]))
    with open(dst / "PEPC.entropy.tsv", "w") as fh:
        fh.write("gene\tposition\tt\tr\tg\tC_trident\tvariability\tn_seqs\n")
        for c, (t, r, g, C, v) in enumerate(per_col):
            fh.write(f"PEPC\t{c + 1}\t{t:.6f}\t{r:.6f}\t{g:.6f}\t{C:.6f}\t{v:.6f}\t{len(seqs)}\n")
    return seqs


@pytest.mark.skipif(not (GOLDEN / "discovery.tab.gz").exists(), reason="PEPC golden fixture missing")
def test_pepc_with_the_entropy_table_alone(tmp_path):
    seqs = _pepc_entropy(tmp_path)
    df = tq.load(GOLDEN / "discovery.tab.gz")
    rows, skipped = tq.collect(df, entropy_dir=tmp_path, entropy_suffix=".entropy.tsv", window=3)
    assert not skipped and set(rows.columns) >= {"g", "g_win", "variability"}
    assert not {"n_removed_flank", "gap_pre"} & set(rows.columns)
    assert len(rows) == len(tq.classify(df)) and set(rows["cls"]) <= set(tq.CLASSES)
    # g at a discovered position is the fraction of '-' or 'X' in that 0-based alignment column
    r = rows.iloc[len(rows) // 2]
    col = int(r["position"])
    direct = sum(s[col] in "-X" for _, s in seqs) / len(seqs)
    assert r["g"] == pytest.approx(direct, abs=1e-6)
    s = tq.summarize(rows)
    assert set(s["metric"]) == {"g", "g_win", "variability"} and (s["n_genes"] <= 1).all()
