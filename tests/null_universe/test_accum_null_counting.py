#!/usr/bin/env python3
# test_accum_null_counting.py: the accumulation permulation null counts hits with the rule of the observed side,
# from the filtered null table SCORING also uses, and takes N from the cycle roster.
# Run: python3 tests/null_universe/test_accum_null_counting.py   (scratch files go to ~/.cc/null_universe_test, removed afterwards)

import argparse
import gzip
import importlib.util
import logging
import shutil
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RANDOMIZE = HERE.parents[1] / "subworkflows" / "CT_ACCUMULATION" / "local" / "src" / "randomization" / "randomize.py"
WORK = Path.home() / ".cc" / "null_universe_test"

# randomize.py imports pyarrow for its parquet outputs of the occupancy nulls, which the permulation null does not use
for name in ("pyarrow", "pyarrow.parquet"):
    sys.modules.setdefault(name, types.ModuleType(name))
spec = importlib.util.spec_from_file_location("randomize", RANDOMIZE)
rz = importlib.util.module_from_spec(spec)
sys.modules["randomize"] = rz       # the worker functions of the occupancy nulls are pickled by module name
spec.loader.exec_module(rz)
logging.getLogger().setLevel(logging.ERROR)


def write_null(path, rows):
    """perm_pos_cycle_caas.tsv.gz rows: (Gene, Position, side, cycle, caas_score)."""
    with gzip.open(path, "wt") as fh:
        fh.write("Gene\tPosition\tside\tcycle\tcaas_score\tn_schemes\tscheme_set\tscore_aggregation\n")
        for g, p, s, c, sc in rows:
            fh.write(f"{g}\t{p}\t{s}\t{c}\t{'' if sc is None else sc}\t1\tUS\tus_pos_species\n")


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"ok   {msg}")


def test_rule_equivalence():
    """Observed rows replayed as one null cycle give the same per-gene counts under every direction."""
    # (gene, position, side, scheme): top, bottom, both-sides, side none, GS-only, duplicates
    observed = pd.DataFrame([
        ("G1", 1, "top", "US"), ("G1", 2, "bottom", "US"), ("G1", 3, "top", "US"), ("G1", 3, "bottom", "US"),
        ("G1", 4, "none", "US"), ("G1", 5, "top", "GS2"), ("G2", 1, "top", "US"), ("G2", 1, "top", "US"),
        ("G3", 7, "bottom", "US"),
    ], columns=["gene", "msa_pos", "side", "caap_group"])
    observed = observed.drop_duplicates(["gene", "msa_pos", "caap_group", "side"])   # the dedup of main()
    for direction in ("top", "bottom", "both"):
        us = observed["caap_group"].eq("US") & rz._side_mask(observed["side"], direction)
        obs_counts = observed[us].groupby("gene").size().to_dict()
        # the null table of one cycle: the same detections as US rows with a positive score, plus noise
        # the filter must drop (score 0, empty score, side none, a GS-only detection has no score)
        rows = [(r.gene, r.msa_pos, r.side, "c1", 0.5) for r in observed[observed["caap_group"] == "US"].itertuples()]
        rows += [("G1", 9, "top", "c1", 0.0), ("G1", 10, "top", "c1", None), ("G1", 11, "none", "c1", 0.4),
                 ("G1", 1, "top", "c1", 0.5)]                                         # exact duplicate of a kept row
        path = WORK / "null_equiv.tsv.gz"
        write_null(path, rows)
        counts, present = rz.null_counts_by_gene_cycle(path, direction)
        null_counts = {g: int(v) for (g, c), v in counts.items()}
        check(null_counts == obs_counts, f"direction {direction}: null counts {null_counts} == observed counts {obs_counts}")
    # 'both' counts a position detected on both sides twice, on the observed and on the null side alike
    check(obs_counts["G1"] == 4, "direction both: G1 has 4 rows (top 1, bottom 1, position 3 on both sides)")


def test_roster_and_n():
    path = WORK / "null_n.tsv.gz"
    write_null(path, [("G1", 1, "top", "b_1", 0.5), ("G1", 1, "top", "b_2", 0.5), ("G1", 2, "top", "b_3", None)])
    _, present = rz.null_counts_by_gene_cycle(path, "both")
    check(present == {"b_1", "b_2", "b_3"}, "cycles present are taken before the score filter")
    roster = {"b_1", "b_2", "b_3", "b_4", "b_5"}
    check(rz.resolve_n_cycles(roster, present) == 5, "roster with cycles that detected nothing gives N = 5")
    check(rz.resolve_n_cycles(set(), present) == 3, "no roster: N = cycles present")
    check(rz.resolve_n_cycles({"b_1"}, present) == 3, "roster that misses a present cycle is ignored")
    (WORK / "roster.txt").write_text("b_1\nb_2\nb_3\nb_4\nb_5\n")
    check(rz.read_cycle_roster(WORK / "roster.txt") == roster, "roster file is read")
    check(rz.read_cycle_roster(None) == set(), "no roster path gives the empty set")


def test_p_values_vs_bruteforce():
    """p = (k + 1) / (N + 1) with k the cycles whose null count reaches the observed count, empty cycles included."""
    gene_to_id = {"G1": 0, "G2": 1, "G3": 2}
    rows = []
    # G1: cycles b_1..b_3 have 1, 2, 3 hits; b_4, b_5 none. G2: b_1 has 1 hit. G3: no hit.
    for c, n in (("b_1", 1), ("b_2", 2), ("b_3", 3)):
        rows += [("G1", p, "top", c, 0.3) for p in range(n)]
    rows += [("G2", 0, "bottom", "b_1", 0.2)]
    path = WORK / "null_p.tsv.gz"
    write_null(path, rows)
    counts, present = rz.null_counts_by_gene_cycle(path, "both")
    N = 5
    for obs in (0, 1, 2, 3, 4):
        res = rz.run_permulation_null(counts, N, gene_to_id, 3, {"us": np.array([obs, obs, obs])})[0]
        got = res["us"]["count_above"]
        full = {"G1": [1, 2, 3, 0, 0], "G2": [1, 0, 0, 0, 0], "G3": [0] * 5}
        expect = [sum(1 for v in full[g] if v >= obs) for g in ("G1", "G2", "G3")]
        check(list(got) == expect, f"observed count {obs}: cycles reaching it {list(got)} == {expect} (N = {N})")
    check(res["n_rands"] == N, "n_rands is N, the cycles without a detection included")


def test_main_end_to_end():
    """main() in permulation mode: observed counts from the filtered_discovery rows, null from the null table."""
    out = WORK / "e2e"
    out.mkdir(parents=True)
    # global table: three genes, five positions each (position is the global running index)
    g = []
    pos = 0
    for gene in ("G1", "G2", "G3"):
        for m in range(5):
            g.append(dict(gene=gene, msa_pos=m, position=pos, cons_idx=0.5, masked=False,
                          iscaas=(gene, m) in {("G1", 0), ("G1", 1), ("G2", 2)}))
            pos += 1
    pd.DataFrame(g).to_csv(out / "global.csv", index=False)
    # filtered_discovery: G1 has two positions (one on both sides), G2 one, plus a GS row and a side none row
    fd = [("G1", 0, "top", "US"), ("G1", 1, "top", "US"), ("G1", 1, "bottom", "US"), ("G2", 2, "bottom", "US"),
          ("G2", 3, "top", "GS3"), ("G3", 4, "none", "US")]
    pd.DataFrame(fd, columns=["Gene", "Position", "side", "caap_group"]).to_csv(out / "fd.tsv", sep="\t", index=False)
    # null: 4 replayed cycles, G1 counts per cycle (both directions) 3, 1, 0, 0; G2 1 in every cycle of b_1 and b_2
    rows = [("G1", 0, "top", "b_1", 0.4), ("G1", 1, "top", "b_1", 0.4), ("G1", 1, "bottom", "b_1", 0.4),
            ("G1", 0, "top", "b_2", 0.4), ("G2", 2, "bottom", "b_1", 0.4), ("G2", 2, "bottom", "b_2", 0.4),
            ("G3", 4, "none", "b_3", 0.4)]
    write_null(out / "null.tsv.gz", rows)
    (out / "roster.txt").write_text("b_1\nb_2\nb_3\nb_4\n")
    for direction, exp_actual, exp_above in (("both", {"G1": 3, "G2": 1, "G3": 0}, {"G1": 1, "G2": 2, "G3": 4}),
                                             ("top", {"G1": 2, "G2": 0, "G3": 0}, {"G1": 1, "G2": 4, "G3": 4})):
        args = argparse.Namespace(
            global_csv=str(out / "global.csv"), caas_csv=str(out / "fd.tsv"), background_positions=None, bg_caas=None,
            output_prefix=str(out / f"acc_{direction}"), randomization_type="permulation",
            null_cycle_caas=str(out / "null.tsv.gz"), cycle_roster=str(out / "roster.txt"), n_randomizations=1,
            workers=1, compress=False, export_individual_rand=False, decile_bins=None, global_seed=1,
            precompute_masks=True, change_side=direction, log_level="ERROR")
        rz.main(args)
        res = pd.read_csv(out / f"acc_{direction}_us_aggregated_results.csv").set_index("Gene")
        check(res["ActualCount_us"].to_dict() == exp_actual, f"{direction}: observed counts {res['ActualCount_us'].to_dict()}")
        check(res["NumRands_us"].unique().tolist() == [4], f"{direction}: N = 4 from the roster")
        check(res["RandsAbove_us"].to_dict() == exp_above, f"{direction}: cycles reaching the observed count {res['RandsAbove_us'].to_dict()}")
        # p = (k + 1) / (N + 1)
        check(np.allclose(res["PValueEmpirical_us"], (res["RandsAbove_us"] + 1) / 5), f"{direction}: p = (k + 1) / (N + 1)")


def test_occupancy_modes_unchanged_observed_counts():
    """The naive occupancy null (untouched) still runs and reports the same observed counts as the permulation mode."""
    out = WORK / "e2e"                    # fixture left by test_main_end_to_end
    args = argparse.Namespace(
        global_csv=str(out / "global.csv"), caas_csv=str(out / "fd.tsv"), background_positions=None, bg_caas=None,
        output_prefix=str(out / "acc_naive"), randomization_type="naive", null_cycle_caas=None, cycle_roster=None,
        n_randomizations=20, workers=1, compress=False, export_individual_rand=False, decile_bins=None,
        global_seed=1, precompute_masks=True, change_side="both", log_level="ERROR")
    rz.main(args)
    res = pd.read_csv(out / "acc_naive_us_aggregated_results.csv").set_index("Gene")
    check(res["ActualCount_us"].to_dict() == {"G1": 3, "G2": 1, "G3": 0}, "naive mode: observed counts as in permulation mode")


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    try:
        test_rule_equivalence()
        test_roster_and_n()
        test_p_values_vs_bruteforce()
        test_main_end_to_end()
        test_occupancy_modes_unchanged_observed_counts()
    finally:
        shutil.rmtree(WORK, ignore_errors=True)
    print("all accumulation null-counting tests passed")


if __name__ == "__main__":
    main()
