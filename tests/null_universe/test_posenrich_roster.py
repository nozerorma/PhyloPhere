#!/usr/bin/env python3
# test_posenrich_roster.py: POSENRICH takes N of the null from the cycle roster (cycles with no detection included),
# and falls back to the cycles present with a log line; the prepped null and the direct reader agree.
# Run: python3 tests/null_universe/test_posenrich_roster.py   (scratch files go to ~/.cc/null_universe_test_pe, removed afterwards)

import gzip
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "subworkflows" / "ENRICHMENT" / "local" / "src"
WORK = Path.home() / ".cc" / "null_universe_test_pe"
sys.path.insert(0, str(SRC))

import posenrich_enrich as pe          # noqa: E402
import posenrich_prep_caas_null as pp  # noqa: E402


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"ok   {msg}")


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    try:
        null = WORK / "null.tsv.gz"
        with gzip.open(null, "wt") as fh:
            fh.write("Gene\tPosition\tside\tcycle\tcaas_score\tn_schemes\n")
            for g, p, s, c, sc in (("G1", 1, "top", "b_1", 0.5), ("G1", 2, "bottom", "b_2", 0.4), ("G2", 3, "top", "b_3", "")):
                fh.write(f"{g}\t{p}\t{s}\t{c}\t{sc}\t1\n")
        roster = WORK / "roster.txt"
        roster.write_text("".join(f"b_{i}\n" for i in range(1, 7)))          # b_4..b_6 detected nothing
        short = WORK / "short.txt"
        short.write_text("b_1\nb_2\n")                                       # misses b_3, a cycle present

        lv = pp.resolve_cycle_levels(str(roster), ["b_1", "b_2", "b_3"])
        check(len(lv) == 6, "roster of 6 cycles gives N = 6")
        check(len(pp.resolve_cycle_levels(None, ["b_1", "b_2", "b_3"])) == 3, "no roster: N = cycles present")
        check(len(pp.resolve_cycle_levels(str(short), ["b_1", "b_2", "b_3"])) == 3, "roster missing a present cycle is ignored")
        check(len(pp.resolve_cycle_levels(str(WORK / "NO_FILE"), ["b_1"])) == 1, "a NO_ sentinel counts as no roster")

        # direct reader: columns of the null matrix span the roster, so the empty cycles are draws of zeros
        long_df, levels, positions = pe.load_caas_cycle_null(str(null), str(roster))
        check(len(levels) == 6, "load_caas_cycle_null: N = 6 with the roster")
        check(positions == {"G1:1", "G1:2"}, "null positions are the ones with caas_score > 0")
        sub = pe.null_direction_subset(long_df, "global")
        m = pe.caas_null_matrix({"G1:1": 0, "G1:2": 1}, 2, sub, levels)
        check(m.shape == (2, 6), "null matrix has one column per roster cycle")
        _, levels_nr, _ = pe.load_caas_cycle_null(str(null), None)
        check(len(levels_nr) == 3, "load_caas_cycle_null without a roster: cycles present")

        # prepped artifact
        out = WORK / "prepped.pkl"
        r = subprocess.run([sys.executable, str(SRC / "posenrich_prep_caas_null.py"), "--caas-cycle-null", str(null),
                            "--cycle-roster", str(roster), "--output", str(out)], capture_output=True, text=True, check=True)
        with open(out, "rb") as fh:
            prepped = pickle.load(fh)
        check(len(prepped["cycle_levels"]) == 6 and "replay roster" in r.stdout, "prepped null carries the 6 roster cycles")
        r = subprocess.run([sys.executable, str(SRC / "posenrich_prep_caas_null.py"), "--caas-cycle-null", str(null),
                            "--output", str(out)], capture_output=True, text=True, check=True)
        check("no cycle roster" in r.stdout, "without a roster the prep step says N = cycles present")
    finally:
        shutil.rmtree(WORK, ignore_errors=True)
    print("all POSENRICH roster tests passed")


if __name__ == "__main__":
    main()
