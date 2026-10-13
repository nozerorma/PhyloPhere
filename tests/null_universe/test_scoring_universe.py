#!/usr/bin/env python3
# test_scoring_universe.py: scoring_compute.R on a synthetic fixture, checked against an independent Python recomputation:
# N from the cycle roster (cycles with no detection included), the positional universe U (null-positive positions inside
# the tested set and the cleaned background genes, plus every observed positive position), nd and null_unseen,
# null_universe_summary.tsv, and the fallbacks (no roster, no tested positions, coordinate mismatch).
# Needs Rscript with dplyr, tidyr and readr. Run: python3 tests/null_universe/test_scoring_universe.py
# Scratch files go to ~/.cc/null_universe_test_scoring, removed afterwards. Nothing here runs Nextflow.

import csv
import gzip
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "subworkflows" / "SCORING" / "local" / "src"
CORE = REPO / "subworkflows" / "CT_DISAMBIGUATION" / "local" / "src" / "core" / "scores.py"
WORK = Path.home() / ".cc" / "null_universe_test_scoring"
TIE = 1e-12
N_CYCLES, SILENT = 30, 3          # the last SILENT cycles of the roster detect nothing


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"ok   {msg}")


def build_fixture(d):
    rnd = random.Random(1)
    genes = [f"G{i}" for i in range(1, 41)]
    obs = {}
    with open(d / "fd.tsv", "w") as fh:
        fh.write("Gene\tPosition\tside\tcaap_group\tasr_path_score\tn_hypotheses\tparticipating_hypotheses\n")
        for g in genes:
            for p in rnd.sample(range(55), 4):
                side, sc = rnd.choice(["top", "bottom"]), round(rnd.uniform(0.2, 1), 3)
                obs[(g, p)] = sc
                fh.write(f"{g}\t{p}\t{side}\tUS\t{sc}\t1\t\n")
        fh.write("G1\t70\ttop\tGS2\t0.5\t1\t\n")                         # a GS-only detection: no score
    cycles = [f"b_{i}" for i in range(1, N_CYCLES + 1)]
    null = []
    for g in genes:
        for p in range(60):
            rate = 0.5 if (g, p) in obs and rnd.random() < 0.7 else 0.1
            if p >= 55:
                rate = 0.3                                                 # positions outside the tested set
            for c in cycles[:-SILENT]:
                if rnd.random() < rate * 0.5:
                    null.append((g, p, rnd.choice(["top", "bottom"]), c, round(rnd.uniform(0.1, 1), 3)))
    with gzip.open(d / "null.tsv.gz", "wt") as fh:
        fh.write("Gene\tPosition\tside\tcycle\tcaas_score\tn_schemes\tscheme_set\tscore_aggregation\n")
        for g, p, s, c, sc in null:
            fh.write(f"{g}\t{p}\t{s}\t{c}\t{sc}\t1\tUS\tus_pos\n")
    (d / "roster.txt").write_text("".join(f"{c}\n" for c in cycles))
    with open(d / "background.output", "w") as fh:
        for g in genes:
            fh.write(g + "\t" + ",".join(str(p) for p in range(55)) + "\n")
    with open(d / "background_shifted.output", "w") as fh:               # a coordinate system the null does not use
        for g in genes:
            fh.write(g + "\t" + ",".join(str(p) for p in range(1000, 1055)) + "\n")
    (d / "universe.txt").write_text("".join(f"{g}\n" for g in genes[:-2]))   # the last two genes are outside it
    return genes, obs, null


def run(d, out, extra):
    out.mkdir(parents=True, exist_ok=True)
    for f in ("scoring_compute.R", "aa_grouping.R", "observed_core_scores.py"):
        shutil.copy(SRC / f, out / f)
    (out / "src" / "core").mkdir(parents=True, exist_ok=True)
    shutil.copy(CORE, out / "src" / "core" / "scores.py")
    (out / "src" / "__init__.py").touch()
    (out / "src" / "core" / "__init__.py").touch()
    subprocess.run([sys.executable, "observed_core_scores.py", "--input", str(d / "fd.tsv"), "--positions-out",
                    "core_positions.tsv", "--genes-out", "core_genes.tsv"], cwd=out, check=True, capture_output=True)
    args = ["Rscript", "scoring_compute.R", "--postproc", str(d / "fd.tsv"), "--core_positions", "core_positions.tsv",
            "--core_genes", "core_genes.tsv", "--caas_pos_cycle_caas", str(d / "null.tsv.gz")] + extra
    r = subprocess.run(args, cwd=out, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"scoring_compute.R failed ({Path(out).name}):\n{r.stderr[-1500:]}")
    print(f"ok   scoring_compute.R exits 0 ({Path(out).name})")
    return r.stdout + r.stderr


def table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    try:
        genes, obs, null = build_fixture(WORK)
        # independent recomputation
        in_univ_genes = set(genes[:-2])
        null_pos = {(g, p) for g, p, s, c, sc in null}
        tested = {(g, p) for g in genes for p in range(55)}
        # an observed position stays in U whether or not it lies in the tested set and the gene universe
        U_null = {k for k in null_pos if (k in tested and k[0] in in_univ_genes) or k in obs}
        nd = {}
        best = {}                                                          # (gene, pos, cycle) -> max side score
        for g, p, s, c, sc in null:
            best[(g, p, c)] = max(best.get((g, p, c), 0), sc)
        for (g, p, c) in best:
            nd[(g, p)] = nd.get((g, p), 0) + 1
        unseen = {k for k in obs if nd.get(k, 0) == 0}
        U = U_null | set(obs)
        n_present = len({c for _, _, _, c, _ in null})

        # ── roster + tested positions + universe ──────────────────────────────
        out = WORK / "full"
        log = run(WORK, out, ["--cycle_roster", str(WORK / "roster.txt"), "--tested_positions", str(WORK / "background.output"),
                              "--universe", str(WORK / "universe.txt")])
        check(f"N = {N_CYCLES} cycles (source: cycle_roster.txt)" in log, "N comes from the roster")
        check(n_present == N_CYCLES - SILENT, f"the fixture has {SILENT} cycles with no detection (present: {n_present})")
        s = {r["metric"]: r for r in table(out / "null_universe_summary.tsv") if r["metric"] != "null_positive_positions_cumulative"}
        check(int(float(s["N_cycles"]["value"])) == N_CYCLES, "summary: N = roster size")
        check(int(float(s["N_cycles_present_in_null"]["value"])) == n_present, "summary: cycles present in the null file")
        check(int(float(s["observed_positive_positions"]["value"])) == len(obs), "summary: observed positive positions")
        check(int(float(s["observed_unseen_positions"]["value"])) == len(unseen), f"summary: {len(unseen)} unseen positions")
        check(int(float(s["null_positive_positions_in_universe"]["value"])) == len(U_null), "summary: null-positive positions in U")
        check(int(float(s["universe_size"]["value"])) == len(U), f"summary: |U| = {len(U)}")
        check(f"BH over {len(U)} positions ({len(obs)} observed, {len(U) - len(obs)} null-only at p = 1)" in log,
              "the BH family is U")
        curve = [r for r in table(out / "null_universe_summary.tsv") if r["metric"] == "null_positive_positions_cumulative"]
        vals = [float(r["value"]) for r in curve]
        check(len(curve) == 10 and vals == sorted(vals) and int(vals[-1]) == len(U_null),
              "saturation curve: 10 checkpoints, cumulative, ends at the null-positive positions of U")
        pos = {(r["Gene"], int(r["Position"])): r for r in table(out / "position_scores.tsv")}
        check(all(int(pos[k]["nd"]) == nd.get(k, 0) for k in obs), "nd of every observed position = cycles that score it")
        check(all((pos[k]["null_unseen"] == "TRUE") == (k in unseen) for k in obs), "null_unseen = (nd == 0)")
        check(all(float(pos[k]["p.emp"]) >= 1 / (N_CYCLES + 1) - 1e-12 for k in obs), "p.emp floor is 1/(N + 1) with N from the roster")
        for k in unseen:                                                   # kept, at the floor, not demoted or dropped
            check(abs(float(pos[k]["p.emp"]) - 1 / (N_CYCLES + 1)) < 1e-9 and pos[k]["p.adj_bh"] != "NA",
                  f"unseen position {k[0]}:{k[1]} keeps p = 1/(N + 1) and a BH value")
            break
        k_emp = {k: sum(1 for (g, p, c), v in best.items() if (g, p) == k and v >= obs[k] - TIE) for k in obs}
        check(all(abs(float(pos[k]["p.emp"]) - (k_emp[k] + 1) / (N_CYCLES + 1)) < 1e-9 for k in obs),
              "p.emp = (k + 1) / (N + 1) against brute force")

        # ── fallbacks ────────────────────────────────────────────────────────
        out = WORK / "noroster"
        log = run(WORK, out, ["--tested_positions", str(WORK / "background.output"), "--universe", str(WORK / "universe.txt")])
        check(f"N = {n_present} cycles (source: cycles present in the null file)" in log and "understates N" in log,
              "no roster and no caas_perms.rds: N = cycles present, said in the log")
        out = WORK / "notested"
        log = run(WORK, out, ["--cycle_roster", str(WORK / "roster.txt")])
        check("no tested positions given" in log, "no tested positions: the family is every null-positive position, said in the log")
        U_all = null_pos | set(obs)
        check(f"BH over {len(U_all)} positions" in log, f"fallback family = every null-positive position plus the observed ones ({len(U_all)})")
        out = WORK / "shifted"
        log = run(WORK, out, ["--cycle_roster", str(WORK / "roster.txt"), "--tested_positions", str(WORK / "background_shifted.output")])
        check("coordinates likely differ" in log and f"BH over {len(U_all)} positions" in log,
              "tested positions on other coordinates are ignored with a warning")
        # the summary of the roster-only run still carries N and the curve
        check((out / "null_universe_summary.tsv").exists(), "the summary is written in every run that has a null")
    finally:
        shutil.rmtree(WORK, ignore_errors=True)
    print("all scoring-universe tests passed")


if __name__ == "__main__":
    main()
