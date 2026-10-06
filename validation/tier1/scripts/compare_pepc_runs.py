#!/usr/bin/env python3
"""Two Tier 1 PEPC runs of the same trait, compared at the position level.

For each trait it reports the headline counts of `scoring/position_scores.tsv` (rows, positions, p.emp < 0.05, BH and SAM
counts), the positions held by only one run, how far the scores and the empirical p-values of the shared positions are apart,
the 10 truth-set sites side by side, and the cycles of the permulation null each run holds. Exit code 1 when `--require-equal`
is given and the runs differ beyond `--tol`.

    python3 validation/tier1/scripts/compare_pepc_runs.py --a validation/tier1/previous_work/output/pepc_pre_unification/results \\
        --b validation/tier1/output/pepc/results --require-equal
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TRUTH = ROOT / "validation/truthsets/tier1/pepc_c4.sites.tsv"
TRAITS = ("c4", "c4_phenotypic")


def positions(results, trait):
    """Per position: best score over the sides, smallest p.emp, smallest adjusted p-values (position_scores.tsv)."""
    d = pd.read_csv(Path(results) / f"{trait}_complete/scoring/position_scores.tsv", sep="\t", keep_default_na=False, na_values=["NA"])
    u = d.groupby("Position").agg(score=("CAAS_score", "max"), p=("p.emp", "min"), bh=("p.adj_bh", "min"), sam=("p.adj_sam", "min"))
    return d, u


def headline(d, u):
    return {"rows": len(d), "positions": len(u), "p.emp < 0.05": int((u.p < .05).sum()), "BH < 0.05": int((u.bh < .05).sum()), "BH < 0.1": int((u.bh < .1).sum()),
            "SAM < 0.05": int((u.sam < .05).sum()), "SAM < 0.1": int((u.sam < .1).sum())}


def null_cycles(results, trait):
    """Cycles of the permulation null that left a row in gene_cycle_scores.tsv."""
    return len(pd.read_csv(Path(results) / f"{trait}_complete/caas_permulation/gene_cycle_scores.tsv", sep="\t"))


def compare(a, b, trait, tol=1e-9):
    (da, ua), (db, ub) = positions(a, trait), positions(b, trait)
    shared = ua.index.intersection(ub.index)
    out = {"trait": trait, "a": headline(da, ua), "b": headline(db, ub), "only_a": sorted(ua.index.difference(ub.index)), "only_b": sorted(ub.index.difference(ua.index)),
           "n_shared": len(shared), "max_dscore": float((ub.loc[shared, "score"] - ua.loc[shared, "score"]).abs().max()),
           "p_emp_equal": int((ua.loc[shared, "p"] == ub.loc[shared, "p"]).sum()), "max_dp": float((ub.loc[shared, "p"] - ua.loc[shared, "p"]).abs().max()),
           "max_dadj": float(max((ub.loc[shared, c] - ua.loc[shared, c]).abs().max() for c in ("bh", "sam"))),
           "cycles_a": null_cycles(a, trait), "cycles_b": null_cycles(b, trait)}
    out["equal"] = (not out["only_a"] and not out["only_b"] and out["max_dscore"] <= tol and out["max_dp"] <= tol and out["max_dadj"] <= tol
                    and out["cycles_a"] == out["cycles_b"])
    truth = pd.read_csv(TRUTH, sep="\t", comment="#")
    rows = []
    for t in truth.itertuples():          # Position is the 0-based alignment column: maize position minus 1
        col = t.position - 1
        cell = lambda u, c: "-" if col not in u.index else f"{u.loc[col, c]:.4g}"
        rows.append({"position": t.position, "change": f"{t.ref_aa}>{t.alt_aa}", "tier": t.tier, "score": (cell(ua, "score"), cell(ub, "score")),
                     "p.emp": (cell(ua, "p"), cell(ub, "p")), "BH": (cell(ua, "bh"), cell(ub, "bh")), "SAM": (cell(ua, "sam"), cell(ub, "sam"))})
    out["truth"] = rows
    return out


def show(r):
    print(f"== {r['trait']}")
    print(pd.DataFrame({"a": r["a"], "b": r["b"]}).to_string())
    print(f"positions only in a: {r['only_a']} | only in b: {r['only_b']} | shared: {r['n_shared']}")
    print(f"shared: max |d score| {r['max_dscore']:.3g}; p.emp equal {r['p_emp_equal']}/{r['n_shared']} (max |d| {r['max_dp']:.3g}); max |d adjusted p| {r['max_dadj']:.3g}")
    print(f"null cycles with a row: a {r['cycles_a']}, b {r['cycles_b']}")
    print(pd.DataFrame(r["truth"]).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="results directory of the first run (holds c4_complete, c4_phenotypic_complete)")
    ap.add_argument("--b", required=True, help="results directory of the second run")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--require-equal", action="store_true")
    args = ap.parse_args()
    results = [compare(args.a, args.b, t, args.tol) for t in TRAITS]
    for r in results:
        show(r)
    print("\nruns equal within tolerance:", {r["trait"]: r["equal"] for r in results})
    sys.exit(0 if all(r["equal"] for r in results) or not args.require_equal else 1)


if __name__ == "__main__":
    main()
