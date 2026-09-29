#!/usr/bin/env python3
"""Cluster trains: union of hypotheses versus per hypothesis (Phase 1T, measurement only).

A train is a position inside an interval whose density count/span reaches maxcaas (span >= minlen),
computed by core.postproc.ctrain over the positions filed under one (gene, caap_group). The two
grains differ in which positions are filed together:

* union: the positions of every hypothesis pooled, one train universe per (gene, caap_group);
* hypothesis: the positions of each hypothesis on its own.

A train of a subset of positions is a train of any superset (adding positions only raises the
density of an interval), so hypothesis flags are contained in the union's; the union can only flag more.
This script measures how much more, on real discovery tables, and how it grows with the number of
hypotheses (random subsets of size H). It changes nothing in the pipeline.

Usage:
  python trains_grain.py --discovery discovery.tab[.gz] --name cancer --sizes 1,2,3,6,9,12 --reps 5 \
      [--maxcaas 0.7 --minlen 3 --seed 1998 --out-tsv sweep.tsv --png sweep.png]
"""

import argparse
import random
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.postproc import ctrain  # noqa: E402

METRICS = ["positions", "flagged_union", "flagged_any_hyp", "only_union", "only_hyp", "records",
           "records_removed_union", "records_removed_hyp", "units", "units_with_train_union", "units_with_train_hyp"]


class Sweep:
    """Discovery positions per (gene, caap_group, hypothesis), with the per-hypothesis trains cached."""

    def __init__(self, df, maxcaas=0.7, minlen=3):
        self.maxcaas, self.minlen = maxcaas, minlen
        self.units = {}
        for (g, grp, h), pos in df.groupby(["gene", "caap_group", "trait"], observed=True)["position"]:
            self.units.setdefault((g, grp), {})[h] = sorted(set(int(p) for p in pos))
        self.hyp_flags = {(u, h): frozenset(ctrain(p, maxcaas, minlen))
                          for u, byh in self.units.items() for h, p in byh.items()}

    def measure(self, hyps):
        hyps = set(hyps)
        tot = dict.fromkeys(METRICS, 0)
        for u, byh in self.units.items():
            sub = {h: p for h, p in byh.items() if h in hyps}
            if not sub:
                continue
            union = set().union(*map(set, sub.values()))
            fu = set(ctrain(sorted(union), self.maxcaas, self.minlen))
            fh = set().union(*(self.hyp_flags[(u, h)] for h in sub))
            tot["units"] += 1
            tot["positions"] += len(union)
            tot["flagged_union"] += len(fu)
            tot["flagged_any_hyp"] += len(fh)
            tot["only_union"] += len(fu - fh)
            tot["only_hyp"] += len(fh - fu)
            tot["records"] += sum(len(p) for p in sub.values())
            tot["records_removed_union"] += sum(len(fu.intersection(p)) for p in sub.values())
            tot["records_removed_hyp"] += sum(len(self.hyp_flags[(u, h)]) for h in sub)
            tot["units_with_train_union"] += bool(fu)
            tot["units_with_train_hyp"] += bool(fh)
        return tot


def measure(df, hyps, maxcaas=0.7, minlen=3):
    return Sweep(df, maxcaas, minlen).measure(hyps)


def _hyp(label):
    m = re.search(r"H\d+", label)
    return m.group(0) if m else label


def load(path):
    """gene, caap_group, trait (normalised to H<n> when it carries one), position."""
    d = pd.read_csv(path, sep="\t", usecols=["gene", "caap_group", "trait", "position"],
                    dtype={"gene": "category", "caap_group": "category", "trait": "category"})
    d["trait"] = d["trait"].astype(str).map(_hyp)
    return d.astype({"trait": "category"})


def sweep(df, sizes, reps, seed, maxcaas, minlen, name):
    sw = Sweep(df, maxcaas, minlen)
    hyps = sorted(df["trait"].unique())
    rng = random.Random(seed)
    rows = []
    for h in sizes:
        if h > len(hyps):
            continue
        n_rep = 1 if h == len(hyps) else reps
        for r in range(n_rep):
            m = sw.measure(rng.sample(hyps, h))
            assert m["only_hyp"] == 0, "per-hypothesis flags must be contained in the union's"
            rows.append({"fixture": name, "H": h, "rep": r, **m})
    return pd.DataFrame(rows)


def summarize(t):
    g = t.groupby(["fixture", "H"])[METRICS].mean()
    out = pd.DataFrame(index=g.index)
    out["positions"] = g["positions"]
    out["frac_flagged_union"] = g["flagged_union"] / g["positions"]
    out["frac_flagged_hyp"] = g["flagged_any_hyp"] / g["positions"]
    out["only_union"] = g["only_union"]
    out["frac_records_removed_union"] = g["records_removed_union"] / g["records"]
    out["frac_records_removed_hyp"] = g["records_removed_hyp"] / g["records"]
    out["frac_units_train_union"] = g["units_with_train_union"] / g["units"]
    out["frac_units_train_hyp"] = g["units_with_train_hyp"] / g["units"]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--name", default="fixture")
    ap.add_argument("--sizes", default="1,2,3,6,9,12")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--maxcaas", type=float, default=0.7)
    ap.add_argument("--minlen", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1998)
    ap.add_argument("--out-tsv")
    ap.add_argument("--png")
    a = ap.parse_args()
    t = sweep(load(a.discovery), [int(x) for x in a.sizes.split(",")], a.reps, a.seed, a.maxcaas, a.minlen, a.name)
    if a.out_tsv:
        t.to_csv(a.out_tsv, sep="\t", index=False)
    s = summarize(t)
    pd.set_option("display.width", 220)
    print(s.round(4).to_string())
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 3.5))
        for col, lab in (("frac_flagged_union", "union"), ("frac_flagged_hyp", "per hypothesis")):
            ax.plot(s.reset_index()["H"], s[col].to_numpy(), marker="o", label=lab)
        ax.set_xlabel("hypotheses (random subsets)")
        ax.set_ylabel("fraction of positions flagged")
        ax.set_title(a.name)
        ax.legend()
        fig.tight_layout()
        fig.savefig(a.png, dpi=150)


if __name__ == "__main__":
    main()
