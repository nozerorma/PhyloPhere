#!/usr/bin/env python3
"""Cluster trains in trimmed versus untrimmed alignment coordinates (Phase 1T, measurement only).

discovery.tab positions are columns of the trimmed alignment, so `ctrain` measures span and density
between trimmed columns. Columns the trimmer removed cannot hold a CAAS; measured on the untrimmed
alignment they add to the span and leave the count unchanged, so two CAAS that are adjacent only
because columns between them were removed are further apart than the trimmed coordinates say.

For every (gene, caap_group) unit, over the positions of all hypotheses, and for every combination of
--minlen-values x --maxcaas-values (the exploratory post-processing grid), this script computes:

* the trains in trimmed coordinates (what the pipeline does);
* the trains in untrimmed coordinates (positions mapped through the trimmer's MAP table).

A train component is a set of flagged positions linked by overlapping qualifying intervals; its size
is the number of discovered positions it holds. Components are reported by size (2, 3, 4, 5-9, 10+) as
kept whole, partly lost or fully lost when the untrimmed coordinates are used. Positions flagged in only
one of the two coordinates are counted separately. With maxcaas above (minlen - 1) / minlen the span can
only grow and every flagged window already holds at least minlen positions, so the untrimmed flags are a
subset of the trimmed ones; below that they need not be.

Usage:
  python trains_coordinates.py --discovery discovery.tab[.gz] --map-dir DIR [--map-suffix .map.tsv]
      [--minlen-values 3 --maxcaas-values 0.7 --out-tsv grid.tsv --components-tsv components.tsv]
"""

import argparse
import collections
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.postproc import ctrain  # noqa: E402
from trains_grain import load  # noqa: E402
from trains_quality import find_file, index_files, load_map  # noqa: E402

BINS = (("2", 2, 2), ("3", 3, 3), ("4", 4, 4), ("5-9", 5, 9), ("10+", 10, 10 ** 9))
TOTALS = ["units", "positions", "flagged_trimmed", "flagged_untrimmed", "only_trimmed", "only_untrimmed", "records",
          "records_removed_trimmed", "records_removed_untrimmed", "units_with_train_trimmed",
          "units_with_train_untrimmed"]


def load_units(df):
    """(gene, caap_group) -> {hypothesis: sorted distinct positions}."""
    units = {}
    for (g, grp, h), pos in df.groupby(["gene", "caap_group", "trait"], observed=True)["position"]:
        units.setdefault((g, grp), {})[h] = sorted({int(p) for p in pos})
    return units


def train_components(positions, maxcaas=0.7, minlen=3):
    """Flagged positions grouped into components: intervals that share a discovered position are merged.

    Same definition as core.postproc.ctrain: a unit with fewer than minlen positions has none; otherwise
    an interval qualifies when span >= minlen and count/span >= maxcaas.
    """
    uniq = sorted({int(p) for p in positions})
    n = len(uniq)
    if n < minlen:
        return []
    spans = []
    for r in range(n):
        for l in range(r + 1):
            span = uniq[r] - uniq[l] + 1
            if span >= minlen and (r - l + 1) / span >= maxcaas:
                spans.append((l, r))
    spans.sort()
    comps, cur = [], None
    for l, r in spans:
        if cur and l <= cur[1]:
            cur[1] = max(cur[1], r)
        else:
            cur = [l, r]
            comps.append(cur)
    return [uniq[l:r + 1] for l, r in comps]


def untrimmed(positions, ori_of_prot):
    """Untrimmed 1-based columns of 0-based trimmed discovery positions."""
    return {int(p): ori_of_prot[int(p) + 1] for p in positions}


def unit_flags(union, ori_of_prot, maxcaas=0.7, minlen=3):
    """(components in trimmed coordinates, flagged positions in trimmed and in untrimmed coordinates)."""
    comps = train_components(union, maxcaas, minlen)
    flagged_t = set(ctrain(sorted(union), maxcaas, minlen))
    col = untrimmed(union, ori_of_prot)
    back = {c: p for p, c in col.items()}
    flagged_u = {back[c] for c in ctrain(sorted(col.values()), maxcaas, minlen)}
    return comps, flagged_t, flagged_u


def size_bin(n):
    return next(name for name, lo, hi in BINS if lo <= n <= hi)


def measure(units, maps, maxcaas=0.7, minlen=3):
    """Component table and totals; units of genes without a usable map are counted in `skipped`."""
    comp_rows, skipped = [], collections.Counter()
    tot = collections.Counter()
    for (gene, grp), byh in units.items():
        if gene not in maps:
            skipped[gene] += 1
            continue
        union = sorted(set().union(*map(set, byh.values())))
        try:
            comps, ft, fu = unit_flags(union, maps[gene], maxcaas, minlen)
        except KeyError:
            skipped[gene] += 1
            continue
        tot["units"] += 1
        tot["positions"] += len(union)
        tot["flagged_trimmed"] += len(ft)
        tot["flagged_untrimmed"] += len(fu)
        tot["only_trimmed"] += len(ft - fu)
        tot["only_untrimmed"] += len(fu - ft)
        tot["records"] += sum(len(p) for p in byh.values())
        tot["records_removed_trimmed"] += sum(len(ft.intersection(p)) for p in byh.values())
        tot["records_removed_untrimmed"] += sum(len(fu.intersection(p)) for p in byh.values())
        tot["units_with_train_trimmed"] += bool(ft)
        tot["units_with_train_untrimmed"] += bool(fu)
        for c in comps:
            kept = sum(p in fu for p in c)
            comp_rows.append({"gene": gene, "caap_group": grp, "size": len(c), "bin": size_bin(len(c)),
                              "kept": kept, "status": "whole" if kept == len(c) else "lost" if kept == 0 else "partly"})
    return pd.DataFrame(comp_rows), tot, skipped


def grid(units, maps, minlens, maxcaases):
    """One row per (minlen, maxcaas) with the totals, and the component outcomes by size bin."""
    rows, comp_tables, skipped_units = [], [], 0
    for minlen in minlens:
        for maxcaas in maxcaases:
            comps, tot, skipped = measure(units, maps, maxcaas, minlen)
            skipped_units = sum(skipped.values())
            rows.append({"minlen": minlen, "maxcaas": maxcaas, **{k: tot[k] for k in TOTALS}})
            if len(comps):
                t = comps.groupby(["bin", "status"]).size().rename("n").reset_index()
                comp_tables.append(t.assign(minlen=minlen, maxcaas=maxcaas))
    comps = pd.concat(comp_tables, ignore_index=True) if comp_tables else pd.DataFrame(
        columns=["bin", "status", "n", "minlen", "maxcaas"])
    return pd.DataFrame(rows), comps, skipped_units


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--map-suffix", default=".map.tsv")
    ap.add_argument("--minlen-values", default="3")
    ap.add_argument("--maxcaas-values", default="0.7")
    ap.add_argument("--out-tsv")
    ap.add_argument("--components-tsv")
    a = ap.parse_args()
    df = load(a.discovery)
    idx = index_files(a.map_dir, a.map_suffix)
    maps, bad = {}, collections.Counter()
    for g in df["gene"].unique():
        try:
            maps[g] = load_map(find_file(idx, g, a.map_suffix))[1]
        except (OSError, ValueError) as e:
            bad[type(e).__name__] += 1
    g, comps, skipped_units = grid(load_units(df), maps, [int(x) for x in a.minlen_values.split(",")],
                                   [float(x) for x in a.maxcaas_values.split(",")])
    pd.set_option("display.width", 250)
    print(f"genes with a map: {len(maps)}  without or unusable: {dict(bad)}  units skipped: {skipped_units}")
    show = g.assign(frac_records_dropped=(1 - g.records_removed_untrimmed / g.records_removed_trimmed).round(3))
    print(show[["minlen", "maxcaas", "units", "positions", "flagged_trimmed", "flagged_untrimmed", "only_trimmed",
                "only_untrimmed", "records_removed_trimmed", "records_removed_untrimmed", "frac_records_dropped",
                "units_with_train_trimmed", "units_with_train_untrimmed"]].to_string(index=False))
    if len(comps):
        w = comps.pivot_table(index=["minlen", "maxcaas", "bin"], columns="status", values="n", fill_value=0).astype(int)
        for s in ("whole", "partly", "lost"):
            if s not in w:
                w[s] = 0
        w["components"] = w[["whole", "partly", "lost"]].sum(axis=1)
        w["frac_lost"] = (w["lost"] / w["components"]).round(3)
        order = {b[0]: i for i, b in enumerate(BINS)}
        w = w.reset_index().sort_values(["minlen", "maxcaas", "bin"], key=lambda c: c.map(order) if c.name == "bin" else c)
        print(w[["minlen", "maxcaas", "bin", "components", "whole", "partly", "lost", "frac_lost"]].to_string(index=False))
    if a.out_tsv:
        g.to_csv(a.out_tsv, sep="\t", index=False)
    if a.components_tsv:
        comps.to_csv(a.components_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
