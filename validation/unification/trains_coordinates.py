#!/usr/bin/env python3
"""Cluster trains in trimmed versus untrimmed alignment coordinates (Phase 1T, measurement only).

discovery.tab positions are columns of the trimmed alignment, so `ctrain` measures span and density
between trimmed columns. Columns the trimmer removed cannot hold a CAAS; measured on the untrimmed
alignment they add to the span and leave the count unchanged, so two CAAS that are adjacent only
because columns between them were removed are further apart than the trimmed coordinates say.

For every (gene, caap_group) unit, over the positions of all hypotheses, this script computes:

* the trains in trimmed coordinates (what the pipeline does);
* the trains in untrimmed coordinates (positions mapped through the trimmer's MAP table).

A train component is a set of flagged positions linked by overlapping qualifying intervals; its size
is the number of discovered positions it holds. Components are reported by size (3, 4, 5-9, 10+) as
kept whole, partly lost or fully lost when the untrimmed coordinates are used. Positions flagged only
in untrimmed coordinates are counted separately (none with the default parameters).

Usage:
  python trains_coordinates.py --discovery discovery.tab[.gz] --map-dir DIR [--map-suffix .Homo_sapiens.map.tsv]
      [--maxcaas 0.7 --minlen 3 --out-tsv components.tsv]
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
from trains_grain import Sweep, load  # noqa: E402
from trains_quality import load_map  # noqa: E402

BINS = (("3", 3, 3), ("4", 4, 4), ("5-9", 5, 9), ("10+", 10, 10 ** 9))


def train_components(positions, maxcaas=0.7, minlen=3):
    """Flagged positions grouped into components: intervals that share a discovered position are merged.

    Same interval definition as core.postproc.ctrain: span >= minlen and count/span >= maxcaas.
    """
    uniq = sorted({int(p) for p in positions})
    n = len(uniq)
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


def measure(df, maps, maxcaas=0.7, minlen=3):
    """Component table and totals; genes without a usable map are listed in `skipped`."""
    sw = Sweep(df, maxcaas, minlen)
    comp_rows, skipped = [], collections.Counter()
    tot = collections.Counter()
    for (gene, grp), byh in sw.units.items():
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--map-suffix", default=".Homo_sapiens.map.tsv")
    ap.add_argument("--maxcaas", type=float, default=0.7)
    ap.add_argument("--minlen", type=int, default=3)
    ap.add_argument("--out-tsv")
    a = ap.parse_args()
    df = load(a.discovery)
    maps, bad = {}, collections.Counter()
    for g in df["gene"].unique():
        try:
            maps[g] = load_map(Path(a.map_dir) / f"{g}{a.map_suffix}")[1]
        except (OSError, ValueError) as e:
            bad[type(e).__name__] += 1
    comps, tot, skipped = measure(df, maps, a.maxcaas, a.minlen)
    pd.set_option("display.width", 200)
    print(f"genes with a map: {len(maps)}  without or unusable: {dict(bad)}  units skipped: {sum(skipped.values())}")
    f = lambda k: tot[k]  # noqa: E731
    print(f"units {f('units')}  discovered positions {f('positions')}")
    print(f"flagged positions: trimmed {f('flagged_trimmed')}  untrimmed {f('flagged_untrimmed')}  "
          f"only untrimmed {f('only_untrimmed')}")
    print(f"records removed: trimmed {f('records_removed_trimmed')}  untrimmed {f('records_removed_untrimmed')}  "
          f"of {f('records')}")
    print(f"units with a train: trimmed {f('units_with_train_trimmed')}  untrimmed {f('units_with_train_untrimmed')}")
    if len(comps):
        t = comps.groupby(["bin", "status"]).size().unstack(fill_value=0).reindex([b[0] for b in BINS]).fillna(0).astype(int)
        for s in ("whole", "partly", "lost"):
            if s not in t:
                t[s] = 0
        t["components"] = t[["whole", "partly", "lost"]].sum(axis=1)
        t["frac_lost"] = (t["lost"] / t["components"]).round(3)
        print(t[["components", "whole", "partly", "lost", "frac_lost"]].to_string())
    if a.out_tsv:
        comps.to_csv(a.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
