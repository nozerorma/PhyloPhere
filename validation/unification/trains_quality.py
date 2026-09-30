#!/usr/bin/env python3
"""Cluster trains and alignment quality (Phase 1T, measurement only).

Positions that only the union of hypotheses flags as a train ("only_union") are compared with
positions flagged under both grains ("both") and with discovered positions that are not flagged
("none"), on independent alignment-quality measures. The question is whether union-only trains sit
in worse alignment regions than the rest of the discovered positions. Nothing in the pipeline changes.

Coordinates: discovery.tab `position` is the 0-based column of the trimmed (BMGE) alignment, so the
1-based trimmed column is position + 1. That column is `entropy.position` and `MAP.prot_ali_col`.

Inputs (only --discovery is required; each measure is computed when its source is given). Files are
named <gene>[.<version>].<species><tail>; a gene is found by the part of the name before the first '.',
whatever its reference species, and a gene with several files in one directory is skipped:

* --entropy-dir   files ending in --entropy-suffix: per-column table with `position`, `g`, `variability`
                  (bin/compute_variability.py). g is the fraction of '-' or 'X' in the column.
                  Gives g, variability and g_win (mean g over +-window trimmed columns).
* --map-dir       files ending in --map-suffix: one row per column of the untrimmed alignment with `status`
                  (selected | removed) and `prot_ali_col`. Gives n_removed_flank: removed columns
                  within +-window untrimmed columns of the position.
* --raw-dir       files ending in --raw-suffix: the untrimmed codon alignment (FASTA). Needs --map-dir.
                  Gives gap_pre and gap_pre_win: fraction of sequences whose codon has a character
                  outside ACGT, at the column and averaged over +-window untrimmed columns
                  (removed columns included).

Higher g, g_win, n_removed_flank, gap_pre and gap_pre_win mean a worse neighbourhood. `variability`
mixes biology and quality (CAAS lie in variable columns by construction); read it as context only.

Comparison: metrics are averaged per gene within each class, then each pair of classes is compared
across the genes that have both (paired Wilcoxon signed-rank). Positions are distinct (gene, position,
class) triples; a position filed under several caap_groups can fall in more than one class.
Per-gene aggregation keeps between-gene differences out of the contrast. Power is limited by the
number of genes with union-only positions.

Usage:
  python trains_quality.py --discovery discovery.tab[.gz] --entropy-dir DIR [--map-dir DIR [--raw-dir DIR]]
      [--window 3 --maxcaas 0.7 --minlen 3 --out-tsv positions.tsv --summary-tsv summary.tsv]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.columns import find_file, index_files, read_map  # noqa: E402
from src.core.postproc import ctrain  # noqa: E402
from trains_grain import Sweep, load  # noqa: E402

METRICS = ["g", "g_win", "variability", "n_removed_flank", "gap_pre", "gap_pre_win"]
CLASSES = ("only_union", "both", "none")
PAIRS = (("only_union", "none"), ("only_union", "both"), ("both", "none"))
G_TOL = 1e-5  # the entropy table is written with 6 decimals


# ── readers ──────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            seqs.append([])
        elif line:
            seqs[-1].append(line)
    return ["".join(s) for s in seqs]


def load_map(path):
    """core.columns.read_map with the removed flags as a boolean array."""
    removed, ori_of_prot = read_map(path)
    return np.asarray(removed, dtype=bool), ori_of_prot


def load_entropy(path):
    """g and variability indexed by 1-based trimmed column; columns must be 1..n."""
    t = pd.read_csv(path, sep="\t", usecols=["position", "g", "variability"])
    if t["position"].tolist() != list(range(1, len(t) + 1)):
        raise ValueError("entropy positions are not 1..n")
    return t.set_index("position")


def codon_gap_fraction(seqs):
    """Per untrimmed column: fraction of sequences whose codon has a character outside ACGT."""
    lens = {len(s) for s in seqs}
    if len(lens) != 1 or lens.pop() % 3:
        raise ValueError("raw alignment rows differ in length or are not codon multiples")
    arr = np.frombuffer("".join(seqs).upper().encode(), dtype=np.uint8).reshape(len(seqs), -1)
    bad = ~np.isin(arr, np.frombuffer(b"ACGT", dtype=np.uint8))
    return bad.reshape(len(seqs), -1, 3).any(axis=2).mean(axis=0)


def window_mean(values, idx, k):
    """Mean of values[idx-k .. idx+k] (0-based, clipped to the array), ignoring NaN."""
    v = np.asarray(values, dtype=float)[max(idx - k, 0): idx + k + 1]
    return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan


# ── per-gene annotation ──────────────────────────────────────────────────────

def annotate_gene(positions, entropy=None, removed=None, ori_of_prot=None, gap_ori=None, window=3):
    """Quality measures for the 0-based discovery positions of one gene; raises on inconsistent inputs."""
    n_cols = None
    if entropy is not None:
        n_cols = len(entropy)
    if ori_of_prot is not None:
        n_sel = len(ori_of_prot)
        if n_cols is not None and n_cols != n_sel:
            raise ValueError(f"entropy has {n_cols} columns, the map selects {n_sel}")
        n_cols = n_sel
    if gap_ori is not None:
        if removed is None or len(gap_ori) != len(removed):
            raise ValueError("raw alignment and map differ in the number of untrimmed columns")
        if entropy is not None:
            sel = np.flatnonzero(~removed)
            gap_sel = gap_ori[sel]
            d = float(np.abs(gap_sel - entropy["g"].to_numpy()).max())
            if d > G_TOL:
                raise ValueError(f"entropy g differs from the raw gap fraction on selected columns (max {d:.3g})")
    if n_cols is not None and len(positions) and max(positions) + 1 > n_cols:
        raise ValueError("a discovery position lies beyond the trimmed alignment")

    g = entropy["g"].to_numpy() if entropy is not None else None
    rows = []
    for p0 in positions:
        p = int(p0) + 1
        r = {"position": int(p0)}
        if entropy is not None:
            r["g"] = float(g[p - 1])
            r["g_win"] = window_mean(g, p - 1, window)
            r["variability"] = float(entropy["variability"].iloc[p - 1])
        if ori_of_prot is not None:
            o = ori_of_prot[p] - 1
            r["n_removed_flank"] = int(removed[max(o - window, 0): o + window + 1].sum())
            if gap_ori is not None:
                r["gap_pre"] = float(gap_ori[o])
                r["gap_pre_win"] = window_mean(gap_ori, o, window)
        rows.append(r)
    return pd.DataFrame(rows)


# ── classes and contrasts ────────────────────────────────────────────────────

def classify(df, maxcaas=0.7, minlen=3):
    """(gene, caap_group, position, cls): only_union | both | none for every discovered position."""
    sw = Sweep(df, maxcaas, minlen)
    out = []
    for (gene, grp), byh in sw.units.items():
        union = sorted(set().union(*map(set, byh.values())))
        fu = set(ctrain(union, maxcaas, minlen))
        fh = set().union(*(sw.hyp_flags[((gene, grp), h)] for h in byh))
        assert fh <= fu, "per-hypothesis flags must be contained in the union's"
        for p in union:
            out.append((gene, grp, p, "both" if p in fh else "only_union" if p in fu else "none"))
    return pd.DataFrame(out, columns=["gene", "caap_group", "position", "cls"])


def collect(df, entropy_dir=None, map_dir=None, raw_dir=None, entropy_suffix=".prot.entropy.tsv",
            map_suffix=".map.tsv", raw_suffix=".fa", window=3, maxcaas=0.7, minlen=3):
    """Classified positions with their quality measures; genes with missing or inconsistent inputs are skipped."""
    cls = classify(df, maxcaas, minlen)
    idx = {k: index_files(d, t) for k, d, t in (("entropy", entropy_dir, entropy_suffix), ("map", map_dir, map_suffix),
                                               ("raw", raw_dir, raw_suffix)) if d}
    pieces, skipped = [], {}
    for gene, sub in cls.groupby("gene", observed=True):
        try:
            ent = load_entropy(find_file(idx["entropy"], gene, entropy_suffix)) if entropy_dir else None
            removed = ori = gap = None
            if map_dir:
                removed, ori = load_map(find_file(idx["map"], gene, map_suffix))
            if raw_dir:
                gap = codon_gap_fraction(read_fasta(find_file(idx["raw"], gene, raw_suffix)))
            pos = sorted(sub["position"].unique())
            q = annotate_gene(pos, ent, removed, ori, gap, window)
        except (OSError, ValueError, KeyError) as e:
            skipped[gene] = f"{type(e).__name__}: {e}"
            continue
        pieces.append(sub.merge(q, on="position", how="left"))
    rows = pd.concat(pieces, ignore_index=True) if pieces else cls.iloc[0:0]
    return rows, skipped


def summarize(rows):
    """Per metric and class pair: genes with both classes, median and mean paired difference, Wilcoxon p."""
    from scipy.stats import wilcoxon
    d = rows.drop_duplicates(["gene", "position", "cls"])
    out = []
    for m in [c for c in METRICS if c in d.columns]:
        per_gene = d.groupby(["gene", "cls"], observed=True)[m].mean().unstack("cls")
        for a, b in PAIRS:
            if a not in per_gene or b not in per_gene:
                continue
            diff = (per_gene[a] - per_gene[b]).dropna()
            p = np.nan
            if len(diff) and (diff != 0).any():
                p = float(wilcoxon(diff, zero_method="wilcox").pvalue)
            out.append({"metric": m, "a": a, "b": b, "n_genes": len(diff),
                        "median_a_minus_b": float(diff.median()) if len(diff) else np.nan,
                        "mean_a_minus_b": float(diff.mean()) if len(diff) else np.nan, "wilcoxon_p": p})
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--entropy-dir")
    ap.add_argument("--map-dir")
    ap.add_argument("--raw-dir")
    ap.add_argument("--entropy-suffix", default=".prot.entropy.tsv")
    ap.add_argument("--map-suffix", default=".map.tsv")
    ap.add_argument("--raw-suffix", default=".fa")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--maxcaas", type=float, default=0.7)
    ap.add_argument("--minlen", type=int, default=3)
    ap.add_argument("--out-tsv")
    ap.add_argument("--summary-tsv")
    a = ap.parse_args()
    if not (a.entropy_dir or a.map_dir):
        ap.error("give at least one of --entropy-dir and --map-dir")
    if a.raw_dir and not a.map_dir:
        ap.error("--raw-dir needs --map-dir")
    rows, skipped = collect(load(a.discovery), a.entropy_dir, a.map_dir, a.raw_dir, a.entropy_suffix,
                            a.map_suffix, a.raw_suffix, a.window, a.maxcaas, a.minlen)
    pd.set_option("display.width", 220)
    d = rows.drop_duplicates(["gene", "position", "cls"])
    print(f"genes annotated: {d['gene'].nunique()}  skipped: {len(skipped)}")
    kinds = {}
    for g, why in skipped.items():
        kinds.setdefault(why.split(":")[0], []).append((g, why))
    for kind, items in kinds.items():
        print(f"  {kind}: {len(items)} genes; first: {items[0][1][:160]}")
    print("positions per class:", d["cls"].value_counts().reindex(CLASSES, fill_value=0).to_dict())
    s = summarize(rows)
    print(s.round(4).to_string(index=False))
    if a.out_tsv:
        rows.to_csv(a.out_tsv, sep="\t", index=False)
    if a.summary_tsv:
        s.to_csv(a.summary_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
