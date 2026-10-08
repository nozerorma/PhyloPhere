#!/usr/bin/env python3
"""Position tables of a Tier 1 PEPC run, as Markdown.

For each trait of a results directory it prints
  1. the ten truth sites (score, rank, p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact) and the positions outside the truth
     set that pass each adjustment at --alpha;
  2. the number of positions at p <= 0.001 / 0.01 / 0.05 and at the floor 1 / (N + 1), for p.emp and p.emp_fact;
  3. the null of each truth site, the null scale and the multiple-testing families;
  4. with --previous, a comparison with the position tables of another run of the same traits (scores are not comparable
     one to one when the runs aggregate the scheme scores differently, so ranks and p-values are given).

Truth positions are in maize PEPC1 numbering: the `Position` column (0-based) plus 1. A position seen on both sides is
reported once, with its best score and its smallest p.

    python3 validation/tier1/scripts/pepc_pvalue_tables.py --results validation/tier1/output/pepc/results \\
        [--previous <dir with <trait>/position_scores.tsv>] [--alpha 0.05]
"""
import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TRUTH = ROOT / "validation/truthsets/tier1/pepc_c4.sites.tsv"
TRAITS = ["c4_complete", "c4_phenotypic_complete"]
P_COLS = ["p.emp", "p.adj_bh", "p.emp_fact", "p.adj_bh_fact"]


def per_position(path):
    """One row per position: best score over sides, smallest p of each kind, score rank (ties share the best rank)."""
    d = pd.read_csv(path, sep="\t")
    agg = {"score": ("CAAS_score", "max")}
    agg.update({c: (c, "min") for c in P_COLS if c in d})
    g = d.groupby("Position").agg(**agg).reset_index()
    g["site"] = g["Position"] + 1
    g["rank"] = g["score"].rank(ascending=False, method="min").astype(int)
    return g


def n_cycles(results):
    with gzip.open(results / "caas_permulation/perm_pos_cycle_caas.tsv.gz", "rt") as f:
        cyc = pd.read_csv(f, sep="\t", usecols=["cycle"])["cycle"]
    return cyc.str.replace(r"~H.*$", "", regex=True).nunique()


def md(df):
    # plain Markdown table, so the script needs nothing beyond pandas
    if not len(df):
        return "(none)"
    cell = lambda x: "" if pd.isna(x) else str(x)
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(cell(x) for x in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def truth_table(g, truth, alpha):
    t = truth.merge(g, left_on="position", right_on="site", how="left")
    cols = ["position", "tier", "score", "rank"] + [c for c in P_COLS if c in g]
    out = t[cols].rename(columns={"position": "site"}).round(4)
    out["score"] = out["score"].map(lambda x: "absent" if pd.isna(x) else f"{x:.3f}")
    print(md(out.fillna("")), "\n")
    for adj in ("p.adj_bh", "p.adj_bh_fact"):
        if adj not in g:
            continue
        hit = g[g[adj] <= alpha]
        in_truth = hit[hit.site.isin(truth.position)]
        other = hit[~hit.site.isin(truth.position)]
        print(f"- `{adj}` <= {alpha}: {len(in_truth)}/{len(truth)} truth sites ({', '.join(map(str, sorted(in_truth.site)))}); "
              f"{len(other)} other positions: "
              + "; ".join(f"{r['site']} (score {r['score']:.3f}, {adj} {r[adj]:.4f})" for _, r in other.iterrows()))
    print()


def tail_table(g, n):
    floor = 1 / (n + 1)
    rows = []
    for col in ("p.emp", "p.emp_fact"):
        if col not in g:
            continue
        p = g[col]
        rows.append({"p": col, "<=0.001": int((p <= 0.001).sum()), "<=0.01": int((p <= 0.01).sum()),
                     "<=0.05": int((p <= 0.05).sum()), f"at floor 1/{n + 1}": int((p <= floor + 1e-12).sum()),
                     "positions": len(p)})
    print(md(pd.DataFrame(rows)), "\n")


def compare(g, prev, truth, alpha):
    m = g.merge(prev, on="site", suffixes=("", "_prev"), how="outer", indicator=True)
    both = m[m["_merge"] == "both"]
    print(f"- positions: {len(g)} here, {len(prev)} in the other run; shared {len(both)}; "
          f"only here {sorted(m[m['_merge'] == 'left_only'].site)}; only there {sorted(m[m['_merge'] == 'right_only'].site)}")
    print(f"- Spearman, score {both['score'].corr(both['score_prev'], method='spearman'):.2f}; "
          f"p.emp {both['p.emp'].corr(both['p.emp_prev'], method='spearman'):.2f}")
    print(f"- positions with `p.adj_bh` <= {alpha}: {int((g['p.adj_bh'] <= alpha).sum())} here, "
          f"{int((prev['p.adj_bh'] <= alpha).sum())} there\n")
    t = truth.merge(m, left_on="position", right_on="site", how="left")
    out = t[["position", "rank", "rank_prev", "p.emp", "p.emp_prev", "p.adj_bh", "p.adj_bh_prev"]].rename(
        columns={"position": "site", "rank": "rank here", "rank_prev": "rank there", "p.emp": "p.emp here",
                 "p.emp_prev": "p.emp there", "p.adj_bh": "BH here", "p.adj_bh_prev": "BH there"}).round(4)
    print(md(out.fillna("")), "\n")


TIE_TOL = 1e-12                         # same constant as scoring_compute.R


def null_scores(results):
    """Per (position, cycle) best score over sides, from the null file; positions are 0-based, as in the tables."""
    with gzip.open(results / "caas_permulation/perm_pos_cycle_caas.tsv.gz", "rt") as f:
        d = pd.read_csv(f, sep="\t", dtype={"caas_score": "string"})
    d["caas_score"] = pd.to_numeric(d["caas_score"], errors="coerce")
    d = d.dropna(subset=["caas_score"])
    d["cycle"] = d["cycle"].str.replace(r"~H.*$", "", regex=True)
    return d.groupby(["Position", "cycle"], as_index=False)["caas_score"].max()


def storey_discoveries(p, lam=0.5, q=0.1):
    """Storey q-values with pi0 = #{p > lam} / (m (1 - lam)), capped at 1; returns (pi0, number of q <= q)."""
    p = pd.Series(p).sort_values().to_numpy()
    m = len(p)
    pi0 = min(1.0, (p > lam).sum() / (m * (1 - lam)))
    qv = pd.Series(pi0 * m * p / (pd.Series(range(1, m + 1)).to_numpy())).iloc[::-1].cummin().iloc[::-1]
    return pi0, int((qv <= q).sum())


def null_report(res, g, truth, n):
    """Truth sites against their own null, null scale, detection of the null and the multiple-testing families."""
    nul = null_scores(res)
    obs = g.set_index("Position")
    by_pos = {p: x["caas_score"].to_numpy() for p, x in nul.groupby("Position")}
    ann = pd.read_csv(res / "scoring/position_scores.tsv", sep="\t").groupby("Position")["n_hypotheses"].max()
    rows = []
    for t in truth.itertuples():
        pos = t.position - 1
        if pos not in obs.index:
            continue
        s = by_pos.get(pos, pd.Series(dtype=float).to_numpy())
        o = obs.loc[pos]
        q = pd.Series(s).quantile([.5, .9, .99]) if len(s) else pd.Series([float("nan")] * 3)
        rows.append({"site": t.position, "observed CAAS": f"{o['score']:.3f}", "n_hyp /100": int(ann.get(pos, 0)), "null detects": len(s),
                     "null q50": f"{q.iloc[0]:.3f}", "null q90": f"{q.iloc[1]:.3f}", "null q99": f"{q.iloc[2]:.3f}",
                     "k_emp": int((s >= o["score"] - TIE_TOL).sum()), "p.emp": round(o["p.emp"], 4), "p.adj_bh": round(o["p.adj_bh"], 4)})
    print("Null of each truth position (quantiles over the cycles that detect it)\n")
    print(md(pd.DataFrame(rows)), "\n")
    # null scale
    pooled = nul["caas_score"]
    per_cycle = nul.groupby("cycle")["Position"].nunique().reindex(sorted(nul["cycle"].unique()), fill_value=0)
    n_obs = len(g)
    print(f"- observed score quartiles {g['score'].quantile([.25, .5, .75]).round(3).tolist()}; pooled null {pooled.quantile([.25, .5, .75]).round(3).tolist()}")
    print(f"- positions detected: observed {n_obs}; null per cycle median {per_cycle.median():.0f}, q95 {per_cycle.quantile(.95):.0f}; "
          f"P(null >= observed) = {(per_cycle >= n_obs).mean():.3f}")
    # d = fraction of null cycles that detect a candidate position
    d = pd.Series({p: len(by_pos.get(p, [])) / n for p in g["Position"]})
    print(f"- null detection rate d of the candidate positions: median {d.median():.2f}, range {d.min():.2f} to {d.max():.2f}; "
          f"largest observed p.emp {g['p.emp'].max():.2f}")
    universe = set(nul["Position"]) | set(g["Position"])
    bh_cand = pd.Series(g["p.emp"].to_numpy())
    bh_cand = (bh_cand.sort_values().reset_index(drop=True) * len(bh_cand) / pd.Series(range(1, len(bh_cand) + 1))).iloc[::-1].cummin().iloc[::-1]
    pi0, n_st = storey_discoveries(g["p.emp"])
    print(f"- families: candidate set {len(g)}, null universe {len(universe)} (null-only {len(universe - set(g['Position']))} at p = 1)")
    print("\n| procedure | family | positions < 0.1 | minimum |\n|---|---|---|---|")
    print(f"| raw p.emp < 0.05 | none | {int((g['p.emp'] < .05).sum())} | {g['p.emp'].min():.4f} |")
    print(f"| BH of p.emp | candidate set | {int((bh_cand < .1).sum())} | {bh_cand.min():.3f} |")
    print(f"| Storey, lambda = 0.5 (pi0 {pi0:.2f}) | candidate set | {n_st} | |")
    print(f"| BH of p.emp (p.adj_bh) | null universe | {int((g['p.adj_bh'] < .1).sum())} | {g['p.adj_bh'].min():.3f} |")
    print(f"| BH of p.emp_fact (p.adj_bh_fact) | null universe | {int((g['p.adj_bh_fact'] < .1).sum())} | {g['p.adj_bh_fact'].min():.3f} |\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--previous", type=Path, help="directory holding <trait>/position_scores.tsv of another run")
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()
    truth = pd.read_csv(TRUTH, sep="\t", comment="#")[["position", "tier"]]
    for trait in TRAITS:
        res = a.results / trait
        if not (res / "scoring/position_scores.tsv").exists():
            sys.exit(f"missing {res}/scoring/position_scores.tsv")
        g = per_position(res / "scoring/position_scores.tsv")
        n = n_cycles(res)
        print(f"## {trait} ({len(g)} positions, N = {n} null cycles)\n")
        truth_table(g, truth, a.alpha)
        tail_table(g, n)
        null_report(res, g, truth, n)
        if a.previous:
            prev = per_position(a.previous / trait / "position_scores.tsv")
            print(f"### Against {a.previous}\n")
            compare(g, prev, truth, a.alpha)


if __name__ == "__main__":
    main()
