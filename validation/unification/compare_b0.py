#!/usr/bin/env python3
"""Observed results vs the b_0 slice of the permulation-null code path.

The null is calibrated only if its labelings go through the same computation as the
real one. Every run also replays the real labeling (b_0) through the null path
(written to ``<run>/caas_permulation/b0/``); this script compares
that slice with the observed results at five checkpoints:

  A  discovery rows            (gene, hypothesis, caap_group, position, caas, amino_encoded)
  B  post-filter survivors     (gene, position, caap_group, side)
  C  asr_path_score            per (gene, position, caap_group, side)
  D  caas_row / CAAS_score     per (gene, position, side)
  E  gene-level score          gene_caas_score{,_top_all,_bottom_all}

Criterion: exact set equality of keys and |delta| <= tol on values. A checkpoint passes
only if both hold. Exit status is 1 if any selected checkpoint fails.

Comparison notes (what "equal" means at each checkpoint):
  * plain mode (no FOP fan-out in b_0): only the hypotheses present in b_0 are compared at A.
  * B compares the union over hypotheses of the observed survivors with the b_0 rows that
    are neither cluster-flagged nor in a removed (cycle, group, gene) unit.
  * E treats "gene absent / 0 in b_0" and "NA observed" as the same empty-gene state and
    reports how many genes fall in it (``n_empty_convention``), since the two sides encode
    it differently.
"""
import argparse
import glob
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CHECKPOINTS = "ABCDE"


def _hyp(s):
    m = re.search(r"H\d+", str(s))
    return m.group(0) if m else "H1"


def _keyset(df, cols):
    return set(map(tuple, df[cols].astype(str).itertuples(index=False, name=None)))


def compare_sets(a, b, label_a="observed", label_b="b_0", n_examples=5):
    return {
        f"n_{label_a}": len(a), f"n_{label_b}": len(b),
        f"only_{label_a}": len(a - b), f"only_{label_b}": len(b - a),
        "examples_only_" + label_a: sorted(a - b)[:n_examples],
        "examples_only_" + label_b: sorted(b - a)[:n_examples],
    }


def compare_values(df_a, df_b, key, cols_a, cols_b, tol, n_examples=5, bitwise=False):
    """Inner-join two tables on `key`; report set differences and value deltas.

    With `bitwise`, also report how many compared values differ in any bit (informational: only
    meaningful when both tables were read and written without rounding)."""
    res = compare_sets(_keyset(df_a, key), _keyset(df_b, key))
    m = df_a.merge(df_b, on=key, suffixes=("_obs", "_b0"))
    worst, n_bad, n_bits, ex = 0.0, 0, 0, []
    for ca, cb in zip(cols_a, cols_b):
        x = pd.to_numeric(m[ca if ca != cb else ca + "_obs"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(m[cb if ca != cb else cb + "_b0"], errors="coerce").to_numpy(float)
        both_nan = np.isnan(x) & np.isnan(y)
        d = np.where(both_nan, 0.0, np.abs(x - y))  # one-sided NaN -> nan -> counted as bad
        bad = ~(d <= tol)
        n_bad += int(bad.sum())
        n_bits += int((~both_nan & (x != y)).sum())
        if len(d):
            worst = max(worst, float(np.nanmax(np.where(np.isnan(d), np.inf, d))))
        for i in np.flatnonzero(bad)[:n_examples]:
            ex.append({**{k: str(m[k].iloc[i]) for k in key}, "col": ca, "observed": x[i], "b_0": y[i]})
    res.update(n_compared=len(m), n_value_mismatch=n_bad, max_abs_delta=worst, examples_value=ex[:n_examples])
    if bitwise:
        res["n_bitwise_different"] = n_bits
    res["pass"] = res["only_observed"] == 0 and res["only_b_0"] == 0 and n_bad == 0
    return res


def read_detail(b0_dir):
    """b_0 per-(gene, position, group, side) detail shards written by the null path."""
    shards = sorted(glob.glob(str(Path(b0_dir) / "perm_pos_detail" / "*.tsv.gz")))
    if not shards:
        sys.exit(f"ERROR: no b_0 detail shards under {b0_dir}/perm_pos_detail (the run has no b_0 slice)")
    return pd.concat((pd.read_csv(s, sep="\t", float_precision="round_trip") for s in shards), ignore_index=True)


def checkpoint_A(run, b0_dir, tol, perm_disc=None):
    obs = pd.read_csv(Path(run) / "caastools" / "discovery.tab", sep="\t")
    obs = obs.assign(hyp=obs["trait"].map(_hyp), position=obs["position"].astype(int))
    rows = []
    for f in sorted(glob.glob(str(Path(perm_disc or Path(run) / "caas_permulation" / "perm_disc") / "*.perm_replay.discovery.output"))):
        d = pd.read_csv(f, sep="\t")
        rows.append(d[d["cycle"].astype(str).str.match(r"^b_0(~|$)")])
    b0 = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["cycle", "gene", "caap_group", "position", "caas", "amino_encoded"])
    b0 = b0.assign(hyp=b0["cycle"].map(_hyp), position=b0["position"].astype(int))
    obs = obs[obs["hyp"].isin(set(b0["hyp"]))]  # plain mode carries only H1
    key = ["gene", "hyp", "caap_group", "position", "caas", "amino_encoded"]
    res = compare_sets(_keyset(obs, key), _keyset(b0, key))
    res["pass"] = res["only_observed"] == 0 and res["only_b_0"] == 0
    res["hypotheses_compared"] = len(set(b0["hyp"]))
    return res


def survivors(detail, removed):
    d = detail[detail["clust"].astype(int) == 0]
    d = d[~d.apply(lambda r: ("b_0", r["caap_group"], r["Gene"]) in removed, axis=1)]
    return d


def checkpoint_B(run, b0_dir, tol):
    fd = pd.read_csv(Path(run) / "postproc" / "gene_filtering" / "filtered_discovery.tsv", sep="\t")
    rm_path = Path(b0_dir) / "removed_units.tsv"
    removed = set()
    if rm_path.is_file():
        removed = set(map(tuple, pd.read_csv(rm_path, sep="\t")[["cycle", "caap_group", "Gene"]].astype(str).itertuples(index=False, name=None)))
    d = survivors(read_detail(b0_dir), removed)
    d = d.rename(columns={"Gene": "Gene", "Position": "Position"})
    key = ["Gene", "Position", "caap_group", "side"]
    res = compare_sets(_keyset(fd, key), _keyset(d, key))
    res["pass"] = res["only_observed"] == 0 and res["only_b_0"] == 0
    res["n_removed_units_b0"] = len(removed)
    return res


def checkpoint_C(run, b0_dir, tol):
    m = pd.read_csv(Path(run) / "ct_disambiguation" / "caas_convergence_master.csv",
                    usecols=["gene", "msa_pos", "caap_group", "side", "asr_path_score"], float_precision="round_trip")
    m = m.rename(columns={"gene": "Gene", "msa_pos": "Position"})
    d = read_detail(b0_dir)
    return compare_values(m, d[["Gene", "Position", "caap_group", "side", "asr_path_score"]],
                          ["Gene", "Position", "caap_group", "side"],
                          ["asr_path_score"], ["asr_path_score"], tol, bitwise=True)


def checkpoint_D(run, b0_dir, tol):
    ps = pd.read_csv(Path(run) / "scoring" / "position_scores.tsv", sep="\t", usecols=["Gene", "Position", "side", "CAAS_score"])
    pc = Path(b0_dir) / "perm_pos_cycle_caas.tsv.gz"
    d = pd.read_csv(pc, sep="\t", float_precision="round_trip")
    d = d.assign(caas_row=d["caas_score"])
    return compare_values(ps, d[["Gene", "Position", "side", "caas_row"]], ["Gene", "Position", "side"],
                          ["CAAS_score"], ["caas_row"], tol)


def checkpoint_E(run, b0_dir, tol):
    gs = pd.read_csv(Path(run) / "scoring" / "gene_scores.tsv", sep="\t",
                     usecols=["Gene", "gene_caas_score", "gene_caas_score_top_all", "gene_caas_score_bottom_all"])
    g0 = pd.read_csv(Path(b0_dir) / "gene_cycle_scores.tsv", sep="\t",
                     usecols=["Gene", "global_caas", "top_caas", "bottom_caas"])
    pairs = [("gene_caas_score", "global_caas"), ("gene_caas_score_top_all", "top_caas"), ("gene_caas_score_bottom_all", "bottom_caas")]
    # A gene with no scored position in a direction is NA on both sides.
    genes = sorted(set(gs["Gene"]) | set(g0["Gene"]))
    gs, g0 = gs.set_index("Gene").reindex(genes), g0.set_index("Gene").reindex(genes)
    bad, worst, ex = 0, 0.0, []
    for co, cb in pairs:
        x, y = gs[co].to_numpy(float), g0[cb].to_numpy(float)
        d = np.where(np.isnan(x) & np.isnan(y), 0.0, np.abs(x - y))  # one-sided NaN -> nan -> bad
        b = ~(d <= tol)
        bad += int(b.sum())
        worst = max(worst, float(np.max(np.where(np.isnan(d), np.inf, d))) if len(d) else 0.0)
        ex += [{"Gene": genes[i], "col": co, "observed": x[i], "b_0": y[i]} for i in np.flatnonzero(b)[:5]]
    return {"n_genes": len(genes), "n_value_mismatch": bad, "max_abs_delta": worst,
            "examples_value": ex[:5], "pass": bad == 0}


FUNCS = {"A": checkpoint_A, "B": checkpoint_B, "C": checkpoint_C, "D": checkpoint_D, "E": checkpoint_E}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="pipeline results dir (contains caastools/, scoring/, caas_permulation/ ...)")
    ap.add_argument("--b0-dir", help="default: <run>/caas_permulation/b0")
    ap.add_argument("--perm-disc", help="default: <run>/caas_permulation/perm_disc")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--checkpoints", default=CHECKPOINTS)
    ap.add_argument("--out", help="write the full report as JSON")
    a = ap.parse_args()
    b0_dir = a.b0_dir or str(Path(a.run) / "caas_permulation" / "b0")

    report, ok = {}, True
    for c in a.checkpoints.upper():
        try:
            r = FUNCS[c](a.run, b0_dir, a.tol, a.perm_disc) if c == "A" else FUNCS[c](a.run, b0_dir, a.tol)
        except FileNotFoundError as e:
            r = {"pass": False, "error": f"missing input: {e.filename}"}
        report[c] = r
        ok &= bool(r["pass"])
        head = {k: v for k, v in r.items() if not k.startswith("examples")}
        print(f"[{c}] {'PASS' if r['pass'] else 'FAIL'}  {json.dumps(head, default=str)}")
        if not r["pass"]:
            for k, v in r.items():
                if k.startswith("examples") and v:
                    print(f"     {k}: {v}")
    if a.out:
        Path(a.out).write_text(json.dumps(report, indent=2, default=str))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
