#!/usr/bin/env python3
"""Two runs of the permulation null, table by table.

Compares ``<run>/caas_permulation`` (and its ``b0/`` slice) of run A with run B:

  gene_cycle_scores.tsv, perm_pos_sample.tsv, perm_pos_quantiles.tsv, removed_units.tsv, perm_pos_cycle_caas.tsv.gz
  caas_perms.rds                          (R's all.equal with the same tolerance; needs Rscript)
  --extra <path relative to the run>      further tables, e.g. scoring/position_scores.tsv

A table passes when its rows (identified by the identifier columns, with repeated identifiers numbered in order)
are the same set and every numeric column differs by at most ``--tol``; a table missing from one run fails, a table
missing from both is skipped. Exit status is 1 if any comparison fails.

Identifier columns are Gene/gene, Position/position, cycle, side, scheme and caap_group, plus every column that is
not numeric in both tables; the remaining columns are values.
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TABLES = ["gene_cycle_scores.tsv", "perm_pos_sample.tsv", "perm_pos_quantiles.tsv", "removed_units.tsv", "perm_pos_cycle_caas.tsv.gz"]
ID_COLUMNS = {"gene", "position", "cycle", "side", "scheme", "caap_group"}


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[])


def compare_table(path_a, path_b, tol):
    a, b = _read(path_a), _read(path_b)
    if list(a.columns) != list(b.columns):
        return {"pass": False, "columns_a": list(a.columns), "columns_b": list(b.columns)}
    values = []
    for c in a.columns:
        if c.lower() in ID_COLUMNS:
            continue
        na, nb = pd.to_numeric(a[c].replace({"NA": np.nan, "": np.nan}), errors="coerce"), pd.to_numeric(b[c].replace({"NA": np.nan, "": np.nan}), errors="coerce")
        if (na.notna() | a[c].isin(["NA", ""])).all() and (nb.notna() | b[c].isin(["NA", ""])).all():
            a[c], b[c] = na, nb
            values.append(c)
    key = [c for c in a.columns if c not in values]
    for df in (a, b):
        df["_n"] = df.groupby(key, dropna=False).cumcount() if key else range(len(df))
    cols = key + ["_n"]
    m = a.merge(b, on=cols, how="outer", suffixes=("_a", "_b"), indicator=True)
    out = {"n_a": len(a), "n_b": len(b), "only_a": int((m["_merge"] == "left_only").sum()), "only_b": int((m["_merge"] == "right_only").sum())}
    both = m[m["_merge"] == "both"]
    deltas, na_mismatch, bitwise = {}, 0, 0
    for c in values:
        x, y = both[f"{c}_a"].to_numpy(float), both[f"{c}_b"].to_numpy(float)
        na_mismatch += int((np.isnan(x) != np.isnan(y)).sum())
        ok = ~np.isnan(x) & ~np.isnan(y)
        d = np.abs(x[ok] - y[ok])
        deltas[c] = float(d.max()) if d.size else 0.0
        bitwise += int((x[ok] != y[ok]).sum())
    worst = max(deltas.values(), default=0.0)
    out.update({"max_abs_delta": worst, "worst_column": max(deltas, key=deltas.get) if deltas else None, "n_value_cells_not_bitwise_equal": bitwise,
                "na_mismatch": na_mismatch})
    out["pass"] = out["only_a"] == 0 and out["only_b"] == 0 and na_mismatch == 0 and worst <= tol
    return out


def compare_rds(path_a, path_b, tol, rscript):
    exe = shutil.which(rscript) or (rscript if Path(rscript).exists() else None)
    if exe is None:
        return {"pass": None, "skipped": f"{rscript} not found"}
    code = (f'a <- readRDS("{path_a}"); b <- readRDS("{path_b}"); r <- all.equal(a, b, tolerance = {tol}); '
            'if (!isTRUE(r)) { cat(head(r, 5), sep = "\\n"); quit(status = 1) }')
    p = subprocess.run([exe, "-e", code], capture_output=True, text=True)
    return {"pass": p.returncode == 0, "all_equal": "TRUE" if p.returncode == 0 else p.stdout.strip()[:400]}


def compare_runs(run_a, run_b, tol=1e-12, extra=(), rscript="Rscript"):
    run_a, run_b = Path(run_a), Path(run_b)
    results = {}
    for prefix in ("caas_permulation", "caas_permulation/b0"):
        for t in TABLES + (["caas_perms.rds"] if prefix == "caas_permulation" else []):
            results[f"{prefix}/{t}"] = _one(run_a / prefix / t, run_b / prefix / t, tol, rscript)
    for rel in extra:
        results[rel] = _one(run_a / rel, run_b / rel, tol, rscript)
    return results


def _one(a, b, tol, rscript):
    if not a.exists() and not b.exists():
        return {"pass": None, "skipped": "absent from both runs"}
    if a.exists() != b.exists():
        return {"pass": False, "missing_from": "A" if not a.exists() else "B"}
    return compare_rds(a, b, tol, rscript) if a.suffix == ".rds" else compare_table(a, b, tol)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="run directory A")
    ap.add_argument("--b", required=True, help="run directory B")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--extra", action="append", default=[], help="further table, relative to the run directory (repeatable)")
    ap.add_argument("--rscript", default="Rscript")
    args = ap.parse_args(argv)
    results = compare_runs(args.a, args.b, args.tol, args.extra, args.rscript)
    for name, r in results.items():
        status = "SKIP" if r["pass"] is None else ("PASS" if r["pass"] else "FAIL")
        print(f"[{status}] {name}  {json.dumps({k: v for k, v in r.items() if k != 'pass'})}")
    failed = [n for n, r in results.items() if r["pass"] is False]
    print(f"{len(failed)} failed, {sum(r['pass'] is True for r in results.values())} passed, {sum(r['pass'] is None for r in results.values())} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
