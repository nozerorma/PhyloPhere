#!/usr/bin/env python3
"""The observed contract files of two runs, file by file.

Compares run A (the baseline, made by the former observed chain) with run B (the b_0 slice of the core):

  caastools/discovery.tab                      header exact; rows as a set (a row is its 19 columns, `ms` as a set of species);
                                               B ordered by gene
  caastools/background.output                  gene -> tested positions; B ordered by gene
  caastools/background_genes.output            set of genes; B ordered
  meta_caas/meta_caas/<name>_meta_caas.tsv     rows without `tag` as a set; every tag of B is the content id of its row
                                               (recomputed with core.meta.caas_id), so ids are compared by the set of
                                               positions they label, never by value
  ct_disambiguation/caas_convergence_master.csv  rows by (gene, msa_pos, caap_group, side); numeric columns within
                                               --tol; `tag_support` by shape (the same tally of counts) and with its ids
                                               among B's meta ids at that position and scheme; the modal-residue columns
                                               (domain_N_{anc,top,bot}_aa) may differ only where two residues tie for the
                                               maximum support, and those cells are counted; `derived_agreement` and
                                               `convergence_type` are functions of the modal derived residues (fop_pool), so
                                               they may differ only in a row where a top/bot residue cell was tolerated as a
                                               tie (or, when B has `agreement_ambiguous`, in a row B flags True, which
                                               is the exact criterion), and those cells are counted too; a column of NEW_IN_B that only B has
                                               (`agreement_ambiguous`) is left out of the comparison and listed

A file missing from one run fails; a file missing from both is skipped. Exit status is 1 if any comparison fails.
The order of the rows inside a position is not compared: the b_0 export lists entries by position, trait file name and
scheme, the former discovery by file-system order.
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.meta import caas_id  # noqa: E402

DISCOVERY = "caastools/discovery.tab"
BACKGROUND = "caastools/background.output"
BACKGROUND_GENES = "caastools/background_genes.output"
META_DIR = "meta_caas/meta_caas"
MASTER = "ct_disambiguation/caas_convergence_master.csv"
MASTER_KEY = ["gene", "msa_pos", "caap_group", "side"]
_MODAL = re.compile(r"domain_\d+_(anc|top|bot)_aa")
_DERIVED_MODAL = re.compile(r"domain_\d+_(top|bot)_aa")   # the residues agree_num / convergence_type are computed from
DERIVED = ("derived_agreement", "convergence_type")
NEW_IN_B = ("agreement_ambiguous",)   # master columns a baseline made before they existed may lack
_ID = re.compile(r"CAAS_[0-9A-F]{16}")
EXAMPLES = 3


def _tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[], quoting=csv.QUOTE_NONE)


def _row_hashes(df):
    return pd.util.hash_pandas_object(df.reset_index(drop=True), index=False).to_numpy()


def _multiset_diff(a, b):
    """Rows of a not in b and of b not in a, as multisets: (n_only_a, n_only_b, example rows of each)."""
    ca, cb = Counter(_row_hashes(a)), Counter(_row_hashes(b))
    only_a, only_b = ca - cb, cb - ca

    def examples(df, extra):
        h = pd.Series(_row_hashes(df))
        return df[h.isin(set(extra))].head(EXAMPLES).astype(str).values.tolist()

    return sum(only_a.values()), sum(only_b.values()), examples(a, only_a), examples(b, only_b)


def _ordered(values):
    """True when the sequence is in byte order (what `LC_ALL=C sort` gives)."""
    values = list(values)
    return values == sorted(values)


def compare_discovery(pa, pb):
    a, b = _tsv(pa), _tsv(pb)
    if list(a.columns) != list(b.columns):
        return {"pass": False, "columns_a": list(a.columns), "columns_b": list(b.columns)}
    for df in (a, b):
        if "ms" in df:
            df["ms"] = df["ms"].map(lambda s: ",".join(sorted(s.split(","))) if s else s)
    oa, ob, ea, eb = _multiset_diff(a, b)
    sorted_b = _ordered(b["gene"])
    return {"n_a": len(a), "n_b": len(b), "only_a": oa, "only_b": ob, "examples_only_a": ea, "examples_only_b": eb,
            "b_ordered_by_gene": sorted_b, "pass": oa == 0 and ob == 0 and sorted_b}


def _background(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    return {l.split("\t", 1)[0]: l.split("\t", 1)[1] if "\t" in l else "" for l in lines}, [l.split("\t", 1)[0] for l in lines]


def compare_background(pa, pb):
    da, _ = _background(pa)
    db, order_b = _background(pb)
    differ = sorted(g for g in set(da) | set(db) if da.get(g) != db.get(g))
    return {"n_a": len(da), "n_b": len(db), "genes_that_differ": len(differ), "examples": differ[:EXAMPLES],
            "b_ordered_by_gene": _ordered(order_b), "pass": not differ and _ordered(order_b)}


def compare_genes(pa, pb):
    ga, gb = [l.strip() for l in open(pa) if l.strip()], [l.strip() for l in open(pb) if l.strip()]
    only_a, only_b = sorted(set(ga) - set(gb)), sorted(set(gb) - set(ga))
    return {"n_a": len(ga), "n_b": len(gb), "only_a": len(only_a), "only_b": len(only_b), "examples": (only_a + only_b)[:EXAMPLES],
            "b_ordered": _ordered(gb), "pass": not only_a and not only_b and _ordered(gb)}


def _na(v):
    return "" if v == "NA" else v


def compare_meta(pa, pb):
    a, b = _tsv(pa), _tsv(pb)
    if list(a.columns) != list(b.columns):
        return {"pass": False, "columns_a": list(a.columns), "columns_b": list(b.columns)}
    rows = [c for c in a.columns if c != "tag"]
    oa, ob, ea, eb = _multiset_diff(a[rows], b[rows])
    bad = []
    for r in b.itertuples(index=False):
        d = r._asdict()
        expect = caas_id(d["Gene"], d["Position"], _na(d.get("trait", "")), d["caap_group"], _na(d["caas"]), _na(d["amino_encoded"]), d["pattern"])
        if d["tag"] != expect:
            bad.append((d["Gene"], d["Position"], d["tag"], expect))
    ids_ok = bool(b["tag"].map(lambda t: bool(_ID.fullmatch(t))).all()) if len(b) else True
    return {"n_a": len(a), "n_b": len(b), "only_a": oa, "only_b": ob, "examples_only_a": ea, "examples_only_b": eb,
            "ids_checked": len(b), "ids_not_the_content_id": len(bad), "id_examples": bad[:EXAMPLES], "ids_well_formed": ids_ok,
            "pass": oa == 0 and ob == 0 and not bad and ids_ok}


def _tally(cell):
    """'A:3,B:1,' -> {'A': 3, 'B': 1}"""
    out = {}
    for tok in str(cell).split(","):
        if ":" in tok:
            k, n = tok.rsplit(":", 1)
            out[k] = int(n)
    return out


def _numeric(col_a, col_b):
    def conv(s):
        return pd.to_numeric(s.replace({"": np.nan, "NA": np.nan}), errors="coerce")
    na, nb = conv(col_a), conv(col_b)
    ok_a = na.notna() | col_a.isin(["", "NA"])
    ok_b = nb.notna() | col_b.isin(["", "NA"])
    return (na, nb) if ok_a.all() and ok_b.all() and (na.notna().any() or nb.notna().any()) else None


def _id_tallies_shape(cell):
    return sorted(_tally(cell).values())


def compare_master(pa, pb, tol, meta_b=None):
    a = pd.read_csv(pa, dtype=str, keep_default_na=False, na_values=[])
    b = pd.read_csv(pb, dtype=str, keep_default_na=False, na_values=[])
    added = [c for c in b.columns if c not in a.columns and c in NEW_IN_B]
    has_flag = "agreement_ambiguous" in added
    flag_b = b["agreement_ambiguous"].copy() if has_flag else None
    b = b.drop(columns=added)
    if list(a.columns) != list(b.columns):
        return {"pass": False, "columns_a": list(a.columns), "columns_b": list(b.columns)}
    if has_flag:
        b["_flag_b"] = flag_b                          # kept apart: the comparison below runs over A's columns
    for df in (a, b):
        df["_n"] = df.groupby(MASTER_KEY).cumcount()
    m = a.merge(b, on=MASTER_KEY + ["_n"], how="outer", suffixes=("_a", "_b"), indicator=True)
    only_a, only_b = int((m["_merge"] == "left_only").sum()), int((m["_merge"] == "right_only").sum())
    both = m[m["_merge"] == "both"]
    deltas, na_mismatch, other_diff, tie_cells, not_tie = {}, 0, {}, 0, []
    tag_shape_bad = 0
    # first the modal-residue cells: the rows holding a tolerated tie in a derived residue explain the derived columns
    tie_rows = set()
    for c in a.columns:
        if not _MODAL.fullmatch(c):
            continue
        for i in both.index[both[f"{c}_a"] != both[f"{c}_b"]]:
            ta, tb = _tally(both.at[i, f"{c}_support_a"]), _tally(both.at[i, f"{c}_support_b"])
            va, vb = both.at[i, f"{c}_a"], both.at[i, f"{c}_b"]
            if ta and tb and ta.get(va) == max(ta.values()) and tb.get(vb) == max(tb.values()):
                tie_cells += 1
                if _DERIVED_MODAL.fullmatch(c):
                    tie_rows.add(i)
            else:
                not_tie.append((c, *(both.loc[i, MASTER_KEY].tolist()), va, vb))
    if has_flag:
        # B says itself which rows rest on a tie of the encoded derived residue, the one agreement is computed from;
        # the residue cell shown in the master is only a proxy for it.
        tie_rows = set(both.index[both["_flag_b"] == "True"])
    derived_tied, derived_unexplained = {}, []

    def derived_differs(c, rows):
        for i in rows:
            if i in tie_rows:
                derived_tied[c] = derived_tied.get(c, 0) + 1
            else:
                derived_unexplained.append((c, *both.loc[i, MASTER_KEY].tolist()))
    for c in a.columns:
        if c in MASTER_KEY or c == "_n" or _MODAL.fullmatch(c):
            continue
        ca, cb = both[f"{c}_a"], both[f"{c}_b"]
        if c == "tag_support":
            tag_shape_bad += int((ca.map(_id_tallies_shape) != cb.map(_id_tallies_shape)).sum())
            continue
        num = _numeric(ca, cb)
        if num is not None:
            x, y = num[0].to_numpy(float), num[1].to_numpy(float)
            na_mismatch += int((np.isnan(x) != np.isnan(y)).sum())
            ok = ~np.isnan(x) & ~np.isnan(y)
            if c in DERIVED:
                derived_differs(c, both.index[ok & (np.abs(np.where(ok, x - y, 0.0)) > tol)])
                continue
            deltas[c] = float(np.abs(x[ok] - y[ok]).max()) if ok.any() else 0.0
            continue
        diff = ca != cb
        if not diff.any():
            continue
        if c in DERIVED:
            derived_differs(c, both.index[diff])
        else:
            other_diff[c] = int(diff.sum())
    worst = max(deltas.values(), default=0.0)
    ids_outside = None
    if meta_b is not None and Path(meta_b).exists():
        meta = _tsv(meta_b)
        known = {(int(p), g): set() for p, g in zip(meta["Position"], meta["caap_group"])}
        for t, p, g in zip(meta["tag"], meta["Position"], meta["caap_group"]):
            known[(int(p), g)].add(t)
        ids_outside = 0
        for r in b.itertuples(index=False):
            d = r._asdict()
            for tok in _tally(d["tag_support"]):
                ids_outside += tok not in known.get((int(d["msa_pos"]), d["caap_group"]), ())
    out = {"columns_added_in_b": added, "n_a": len(a), "n_b": len(b), "only_a": only_a, "only_b": only_b, "max_abs_delta": worst,
           "worst_column": max(deltas, key=deltas.get) if deltas else None, "na_mismatch": na_mismatch,
           "columns_with_other_differences": other_diff, "modal_residue_tie_cells_tolerated": tie_cells,
           "modal_residue_cells_differing_without_a_tie": len(not_tie), "examples_not_tie": not_tie[:EXAMPLES],
           "rows_flagged_ambiguous_in_b": len(tie_rows) if has_flag else None, "derived_cells_explained_by_a_tie": derived_tied,
           "derived_cells_differing_without_a_tie": len(derived_unexplained), "examples_derived_without_a_tie": derived_unexplained[:EXAMPLES],
           "tag_support_shape_differs": tag_shape_bad, "tag_support_ids_outside_meta_b": ids_outside}
    out["pass"] = (only_a == 0 and only_b == 0 and na_mismatch == 0 and worst <= tol and not other_diff and not not_tie and not derived_unexplained
                   and tag_shape_bad == 0 and not ids_outside)
    return out


def compare_runs(run_a, run_b, tol=1e-12):
    run_a, run_b = Path(run_a), Path(run_b)
    results = {}

    def one(rel, fn, *extra):
        a, b = run_a / rel, run_b / rel
        if not a.exists() and not b.exists():
            results[rel] = {"pass": None, "skipped": "absent from both runs"}
        elif a.exists() != b.exists():
            results[rel] = {"pass": False, "missing_from": "A" if not a.exists() else "B"}
        else:
            results[rel] = fn(a, b, *extra)

    one(DISCOVERY, compare_discovery)
    one(BACKGROUND, compare_background)
    one(BACKGROUND_GENES, compare_genes)
    names = {f.name for d in (run_a, run_b) for f in (d / META_DIR).glob("*_meta_caas.tsv")} if (run_a / META_DIR).exists() or (run_b / META_DIR).exists() else set()
    for name in sorted(names):
        one(f"{META_DIR}/{name}", compare_meta)
    one(MASTER, compare_master, tol, run_b / META_DIR / "global_meta_caas.tsv")
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="baseline run directory (the former observed chain)")
    ap.add_argument("--b", required=True, help="run directory of the b_0 slice")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--report", help="write the results as JSON")
    args = ap.parse_args(argv)
    results = compare_runs(args.a, args.b, args.tol)
    for name, r in results.items():
        status = "SKIP" if r["pass"] is None else ("PASS" if r["pass"] else "FAIL")
        print(f"[{status}] {name}  {json.dumps({k: v for k, v in r.items() if k != 'pass'}, default=str)}")
    if args.report:
        Path(args.report).write_text(json.dumps(results, indent=1, default=str))
    failed = [n for n, r in results.items() if r["pass"] is False]
    print(f"{len(failed)} failed, {sum(r['pass'] is True for r in results.values())} passed, {sum(r['pass'] is None for r in results.values())} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
