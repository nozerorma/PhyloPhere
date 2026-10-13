#!/usr/bin/env python3
# test_accum_occupancy_universe.py: the occupancy nulls of accumulation ('naive', 'cons_decile') draw from the positional
# universe U when the filtered permulation null is given, and from the tested pool otherwise.
#   U = (tested positions of the cleaned-background genes with caas_score > 0 in >= 1 null cycle) union (observed positions)
# Every draw is recorded through a wrapper of the random generator, so the checks look at the draws themselves.
# Run: python3 tests/null_universe/test_accum_occupancy_universe.py   (scratch files go to ~/.cc/occupancy_universe_test, removed afterwards)

import argparse
import gzip
import importlib.util
import json
import logging
import shutil
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RANDOMIZE = HERE.parents[1] / "subworkflows" / "CT_ACCUMULATION" / "local" / "src" / "randomization" / "randomize.py"
WORK = Path.home() / ".cc" / "occupancy_universe_test"

# pyarrow is only needed for the parquet outputs, which are captured below instead of written
for name in ("pyarrow", "pyarrow.parquet"):
    sys.modules.setdefault(name, types.ModuleType(name))
spec = importlib.util.spec_from_file_location("randomize", RANDOMIZE)
rz = importlib.util.module_from_spec(spec)
sys.modules["randomize"] = rz
spec.loader.exec_module(rz)
logging.getLogger().setLevel(logging.ERROR)

GENES = ["G1", "G2", "G3", "G4"]
N_POS = 20


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"ok   {msg}")


class RecRng:
    """Generator wrapper that appends every choice() call to a jsonl file."""
    def __init__(self, rng, path):
        self._rng, self._path = rng, path

    def choice(self, a, size, replace):
        out = self._rng.choice(a, size=size, replace=replace)
        with open(self._path, "a") as fh:
            fh.write(json.dumps({"elig": sorted(int(x) for x in a), "size": int(size), "replace": bool(replace),
                                 "out": [int(x) for x in out]}) + "\n")
        return out


def write_global(path, cons):
    """global table: 4 genes x 20 msa positions; `position` is the running row index; cons maps (gene, msa_pos) to cons_idx."""
    rows, pos = [], 0
    for g in GENES:
        for m in range(N_POS):
            rows.append(dict(gene=g, msa_pos=m, position=pos, cons_idx=cons(g, m), masked=False, iscaas=False))
            pos += 1
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df


def write_null(path, rows):
    """perm_pos_cycle_caas.tsv.gz rows: (Gene, Position, side, cycle, caas_score)."""
    with gzip.open(path, "wt") as fh:
        fh.write("Gene\tPosition\tside\tcycle\tcaas_score\n")
        for g, p, s, c, sc in rows:
            fh.write(f"{g}\t{p}\t{s}\t{c}\t{sc}\n")


def run(out, rtype, null=None, observed=None, tested=None, bg_genes=None, decile_bins=None, change_side="both",
        n_rands=40, cons=lambda g, m: ((GENES.index(g) * N_POS + m) * 7 % 20) / 20 + 0.01):
    """Run main() with every draw recorded. Returns (global df, draws, captured parquet tables, log lines)."""
    out.mkdir(parents=True, exist_ok=True)
    gdf = write_global(out / "global.csv", cons)
    obs = observed if observed is not None else []
    pd.DataFrame(obs, columns=["Gene", "Position", "side", "caap_group"]).to_csv(out / "fd.tsv", sep="\t", index=False)
    if tested is not None:
        (out / "bg.output").write_text("".join(f"{g}\t{','.join(map(str, p))}\n" for g, p in tested.items()))
    (out / "bg_genes.txt").write_text("\n".join(bg_genes or GENES) + "\n")
    if null is not None:
        write_null(out / "null.tsv.gz", null)
    rec = out / "draws.jsonl"
    rec.unlink(missing_ok=True)
    args = argparse.Namespace(
        global_csv=str(out / "global.csv"), caas_csv=str(out / "fd.tsv"),
        background_positions=str(out / "bg.output") if tested is not None else None,
        bg_caas=str(out / "bg_genes.txt") if bg_genes is not None else None,
        output_prefix=str(out / "acc"), randomization_type=rtype,
        null_cycle_caas=str(out / "null.tsv.gz") if null is not None else None, cycle_roster=None,
        n_randomizations=n_rands, workers=1, compress=False, export_individual_rand=False, decile_bins=decile_bins,
        global_seed=3, precompute_masks=True, change_side=change_side, log_level="INFO")
    tables, lines = {}, []

    class H(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())
    h, root = H(), logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.INFO)
    real_rng, real_pq = np.random.default_rng, pd.DataFrame.to_parquet
    np.random.default_rng = lambda seed=None: RecRng(real_rng(seed), rec)
    pd.DataFrame.to_parquet = lambda self, path, **kw: tables.__setitem__(Path(path).name, self.copy())
    try:
        rz.main(args)
    finally:
        np.random.default_rng, pd.DataFrame.to_parquet = real_rng, real_pq
        root.removeHandler(h)
        root.setLevel(logging.ERROR)
    draws = [json.loads(x) for x in rec.read_text().splitlines()] if rec.exists() else []
    # The table the draws index: main() left-joins the CAAS rows on the position table, so a position with k detections
    # (sides, schemes) repeats k times; the draw pool keeps the first row of each position.
    fd = pd.DataFrame(obs, columns=["gene", "msa_pos", "side", "caap_group"]).drop_duplicates()
    mdf = pd.merge(gdf, fd, how="left", on=["gene", "msa_pos"])
    return mdf, draws, tables, lines


def rows_of(gdf, pairs):
    """Pool row (first row of the position) of each (gene, msa_pos) pair; gdf is the merged table returned by run()."""
    key = {}
    for i, r in enumerate(gdf.itertuples()):
        key.setdefault((r.gene, r.msa_pos), i)
    return {key[p] for p in pairs}


# Fixture. Tested: positions 0-14 of G1-G3; G4 is not in the cleaned background.
TESTED = {g: list(range(15)) for g in GENES}
BG_GENES = ["G1", "G2", "G3"]
# Null-positive (caas_score > 0) positions: G1 0-5 (and 17, untested), G2 2-9, G3 10-11, G4 3 (gene outside the
# background); G1:12 has score 0 and G2:12 an empty score, so they are not positive.
NULL = ([("G1", p, "top", "b_1", 0.5) for p in range(6)] + [("G1", 2, "bottom", "b_2", 0.4), ("G1", 17, "top", "b_2", 0.4)]
        + [("G2", p, "bottom", "b_1", 0.3) for p in range(2, 10)] + [("G3", 10, "top", "b_3", 0.2), ("G3", 11, "top", "b_3", 0.2)]
        + [("G4", 3, "top", "b_1", 0.5), ("G1", 12, "top", "b_1", 0), ("G2", 12, "top", "b_1", "")])
# Observed: G1:0 and G1:1 (null-positive), G2:3 (null-positive), G3:13 (tested, never null-positive), G2:18 (untested), a
# GS-only position G3:14 (tested, never null-positive), a side none row (not counted) and the two sides of G1:3.
OBS = [("G1", 0, "top", "US"), ("G1", 1, "bottom", "US"), ("G2", 3, "top", "US"), ("G3", 13, "top", "US"),
       ("G2", 18, "bottom", "US"), ("G3", 14, "top", "GS1"), ("G4", 5, "none", "US"),
       ("G1", 3, "top", "US"), ("G1", 3, "bottom", "US")]
OBS_COUNTED = [o for o in OBS if o[3] == "US" and o[2] in ("top", "bottom")]       # 7 rows under direction both


def expected_universe():
    tested = {(g, p) for g in BG_GENES for p in TESTED[g]}
    null_pos = {(g, p) for g, p, s, c, sc in NULL if sc not in (0, "")}
    observed = {(g, p) for g, p, s, c in OBS if s in ("top", "bottom")}
    return (tested & null_pos) | observed


def test_draws_stay_in_universe():
    U = expected_universe()
    check(("G1", 17) not in U and ("G4", 3) not in U and ("G1", 12) not in U and ("G2", 12) not in U,
          "fixture: untested, outside-gene and non-positive null rows are not in U")
    check(("G2", 18) in U and ("G3", 13) in U and ("G3", 14) in U, "fixture: observed positions stay in U although the null never scored them or they are untested")
    for rtype in ("naive", "cons_decile"):
        gdf, draws, tables, lines = run(WORK / f"u_{rtype}", rtype, null=NULL, observed=OBS, tested=TESTED, bg_genes=BG_GENES)
        u_rows = rows_of(gdf, U)
        check(len(draws) > 0, f"{rtype}: draws were recorded")
        check(all(set(d["elig"]) <= u_rows for d in draws), f"{rtype}: every eligible pool is a subset of U")
        check(all(set(d["out"]) <= u_rows for d in draws), f"{rtype}: no drawn position leaves U")
        check(any("Draw pool" in x and ": U (" in x for x in lines), f"{rtype}: the log names U as the pool")
        if rtype == "naive":
            check(all(set(d["elig"]) == u_rows for d in draws), "naive: the pool is exactly U")
            check(all(d["size"] == len(OBS_COUNTED) for d in draws), f"naive: each draw redraws the {len(OBS_COUNTED)} observed rows")


def test_cons_decile_counts():
    """The number drawn in each decile equals the observed rows of that decile, and each draw stays inside its decile."""
    U = expected_universe()
    gdf, draws, tables, lines = run(WORK / "dec", "cons_decile", null=NULL, observed=OBS, tested=TESTED, bg_genes=BG_GENES)
    bins = np.percentile(gdf["cons_idx"], np.arange(0, 101, 10))
    dec = rz.decile_keys(gdf["cons_idx"].to_numpy(np.float32), bins)   # the workers read cons_idx as float32
    obs_rows = [rows_of(gdf, [(g, p)]).pop() for g, p, s, c in OBS_COUNTED]
    expect = pd.Series(dec[obs_rows]).value_counts().to_dict()
    got = {}
    for d in draws:
        ds = set(dec[d["elig"]])
        assert len(ds) == 1, "a pool spans several deciles"
        got.setdefault(int(ds.pop()), set()).add(d["size"])
    check(True, "cons_decile: each pool lies in one decile")
    check({k: next(iter(v)) for k, v in got.items()} == {int(k): int(v) for k, v in expect.items()},
          f"cons_decile: drawn counts per decile {({k: sorted(v) for k, v in got.items()})} == observed per decile {expect}")
    check(all(len(v) == 1 for v in got.values()), "cons_decile: the count of a decile is the same in every draw")
    check(all(set(d["elig"]) <= rows_of(gdf, U) for d in draws), "cons_decile: pools within U")
    g = tables["acc_cons_decile_global.parquet"]
    check(int(g["NumObserved"].sum()) == len(OBS_COUNTED), "diagnostic: NumObserved sums to the observed rows")
    check(int(g["NumPool"].sum()) == len(U), f"diagnostic: NumPool sums to |U| = {len(U)}")
    pg = tables["acc_cons_decile_pool_per_gene.parquet"].set_index("gene")
    per_gene = {g_: sum(1 for (gg, p) in U if gg == g_) for g_ in BG_GENES}
    check({k: int(pg.loc[k].sum()) for k in per_gene} == per_gene and "G4" not in pg.index,
          "diagnostic: per-gene pool occupancy equals U per gene (G4 absent)")
    check(all(isinstance(c, str) for t_ in tables.values() for c in t_.columns), "diagnostic: all column names are strings (parquet)")
    check("acc_cons_decile_per_gene.parquet" in tables and "NumPositions" in g.columns,
          "diagnostic: the previous tables and columns are still written")


def test_fallback_pool():
    """Without a null table the pool is the tested positions of the cleaned-background genes."""
    old = {(g, p) for g in BG_GENES for p in TESTED[g]}
    for rtype in ("naive", "cons_decile"):
        gdf, draws, tables, lines = run(WORK / f"f_{rtype}", rtype, null=None, observed=OBS, tested=TESTED, bg_genes=BG_GENES)
        check(all(set(d["elig"]) <= rows_of(gdf, old) for d in draws) and len(draws) > 0,
              f"{rtype}, no null: draws stay in the tested pool")
        if rtype == "naive":
            check(all(set(d["elig"]) == rows_of(gdf, old) for d in draws), "naive, no null: the pool is the 45 tested positions")
        check(any("tested pool (no null table given" in x for x in lines), f"{rtype}, no null: the log names the fallback pool")
    # a null without any positive score is a null without cycles: same fallback
    empty = [("G1", 0, "top", "b_1", 0)]
    gdf, draws, tables, lines = run(WORK / "f_empty", "naive", null=empty, observed=OBS, tested=TESTED, bg_genes=BG_GENES)
    check(all(set(d["elig"]) == rows_of(gdf, old) for d in draws), "null without positive scores: the old pool is used")
    check(any("previous pool is used" in x for x in lines), "null without positive scores: the log says so")


def test_replacement_edge():
    """A key with fewer pool positions than observed rows is drawn with replacement, logged and flagged in the table."""
    # U = {G1:0, G1:1} only; G1:1 is observed on both sides (2 rows, 1 position)
    null = [("G1", 0, "top", "b_1", 0.5), ("G1", 1, "top", "b_1", 0.5)]
    obs = [("G1", 0, "top", "US"), ("G1", 1, "top", "US"), ("G1", 1, "bottom", "US")]
    cons = lambda g, m: 0.9 if (g, m) == ("G1", 1) else 0.1
    gdf, draws, tables, lines = run(WORK / "edge", "cons_decile", null=null, observed=obs, tested=TESTED, bg_genes=BG_GENES,
                                    decile_bins="0,0.5,1", cons=cons)
    sizes = {tuple(d["elig"]): (d["size"], d["replace"]) for d in draws}
    hi = tuple(sorted(rows_of(gdf, [("G1", 1)])))
    lo = tuple(sorted(rows_of(gdf, [("G1", 0)])))
    check(sizes[hi] == (2, True), "the decile with 1 pool position and 2 observed rows is drawn with replacement")
    check(sizes[lo] == (1, False), "the other decile is drawn without replacement")
    check(all(set(d["out"]) <= set(d["elig"]) for d in draws), "draws with replacement stay inside the pool")
    check(any("1 of 2 draw keys hold fewer pool positions" in x for x in lines), "the log counts the deciles that hit the branch")
    g = tables["acc_cons_decile_global.parquet"].set_index("decile")
    check(g["Replace"].to_dict() == {1: False, 2: True} and g["NumPool"].to_dict() == {1: 1, 2: 1} and g["NumObserved"].to_dict() == {1: 1, 2: 2},
          f"the decile table flags it: {g[['NumPool', 'NumObserved', 'Replace']].to_dict('index')}")
    # naive: the whole pool (2 positions) is smaller than the 3 observed rows
    gdf, draws, tables, lines = run(WORK / "edge_naive", "naive", null=null, observed=obs, tested=TESTED, bg_genes=BG_GENES)
    check(all(d["replace"] and d["size"] == 3 and len(d["elig"]) == 2 for d in draws), "naive: a pool of 2 positions with 3 observed rows is drawn with replacement")
    check(any("1 of 1 draw keys hold fewer pool positions" in x for x in lines), "naive: logged")


if __name__ == "__main__":
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    try:
        test_draws_stay_in_universe()
        test_cons_decile_counts()
        test_fallback_pool()
        test_replacement_edge()
        print("all occupancy-universe tests passed")
    finally:
        shutil.rmtree(WORK, ignore_errors=True)
