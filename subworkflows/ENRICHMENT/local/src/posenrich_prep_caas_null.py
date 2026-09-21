#!/usr/bin/env python3
# =============================================================================
# posenrich_prep_caas_null.py - one-time reduction of the CAAS permulation
# null for POSENRICH_RUN_BATCHED
# =============================================================================
# perm_pos_cycle_caas.tsv.gz is broadcast identically to every
# POSENRICH_RUN_BATCHED task (same CAAS permulation null; only the GMT term
# sets differ per batch). Parsing the raw long-format file independently in
# every batch task means re-reading a tens-of-millions-row table, rebuilding
# a string pos_id column, and re-deriving the same per-direction
# (global/top/bottom) subsets once per batch, for no benefit -- none of that
# depends on which GMTs a batch is testing.
#
# This script does that parse + derivation exactly ONCE and writes a compact
# pickle (categorical Gene/pos_id, int32 cycle, float32 score) that every
# batch task then just loads. The per-direction subsetting mirrors
# posenrich_enrich.py's load_caas_cycle_null()/null_direction_subset() exactly
# (same side-handling, same drop_duplicates(subset=["pos_id","cycle"],
# keep="last") resolution) -- this is a reduction of redundant work, not a
# change in what is computed.
# =============================================================================

import argparse
import os
import pickle

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Pre-reduce the CAAS permulation null for POSENRICH_RUN_BATCHED.")
    p.add_argument("--caas-cycle-null", required=True,
                   help="perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, "
                        "caas_sum, n_schemes). Omit or pass a NO_FILE* sentinel to "
                        "write an empty artifact (p.perm stays NA downstream).")
    p.add_argument("--output", default="caas_null_prepped.pkl")
    return p.parse_args()


def write_empty(path):
    with open(path, "wb") as fh:
        pickle.dump({"cycle_levels": None, "global": None, "top": None, "bottom": None}, fh)


def main():
    args = parse_args()
    path = args.caas_cycle_null

    if not path or os.path.basename(path).startswith("NO_FILE") or not os.path.exists(path):
        write_empty(args.output)
        print(f"[posenrich_prep] no CAAS null supplied -> wrote empty {args.output}", flush=True)
        return

    df = pd.read_csv(
        path, sep="\t",
        usecols=["Gene", "Position", "side", "cycle", "caas_sum", "n_schemes"],
        dtype={"Gene": "category", "side": "category", "cycle": "int32",
               "caas_sum": "float32", "n_schemes": "int32"},
    )
    if df.empty:
        write_empty(args.output)
        print(f"[posenrich_prep] empty CAAS null file -> wrote empty {args.output}", flush=True)
        return

    # pos_id is built once here instead of once per batch task; immediately
    # cast to category so the ~27M-row string column collapses to integer
    # codes over its (much smaller) set of distinct positions.
    pos_id = (df["Gene"].astype(str) + ":" + df["Position"].astype(str)).astype("category")
    df["pos_id"] = pos_id
    df["score"] = (df["caas_sum"] / df["n_schemes"].replace(0, np.nan)).fillna(0.0).astype("float32")
    cycle_levels = np.sort(df["cycle"].unique())

    long_df = df[["pos_id", "side", "cycle", "score"]]

    out = {"cycle_levels": cycle_levels}
    for direction in ("global", "top", "bottom"):
        if direction == "global":
            sub = long_df
        elif direction == "top":
            sub = long_df[long_df["side"] == "top"]
        else:
            sub = long_df[long_df["side"] == "bottom"]
        sub = sub.drop_duplicates(subset=["pos_id", "cycle"], keep="last")[["pos_id", "cycle", "score"]].copy()
        sub["pos_id"] = sub["pos_id"].astype("category")
        sub = sub.reset_index(drop=True)
        out[direction] = sub
        print(f"[posenrich_prep] {direction}: {len(sub)} (pos_id,cycle) null rows", flush=True)

    with open(args.output, "wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[posenrich_prep] wrote {args.output} ({len(cycle_levels)} cycles)", flush=True)


if __name__ == "__main__":
    main()
