#!/usr/bin/env python3
# posenrich_prep_caas_null.py — Reduces the CAAS permulation null once, for all position-enrichment batches.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
PosenrichPrepCaasNull: parses the long-format CAAS permulation null once and stores, per
direction (global, top, bottom), one score per (position, cycle) in a compact pickle.

Every POSENRICH_RUN_BATCHED task uses the same null (only the gene-set files differ
between batches). Parsing a table of tens of millions of rows, building the position ID
and splitting it by direction in every task would repeat identical work, so it is done
here and each task loads the result (posenrich_enrich.py --caas-null-prepped).

The per-direction reduction is the one of posenrich_enrich.py (load_caas_cycle_null,
null_direction_subset): "top" and "bottom" keep that side, "global" keeps the maximum
over sides, and each (pos_id, cycle) pair has a single score.

Called by:  POSENRICH_PREP_NULL Nextflow process (posenrich.nf → posenrich_prep_caas_null.py)
Inputs:     --caas-cycle-null  perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, caas_score, n_schemes);
                               absent or a NO_FILE* sentinel gives an empty artifact
Outputs:    caas_null_prepped.pkl  dict with cycle_levels (sorted array of every cycle in the file)
                                   and global/top/bottom (DataFrames pos_id, cycle, score; pos_id categorical);
                                   all four entries are None for an empty artifact
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import os
import pickle

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Pre-reduce the CAAS permulation null for POSENRICH_RUN_BATCHED.")
    p.add_argument("--caas-cycle-null", required=True,
                   help="perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, "
                        "caas_score, n_schemes). Omit or pass a NO_FILE* sentinel to "
                        "write an empty artifact (p.perm stays NA downstream).")
    p.add_argument("--output", default="caas_null_prepped.pkl")
    return p.parse_args()


# ── Reduction ─────────────────────────────────────────────────────────────────


def write_empty(path):
    """Write the artifact that means "no CAAS null": every entry None (p-values stay NA downstream)."""
    with open(path, "wb") as fh:
        pickle.dump({"cycle_levels": None, "global": None, "top": None, "bottom": None}, fh)


def main():
    """Parse the null, reduce it per direction and pickle it."""
    args = parse_args()
    path = args.caas_cycle_null

    if not path or os.path.basename(path).startswith("NO_FILE") or not os.path.exists(path):
        write_empty(args.output)
        print(f"[posenrich_prep] no CAAS null supplied -> wrote empty {args.output}", flush=True)
        return

    # `cycle` is a label, not necessarily numeric (e.g. "b_1000"), so it is read as a
    # category rather than coerced to a number. caas_score is read exactly (round_trip);
    # only the repetitive string columns (Gene, side, cycle, pos_id) become categories,
    # which is where the memory goes.
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    if "caas_score" not in header:
        raise ValueError(f"{path} has no caas_score column: it predates the shared position score. "
                         "Regenerate the CAAS permulation null.")
    df = pd.read_csv(
        path, sep="\t",
        usecols=["Gene", "Position", "side", "cycle", "caas_score"],
        dtype={"Gene": "category", "side": "category", "cycle": "category"},
        float_precision="round_trip",
    )
    if df.empty:
        write_empty(args.output)
        print(f"[posenrich_prep] empty CAAS null file -> wrote empty {args.output}", flush=True)
        return

    # pos_id is cast to category at once, so the string column of every row collapses to
    # integer codes over the (much smaller) set of distinct positions.
    pos_id = (df["Gene"].astype(str) + ":" + df["Position"].astype(str)).astype("category")
    df["pos_id"] = pos_id
    df["score"] = df["caas_score"].fillna(0.0)
    # The categories of cycle are exactly its distinct values, so this reads the category
    # list instead of scanning every row.
    cycle_levels = np.sort(np.asarray(df["cycle"].cat.categories, dtype=object))

    long_df = df[["pos_id", "side", "cycle", "score"]]

    out = {"cycle_levels": cycle_levels}
    for direction in ("global", "top", "bottom"):
        if direction == "global":
            sub = long_df
        elif direction == "top":
            sub = long_df[long_df["side"] == "top"]
        else:
            sub = long_df[long_df["side"] == "bottom"]
        sub = (sub.groupby(["pos_id", "cycle"], observed=True, sort=False)["score"]
                  .max().reset_index())
        sub["pos_id"] = sub["pos_id"].astype("category")
        sub = sub.reset_index(drop=True)
        out[direction] = sub
        print(f"[posenrich_prep] {direction}: {len(sub)} (pos_id,cycle) null rows", flush=True)

    with open(args.output, "wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[posenrich_prep] wrote {args.output} ({len(cycle_levels)} cycles)", flush=True)


if __name__ == "__main__":
    main()
