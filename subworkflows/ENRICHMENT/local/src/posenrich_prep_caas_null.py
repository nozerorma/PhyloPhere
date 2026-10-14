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
Inputs:     --caas-cycle-null  perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, caas_score, score_aggregation);
                               absent or a NO_FILE* sentinel gives an empty artifact
            --cycle-roster     cycle_roster.txt of the same CAAS_CORE_MERGE run (optional): N of the null, the
                               cycles with no detection included
Outputs:    caas_null_prepped.pkl  dict with cycle_levels (sorted array of the replay roster, or of the cycles
                                   present in the file when no roster fits)
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
                        "caas_score, score_aggregation). Omit or pass a NO_FILE* sentinel to "
                        "write an empty artifact (p.perm stays NA downstream).")
    p.add_argument("--cycle-roster", default=None,
                   help="cycle_roster.txt (one replayed cycle per line): N of the null. Without it, or when it "
                        "does not list every cycle of the null file, N is the number of cycles present in the file.")
    p.add_argument("--output", default="caas_null_prepped.pkl")
    return p.parse_args()


# ── Reduction ─────────────────────────────────────────────────────────────────


def resolve_cycle_levels(roster_path, cycles_present):
    """Sorted array of the cycles that make N of the null: the roster when it lists every cycle present, else the cycles present.

    The roster counts the cycles that detected nothing, which the null file cannot show; falling back to the
    cycles present can only raise the permulation p-values. The choice is printed either way. The same rule as
    scoring_compute.R and randomize.py (accumulation).
    """
    present = {str(c) for c in cycles_present}
    if roster_path and os.path.exists(roster_path) and not os.path.basename(roster_path).startswith("NO_"):
        with open(roster_path) as fh:
            roster = {line.strip() for line in fh if line.strip()}
        stray = present - roster
        if not stray:
            print(f"[posenrich_prep] N = {len(roster)} cycles (replay roster; {len(present)} present in the null file)",
                  flush=True)
            return np.sort(np.asarray(sorted(roster), dtype=object))
        print(f"[posenrich_prep] WARNING: the roster does not list {len(stray)} cycle(s) of the null file "
              f"(e.g. {sorted(stray)[:3]}): ignored, N = cycles present ({len(present)})", flush=True)
    else:
        print(f"[posenrich_prep] no cycle roster: N = cycles present in the null file ({len(present)}), "
              f"which understates N by the cycles with no detection", flush=True)
    return np.sort(np.asarray(sorted(present), dtype=object))


def write_empty(path):
    """Write the artifact that means "no CAAS null": every entry None (p-values stay NA downstream)."""
    with open(path, "wb") as fh:
        pickle.dump({"cycle_levels": None, "null_positions": None, "global": None, "top": None, "bottom": None}, fh)


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

    # The categories of cycle are exactly its distinct values, so this reads the category
    # list instead of scanning every row.
    cycle_levels = resolve_cycle_levels(args.cycle_roster, df["cycle"].cat.categories)

    # Keep only null hits with caas_score > 0 (U_CAAS possibility space)
    df["caas_score"] = pd.to_numeric(df["caas_score"], errors="coerce")
    df = df[df["caas_score"] > 0].copy()
    if df.empty:
        write_empty(args.output)
        print(f"[posenrich_prep] no null hit with caas_score > 0 -> wrote empty {args.output}", flush=True)
        return

    # pos_id is cast to category at once, so the string column of every row collapses to
    # integer codes over the (much smaller) set of distinct positions.
    pos_id = (df["Gene"].astype(str) + ":" + df["Position"].astype(str)).astype("category")
    df["pos_id"] = pos_id
    df["score"] = df["caas_score"]

    long_df = df[["pos_id", "side", "cycle", "score"]]
    null_positions = np.sort(np.asarray(df["pos_id"].unique(), dtype=object))

    out = {"cycle_levels": cycle_levels, "null_positions": null_positions}
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
    print(f"[posenrich_prep] wrote {args.output} ({len(cycle_levels)} cycles, {len(null_positions)} null positions)", flush=True)


if __name__ == "__main__":
    main()
