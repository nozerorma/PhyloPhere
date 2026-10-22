#!/usr/bin/env python3
# reduce_null_blocks.py — The small per-block outputs of the permulation null joined into the tables of the whole null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/

"""
CAAS_CORE_MERGE_REDUCE: joins what merge_null_block.py wrote for each block of cycles.

Every file is small (a row per cycle or per removed unit), so the join reads them whole. The matrices are joined by
scoring_caas_perms.R and the scored positions stay in their blocks (the position-level null is consumed block by block).

Called by:  CAAS_CORE_MERGE_REDUCE (subworkflows/CT/caas_permulation.nf) → reduce_null_blocks.py
Inputs:     --blocks  the output directories of merge_null_block.py, one per block of permuted cycles
Outputs:    <out>/cycle_roster.txt   every replayed cycle, the ones with no detection included: N of every consumer
            <out>/null_positions.parquet  per scored position of the null (gene, position): nd, the number of cycles that
                                     score it, and first_rank, the position in replay order (0-based) of the first one
            <out>/removed_units.tsv  (cycle, Gene) units dropped by gene removal, when any block applied it
            <out>/perm_pos_sample.tsv, perm_pos_quantiles.tsv, cycle_summary.tsv, cycle_hist.tsv
            <out>/counts.json        units, removed units, scored positions and the cycles with a row
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from src.core.scores import cycle_order_key  # noqa: E402


def read_rows(path: Path):
    """(header, rows) of a TSV."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


def join_tables(blocks, name, out: Path, sort_rows: bool = False):
    """Concatenate the `name` table of every block that has it (same header) into out/name; False when none has it."""
    header, rows = None, []
    for b in blocks:
        path = b / name
        if not path.is_file():
            continue
        h, r = read_rows(path)
        if header is not None and h != header:
            raise ValueError(f"{path} has columns {h}, expected {header}")
        header = h
        rows.extend(r)
    if header is None:
        return False
    if sort_rows:
        rows.sort()
    with open(out / name, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return True


def join_detections(blocks, roster, out: Path):
    """Sum nd and take the earliest first cycle, per position, over the blocks; written as null_positions.parquet.

    The blocks are folded in one at a time, so memory follows the number of distinct scored positions of the null and
    not the number of blocks.
    """
    ordered = sorted(roster, key=cycle_order_key)
    acc = None
    for b in blocks:
        t = pq.read_table(b / "detections.parquet")
        rank = pc.index_in(t["first"], value_set=pa.array(ordered, pa.string())).cast(pa.int32())
        t = pa.table({"gene": t["gene"], "position": t["position"], "nd": t["nd"], "first_rank": rank})
        acc = t if acc is None else pa.concat_tables([acc, t])
        acc = acc.group_by(["gene", "position"]).aggregate([("nd", "sum"), ("first_rank", "min")])
        acc = pa.table({"gene": acc["gene"], "position": acc["position"], "nd": acc["nd_sum"].cast(pa.int32()),
                        "first_rank": acc["first_rank_min"]})
    if acc is None:
        acc = pa.table({"gene": pa.array([], pa.string()), "position": pa.array([], pa.int32()),
                        "nd": pa.array([], pa.int32()), "first_rank": pa.array([], pa.int32())})
    acc = acc.sort_by([("gene", "ascending"), ("position", "ascending")])
    pq.write_table(acc, out / "null_positions.parquet", compression="zstd")
    return acc.num_rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blocks", required=True, nargs="+", type=Path, help="merge_null_block.py output directories")
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    blocks = sorted(args.blocks, key=lambda p: p.name)

    roster = set()
    for b in blocks:
        roster.update(c for c in (b / "cycles.txt").read_text().splitlines() if c)
    (args.output_dir / "cycle_roster.txt").write_text("".join(f"{c}\n" for c in sorted(roster)))

    join_tables(blocks, "removed_units.tsv", args.output_dir, sort_rows=True)
    for name in ("perm_pos_sample.tsv", "perm_pos_quantiles.tsv", "cycle_summary.tsv", "cycle_hist.tsv"):
        if not join_tables(blocks, name, args.output_dir):
            print(f"no block wrote {name}", file=sys.stderr)
            return 1

    n_positions = join_detections(blocks, roster, args.output_dir)

    total = {"units": 0, "removed_units": 0, "scored_positions": 0, "cycles": 0, "cycles_present": 0}
    rules = set()
    for b in blocks:
        c = json.loads((b / "counts.json").read_text())
        rules.add(c.get("score_rule"))
        for k in ("units", "removed_units", "scored_positions", "cycles"):
            total[k] += c[k]
        total["cycles_present"] += len(c["cycles_present"])
    if len(rules) != 1:
        print(f"the blocks were scored with different rules: {sorted(map(str, rules))}", file=sys.stderr)
        return 1
    total["score_rule"] = rules.pop()      # the rule of the position score of the null: SCORING refuses a null built with another one
    (args.output_dir / "counts.json").write_text(json.dumps(total))
    silent = total["cycles"] - total["cycles_present"]
    if silent:
        print(f"[reduce] {silent} of {total['cycles']} replayed cycles left no row: they count in N as cycles with no signal",
              file=sys.stderr)
    print(f"[reduce] {len(blocks)} blocks, {total['cycles']} cycles, {total['units']} units "
          f"({total['removed_units']} removed), {total['scored_positions']} scored positions, "
          f"{n_positions} distinct scored positions", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
