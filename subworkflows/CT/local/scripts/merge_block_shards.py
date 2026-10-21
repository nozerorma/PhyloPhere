#!/usr/bin/env python3
# merge_block_shards.py — Join the shards of one gene that several cycle blocks wrote into a single shard.
# PhyloPhere | subworkflows/CT/local/scripts/

"""
union_shards.sh calls this script for the shard paths that more than one source holds: a gene replayed against
several cycle blocks has one perm_pos_detail/<Gene>.tsv.gz per block.

The merged shard holds the rows of every block in the order of an unblocked replay: a replay writes the cycles of a
gene in the ascending order of their tag as text (`b_100` before `b_11` before `b_3`), the rows of one cycle together
in the order they were produced. Each block's shard already follows that order for its own cycles, so the shards are
merged as sorted streams, never loaded whole. `n_detected` (the number of cycles that detected a position) is a
count over the cycles of one block in each shard; the merged value is the sum of the blocks' counts, which is the
count over all the cycles because the blocks hold disjoint cycles.

Usage: merge_block_shards.py --dest DIR --shared LIST SRC...
  --shared     file with one shard path relative to the sources per line (the paths found in more than one source)
"""
import argparse
import csv
import gzip
import heapq
import os
import sys


def _ordered_rows(path, idx):
    """The rows of a shard keyed by their cycle, checking that they follow the cycle order."""
    last = ""
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if row[idx] < last:
                raise SystemExit(f"[merge_block_shards] {path}: rows are not in cycle order at cycle {row[idx]}")
            last = row[idx]
            yield row[idx], row


def merge_shard(sources, dest):
    """Write the merged shard `dest` from the shards `sources` (paths of the same gene in different blocks)."""
    header = None
    totals = {}
    for path in sources:
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            head = next(reader)
            if header is None:
                header = head
            elif head != header:
                raise SystemExit(f"[merge_block_shards] {path}: header differs from {sources[0]}")
            i_pos, i_nd = head.index("Position"), head.index("n_detected")
            seen = set()
            for row in reader:
                if row[i_pos] not in seen:  # one value per position within a block
                    seen.add(row[i_pos])
                    totals[row[i_pos]] = totals.get(row[i_pos], 0) + int(row[i_nd])
    i_cyc, i_pos, i_nd = header.index("cycle"), header.index("Position"), header.index("n_detected")

    streams = [_ordered_rows(path, i_cyc) for path in sources]
    if os.path.lexists(dest):
        os.remove(dest)  # the destination may be a hard link to one source's shard
    with gzip.open(dest, "wt", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for _, row in heapq.merge(*streams, key=lambda item: item[0]):
            row[i_nd] = str(totals[row[i_pos]])
            writer.writerow(row)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", required=True)
    ap.add_argument("--shared", required=True)
    ap.add_argument("sources", nargs="+")
    args = ap.parse_args(argv)

    with open(args.shared) as handle:
        shared = [line.rstrip("\n") for line in handle if line.strip()]
    for rel in shared:
        sources = [os.path.join(d, rel) for d in args.sources if os.path.isfile(os.path.join(d, rel))]
        merge_shard(sources, os.path.join(args.dest, rel))
    print(f"[merge_block_shards] {len(shared)} shards held by several blocks were joined", file=sys.stderr)


if __name__ == "__main__":
    main()
