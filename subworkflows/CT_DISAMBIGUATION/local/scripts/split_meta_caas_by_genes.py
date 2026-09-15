#!/usr/bin/env python3
"""Split a CT_META_CAAS metadata table into per-batch files, chunked by gene.

Mirrors disambiguation_main.py's own gene-extraction rule exactly
(``GenePos.split("_")[0]``) so a batch's file contains precisely the rows
disambiguation_main.py would itself attribute to the genes in that batch --
running it unbatched on one such file must be byte-identical to running it on
the corresponding gene slice of the full run.
"""

import argparse
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a CAAS metadata table into per-batch files by gene."
    )
    parser.add_argument("--meta-caas", required=True, help="Path to global_meta_caas.tsv (or meta_caas.tsv)")
    parser.add_argument("--batch-size", type=int, required=True, help="Number of genes per batch")
    parser.add_argument("--outdir", required=True, help="Directory to write batch_*.meta_caas.tsv files into")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.meta_caas, "r") as fh:
        header = fh.readline()
    sep = "\t" if "\t" in header else ","

    df = pd.read_csv(args.meta_caas, sep=sep)
    if "GenePos" not in df.columns:
        raise ValueError(f"Expected a GenePos column in {args.meta_caas}, found: {list(df.columns)}")

    gene_of_row = df["GenePos"].apply(
        lambda x: x.split("_")[0] if "_" in str(x) else None
    )
    unique_genes = sorted(g for g in gene_of_row.dropna().unique().tolist())

    if not unique_genes:
        raise ValueError(f"No genes could be parsed from GenePos in {args.meta_caas}")

    batch_size = max(1, args.batch_size)
    n_batches = 0
    for start in range(0, len(unique_genes), batch_size):
        n_batches += 1
        batch_genes = set(unique_genes[start:start + batch_size])
        batch_df = df[gene_of_row.isin(batch_genes)]
        out_path = f"{args.outdir}/batch_{n_batches:05d}.meta_caas.tsv"
        batch_df.to_csv(out_path, sep=sep, index=False)
        print(f"[split_meta_caas_by_genes] {out_path}: {len(batch_genes)} genes, {len(batch_df)} rows")

    print(f"[split_meta_caas_by_genes] Wrote {n_batches} batches ({len(unique_genes)} genes total) to {args.outdir}")


if __name__ == "__main__":
    sys.exit(main())
