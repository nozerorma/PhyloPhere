#!/usr/bin/env python3
# gene_lengths.py — One table of gene lengths for every gene of the alignment manifest.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Writes gene_lengths.tsv, the single source of gene lengths of a run, read by the observed gene filter
(filter_caas_genes.py) and by every block merge of the permulation null (merge_null_block.py), so the
gene-removal rule evaluates the same genes with the same lengths on both sides.

Called by:  GENE_LENGTHS process (ctpp_clustfilter.nf)
Inputs:     --genes       manifest of gene IDs, one per line (the alignments of the run); without it, every gene of
                          the annotation and of the MAP directory
            --annotation  gene annotation TSV (gene, length ...), the fallback source
            --map-dir     the trimmer's MAP tables: untrimmed alignment columns, the preferred source
Outputs:    -o            gene, length, length_source for every gene; length is NA and length_source is "missing"
                          for a gene with no length in any source
            --summary     JSON with the number of genes, genes with a length, genes without one and counts per source
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION" / "local"))
from src.core.gene_rules import LENGTH_TABLE_COLUMNS, load_gene_lengths  # noqa: E402


def build_table(genes, annotation, map_dir, map_suffix=".map.tsv"):
    """Rows (gene, length, length_source) of every gene, sorted by gene; genes=None takes the genes of the sources."""
    sources = {}
    lengths = load_gene_lengths(path=annotation, map_dir=map_dir, map_suffix=map_suffix,
                                genes=genes, sources=sources)
    names = sorted(set(genes) if genes is not None else lengths)
    return [(g, lengths.get(g), sources.get(g, "missing")) for g in names]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=Path, default=None, help="manifest of gene IDs, one per line")
    ap.add_argument("--annotation", default=None, help="gene annotation TSV (gene, length)")
    ap.add_argument("--map-dir", default=None, help="directory of the trimmer's MAP tables")
    ap.add_argument("--map-suffix", default=".map.tsv")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()
    if not args.annotation and not args.map_dir:
        ap.error("a length source is needed: --annotation or --map-dir")
    genes = None
    if args.genes is not None:
        genes = [g.strip() for g in args.genes.read_text().splitlines() if g.strip()]
        if len(set(genes)) != len(genes):
            ap.error("duplicate gene IDs in the manifest")
    rows = build_table(genes, args.annotation, args.map_dir, args.map_suffix)
    with open(args.output, "w") as out:
        out.write("\t".join(LENGTH_TABLE_COLUMNS) + "\n")
        for gene, length, source in rows:
            out.write(f"{gene}\t{'NA' if length is None else repr(float(length))}\t{source}\n")
    n_missing = sum(1 for _, length, _ in rows if length is None)
    counts = {}
    for _, length, source in rows:
        counts[source] = counts.get(source, 0) + 1
    summary = {"n_genes": len(rows), "n_with_length": len(rows) - n_missing, "n_missing": n_missing, "by_source": counts}
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=1) + "\n")
    print(f"gene lengths: {summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
