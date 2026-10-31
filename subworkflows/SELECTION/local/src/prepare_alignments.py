#!/usr/bin/env python3
"""Prepare one gene or a manifest with identical validation and exclusion rules.

Biological exclusions (no informative on-tree sequence) are audited; malformed
inputs and worker failures make the entire task fail. A batch parses its tree once.
"""
import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path
from Bio import SeqIO
from alignment_io import informative, read_alignment, tree_taxa


def prepare(gene, alignment, taxa):
    output = Path(f"{gene}.filtered.fa")
    audit = Path(f"{gene}.prep.json")
    output.unlink(missing_ok=True)
    audit.unlink(missing_ok=True)
    result = {"gene": gene, "alignment": Path(alignment).name}
    try:
        records = read_alignment(alignment)
        no_information = [r.id for r in records if not informative(r)]
        off_tree = [r.id for r in records if r.id not in taxa]
        kept = [r for r in records if informative(r) and r.id in taxa]
        result.update(input_taxa=len(records), retained_taxa=[r.id for r in kept],
                      excluded_uninformative=no_information, excluded_off_tree=off_tree,
                      alignment_columns=len(records[0].seq),
                      status="prepared" if kept else "excluded_no_usable_tree_taxa")
        if kept:
            SeqIO.write(kept, output, "fasta")
    except Exception as exc:
        output.unlink(missing_ok=True)
        result.update(status="technical_failure", error=str(exc))
    audit.write_text(json.dumps(result, indent=2) + "\n")
    return result


def validate_gene(gene):
    if (not gene or gene in {".", ".."} or gene != gene.strip()
            or any(c in gene for c in "/\\\t\r\n")):
        raise ValueError(f"Invalid gene identifier: {gene!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--workers", type=int, default=1)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--manifest")
    group.add_argument("--alignment")
    parser.add_argument("--gene")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    taxa = tree_taxa(args.tree)
    if args.manifest:
        jobs, seen = [], set()
        with open(args.manifest, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                if len(row) != 2:
                    raise ValueError(f"Invalid manifest row: {row!r}")
                gene, name = row
                validate_gene(gene)
                if gene in seen or not name or Path(name).name != name:
                    raise ValueError(f"Duplicate gene or invalid filename: {row!r}")
                seen.add(gene)
                jobs.append((gene, Path("alignments") / name))
        if not jobs:
            raise ValueError("Empty batch manifest")
    else:
        validate_gene(args.gene)
        jobs = [(args.gene, args.alignment)]
    # Parsing is CPU-bound Python work. Processes use the allocated CPUs rather
    # than competing under the GIL; a singleton needs no extra worker process.
    workers = min(args.workers, len(jobs))
    if workers == 1:
        results = [prepare(*job, taxa) for job in jobs]
    else:
        genes, alignments = zip(*jobs)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(prepare, genes, alignments, repeat(taxa)))
    failures = [r for r in results if r["status"] == "technical_failure"]
    if failures:
        raise RuntimeError("Alignment preparation failed: " + "; ".join(
            f"{r['gene']}: {r['error']}" for r in failures))


if __name__ == "__main__":
    main()
