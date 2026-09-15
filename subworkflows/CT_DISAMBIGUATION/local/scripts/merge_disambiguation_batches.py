#!/usr/bin/env python3
"""Merge N per-batch ct_disambiguation/ output directories into one.

Each Nextflow batch task (CT_DISAMBIGUATION_RUN_BATCHED) processes a disjoint
gene subset of the same run, so every piece here is either a straight
row-concat or a directory union -- no gene-level aggregation logic, because no
output here is computed across genes (unlike CAAS_PERMS_DISAMBIGUATE's
genome-wide gene_cycle_scores.tsv/perm_pos_sample.tsv/perm_pos_quantiles.tsv,
which is why that stage needs CAAS_PERMS_REBUILD instead of a plain merge).

caas_convergence_master.csv/no_change_debug.csv column schemas are identical
across batches by construction: max_pairs is computed once from the shared,
unpartitioned trait file (see disambiguation_main.py::_compute_max_pairs_from_trait)
and threaded through to every batch's export_from_db call, so every batch
emits the same domain_N_* column set. This script asserts that invariant
rather than silently reconciling a mismatch.
"""

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Merge batched ct_disambiguation/ output directories.")
    parser.add_argument("--batch-dirs", nargs="+", required=True, help="Batch ct_disambiguation/ directories, in order")
    parser.add_argument("--output-dir", required=True, help="Merged ct_disambiguation/ output directory")
    return parser.parse_args()


def _concat_csv(batch_dirs, relpath, out_path):
    header = None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with open(out_path, "w", newline="") as out_f:
        writer = None
        for bdir in batch_dirs:
            src = bdir / relpath
            if not src.exists():
                continue
            with open(src, "r", newline="") as in_f:
                reader = csv.reader(in_f)
                rows = list(reader)
            if not rows:
                continue
            batch_header, batch_rows = rows[0], rows[1:]
            if header is None:
                header = batch_header
                writer = csv.writer(out_f)
                writer.writerow(header)
            elif batch_header != header:
                raise ValueError(
                    f"Schema mismatch merging {relpath}: {src} has header "
                    f"{batch_header}, expected {header} (see max_pairs fix "
                    f"in disambiguation_main.py/gene_wrapper.py)"
                )
            for row in batch_rows:
                writer.writerow(row)
                n_rows += 1
    if header is None:
        # No batch produced this file -- leave nothing behind, matching the
        # unbatched path where an empty result still writes a header-only file.
        return 0
    return n_rows


def _union_dir(batch_dirs, relpath, out_dir):
    n_files = 0
    for bdir in batch_dirs:
        src = bdir / relpath
        if not src.is_dir():
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        for f in sorted(src.iterdir()):
            dest = out_dir / f.name
            if dest.exists():
                raise ValueError(
                    f"Unexpected filename collision merging {relpath}: {f} "
                    f"already present as {dest} -- batches should partition "
                    f"genes disjointly"
                )
            shutil.copy2(f, dest)
            n_files += 1
    return n_files


def _merge_skipped_genes(batch_dirs, out_path):
    genes = []
    seen = set()
    for bdir in batch_dirs:
        src = bdir / "diagnostics" / "skipped_genes.txt"
        if not src.exists():
            continue
        with open(src, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line not in seen:
                    seen.add(line)
                    genes.append(line)
    if not genes:
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"# Skipped {len(genes)} genes (missing alignment or no CAAS positions)\n")
        for g in sorted(genes):
            f.write(f"{g}\n")
    return len(genes)


def _merge_summary_json(batch_dirs, out_path):
    by_gene_counts = {}
    num_pairs = None
    schema_version = None
    for bdir in batch_dirs:
        src = bdir / "caas_convergence_summary.json"
        if not src.exists():
            continue
        with open(src, "r", encoding="utf-8") as f:
            obj = json.load(f)
        meta = obj.get("metadata", {})
        batch_num_pairs = meta.get("num_pairs")
        if num_pairs is None:
            num_pairs = batch_num_pairs
            schema_version = meta.get("schema_version")
        elif batch_num_pairs != num_pairs:
            raise ValueError(
                f"Schema mismatch merging caas_convergence_summary.json: {src} "
                f"has num_pairs={batch_num_pairs}, expected {num_pairs}"
            )
        for gene, count in obj.get("by_gene_counts", {}).items():
            if gene in by_gene_counts:
                raise ValueError(
                    f"Unexpected gene overlap merging caas_convergence_summary.json: "
                    f"{gene} counted in more than one batch"
                )
            by_gene_counts[gene] = count

    if not by_gene_counts:
        return

    summary_obj = {
        "metadata": {
            "num_genes": len(by_gene_counts),
            "total_positions": sum(by_gene_counts.values()),
            "num_pairs": num_pairs,
            "schema_version": schema_version,
        },
        "by_gene_counts": by_gene_counts,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(summary_obj, f, indent=2, ensure_ascii=False)


def main():
    args = parse_args()
    batch_dirs = [Path(p) for p in args.batch_dirs]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_master = _concat_csv(batch_dirs, "caas_convergence_master.csv", out_dir / "caas_convergence_master.csv")
    print(f"[merge_disambiguation_batches] caas_convergence_master.csv: {n_master} rows")

    n_no_change = _concat_csv(
        batch_dirs, "diagnostics/no_change_debug.csv", out_dir / "diagnostics" / "no_change_debug.csv"
    )
    print(f"[merge_disambiguation_batches] diagnostics/no_change_debug.csv: {n_no_change} rows")

    n_json = _union_dir(batch_dirs, "json_summaries", out_dir / "json_summaries")
    print(f"[merge_disambiguation_batches] json_summaries/: {n_json} files")

    n_nodes = _union_dir(batch_dirs, "diagnostics/node_dumps", out_dir / "diagnostics" / "node_dumps")
    print(f"[merge_disambiguation_batches] diagnostics/node_dumps/: {n_nodes} files")

    n_skipped = _merge_skipped_genes(batch_dirs, out_dir / "diagnostics" / "skipped_genes.txt")
    print(f"[merge_disambiguation_batches] diagnostics/skipped_genes.txt: {n_skipped} genes")

    _merge_summary_json(batch_dirs, out_dir / "caas_convergence_summary.json")
    print("[merge_disambiguation_batches] caas_convergence_summary.json merged")

    # aggregation.sqlite3 is per-batch only; nothing downstream reads the DB
    # directly (only the exported files above), so keep each batch's DB
    # alongside the merged outputs rather than merging it.
    for i, bdir in enumerate(batch_dirs, start=1):
        src = bdir / "aggregation.sqlite3"
        if src.exists():
            shutil.copy2(src, out_dir / f"aggregation_batch_{i:05d}.sqlite3")


if __name__ == "__main__":
    sys.exit(main())
