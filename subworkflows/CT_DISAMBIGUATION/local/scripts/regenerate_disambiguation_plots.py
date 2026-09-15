#!/usr/bin/env python3
"""Regenerate ct_disambiguation/plots/ from the merged master CSV.

disambiguation_main.py calls generate_bulk_plots() once per batch when run
unbatched; when CT_DISAMBIGUATION_RUN is chunked by gene, running that same
per-batch call would only ever see one batch's gene subset, so it is skipped
in CT_DISAMBIGUATION_RUN_BATCHED and re-run here, once, against the merged
master CSV -- same call, same output shape, just after the merge instead of
before it.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from plots.plotter import generate_bulk_plots  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description="Regenerate bulk plots from a merged ct_disambiguation master CSV.")
    parser.add_argument("--caas-csv", required=True, help="Merged caas_convergence_master.csv")
    parser.add_argument("--output-dir", required=True, help="Directory to write plots/ into")
    parser.add_argument("--asr-cache-dir", default=None, help="Shared ASR cache directory")
    parser.add_argument("--node-dumps-root", default=None, help="Merged diagnostics/node_dumps directory")
    parser.add_argument("--ensembl-genes-file", default=None, help="TSV/CSV with gene column")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    caas_csv = Path(args.caas_csv)

    if not caas_csv.exists():
        print(f"[regenerate_disambiguation_plots] Skipping; master CSV not found: {caas_csv}")
        return 0

    asr_root = Path(args.asr_cache_dir) if args.asr_cache_dir else output_dir / "asr"
    node_dumps_root = Path(args.node_dumps_root) if args.node_dumps_root else output_dir / "diagnostics" / "node_dumps"
    ensembl_path = Path(args.ensembl_genes_file) if args.ensembl_genes_file else None

    generate_bulk_plots(
        caas_csv=caas_csv,
        output_dir=output_dir / "plots",
        ensembl_csv=ensembl_path,
        asr_root=asr_root,
        node_dumps_root=node_dumps_root,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
