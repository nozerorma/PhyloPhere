#!/usr/bin/env python3
"""CLI for the CAAS permulation-excess null (genome-wide *excess* null for FCS).

Loads each gene's precomputed ASR posteriors ONCE and replays N permuted phenotype
labelings (the full-pool `export_perm_discovery` from a bootstrap run + the matching
`resample_*.tab` labelings) over the cached posteriors, scoring each via
analyze_gene_disambiguation / compute_asr_path_score VERBATIM. Emits a long
per-(gene, cycle, position, scheme) asr_path_score table that
scoring_caas_perms.R turns into a genes×N null matrix → FCS p.perm.

See docs/CAAS_PERMULATION_EXCESS.md.
"""

import sys
import csv
import argparse
import logging
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.utils.gene_wrapper import process_all_genes_perms
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


def parse_arguments():
    p = argparse.ArgumentParser(
        description="CAAS permulation-excess null (load-once ASR, replay-N labelings)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--alignment-dir", required=True, help="Directory with alignment files")
    p.add_argument("--tree", required=True, help="Phylogenetic tree (Newick)")
    p.add_argument(
        "--perm-discovery", required=True,
        help="Full-pool export_perm_discovery file (canonical headers, Cycle column)",
    )
    p.add_argument(
        "--resample-dir", required=True,
        help="Directory with resample_*.tab labelings (cycle, fg_csv, bg_csv)",
    )
    p.add_argument("--output-dir", required=True, help="Output directory")
    p.add_argument("--asr-model", default="lg", help="ASR substitution model (default: lg)")
    p.add_argument("--asr-cache-dir", required=True, help="Precomputed ASR cache dir")
    p.add_argument("--posterior-threshold", type=float, default=0.0)
    p.add_argument(
        "--convergence-mode", choices=["focal_clade", "mrca"], default="mrca",
        help="Convergence analysis mode (must match the real run; default: mrca)",
    )
    p.add_argument("--taxid-mapping", default=None)
    p.add_argument("--ensembl-genes-file", default=None)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--max-tasks-per-child", type=int, default=None)
    p.add_argument(
        "--cycles", default=None,
        help="Comma-separated cycle tags to process (default: all cycles in the export)",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args()


def _genes_from_export(perm_discovery_file: Path) -> list:
    """Genes that produced ≥1 CAAS in any cycle (the only genes worth replaying)."""
    genes = set()
    with open(perm_discovery_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            g = (row.get("gene") or "").strip()
            if g:
                genes.add(g)
    return sorted(genes)


def main():
    args = parse_arguments()
    configure_logging(verbose=args.verbose, log_file=args.log_file)

    logger.info("=" * 80)
    logger.info("CAAS Permulation-Excess Null")
    logger.info("=" * 80)
    for k, v in vars(args).items():
        logger.info(f"  {k}: {v}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = _genes_from_export(Path(args.perm_discovery))
    logger.info(f"Genes with CAAS hits across cycles: {len(genes)}")
    if not genes:
        logger.warning("No genes in export_perm_discovery; writing empty null table")

    cycles = None
    if args.cycles:
        cycles = [c.strip() for c in args.cycles.split(",") if c.strip()]

    t0 = time.time()
    out_path = process_all_genes_perms(
        genes=genes,
        alignment_dir=args.alignment_dir,
        tree_file=args.tree,
        perm_discovery_file=args.perm_discovery,
        resample_dir=args.resample_dir,
        taxid_mapping_path=args.taxid_mapping,
        asr_model=args.asr_model,
        asr_cache_dir=args.asr_cache_dir,
        posterior_threshold=args.posterior_threshold,
        convergence_mode=args.convergence_mode,
        workers=args.workers,
        output_dir=output_dir,
        ensembl_genes_file=args.ensembl_genes_file,
        cycles=cycles,
        max_tasks_per_child=args.max_tasks_per_child,
    )
    logger.info(f"Done in {time.time() - t0:.1f}s → {out_path}")


if __name__ == "__main__":
    main()
