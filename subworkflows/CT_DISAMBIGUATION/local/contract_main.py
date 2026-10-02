#!/usr/bin/env python3
"""Write the observed contract files (discovery.tab, background.output, background_genes.output, meta_caas/,
caas_convergence_master.csv) from the b_0 directories of the permulation batches. See src/core/contract.py."""
import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.core import contract
from src.core.labelings import design_max_pairs
from src.reporting.disambiguation_writers import _generate_dynamic_fields
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--b0-dirs", nargs="+", required=True, help="Directories written by the batches (b0_observed/)")
    p.add_argument("--design", required=True, help="Observed design (trait file or traitfile_H*.tab directory): the master columns")
    p.add_argument("--discovery-file", default=None,
                   help="A discovery.tab that already exists: write only the meta tables (from it) and the master (from the shards)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()
    configure_logging(verbose=args.verbose)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    dirs = args.b0_dirs  # a sentinel file among them has no batch files
    shards = contract.batch_files(dirs, contract.MASTER_SUFFIX)
    fields = _generate_dynamic_fields(design_max_pairs(Path(args.design)))
    master_path = out / "ct_disambiguation" / "caas_convergence_master.csv"
    master_path.parent.mkdir(parents=True, exist_ok=True)
    if args.discovery_file:
        counts = contract.write_meta(Path(args.discovery_file), out / "meta_caas")
        n_master = contract.write_master(shards, master_path, fields)
        logger.info(f"[contract] {len(shards)} master shards, {n_master} master rows; meta_caas: {counts}")
        return
    disc = contract.batch_files(dirs, contract.DISCOVERY_SUFFIX)
    bg = contract.batch_files(dirs, contract.BACKGROUND_SUFFIX)
    if not (disc or bg or shards):
        logger.info("[contract] the batches carry no b_0 slice: no observed file written")
        master_path.parent.rmdir()
        return
    n_disc = contract.write_discovery(disc, out / "discovery.tab")
    n_bg = contract.write_background(bg, out / "background.output", out / "background_genes.output")
    counts = contract.write_meta(out / "discovery.tab", out / "meta_caas")
    n_master = contract.write_master(shards, master_path, fields)
    logger.info(f"[contract] {len(disc)} discovery files, {n_disc} discovery rows, {n_bg} background genes, "
                f"{len(shards)} master shards, {n_master} master rows; meta_caas: {counts}")


if __name__ == "__main__":
    main()
