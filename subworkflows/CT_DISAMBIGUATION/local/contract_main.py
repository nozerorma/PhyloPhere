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
    p.add_argument("--output-dir", required=True)
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()
    configure_logging(verbose=args.verbose)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    dirs = args.b0_dirs  # a sentinel file among them has no batch files
    disc = contract.batch_files(dirs, contract.DISCOVERY_SUFFIX)
    n_disc = contract.write_discovery(disc, out / "discovery.tab")
    n_bg = contract.write_background(contract.batch_files(dirs, contract.BACKGROUND_SUFFIX), out / "background.output",
                                     out / "background_genes.output")
    counts = contract.write_meta(out / "discovery.tab", out / "meta_caas")
    fields = _generate_dynamic_fields(design_max_pairs(Path(args.design)))
    shards = contract.batch_files(dirs, contract.MASTER_SUFFIX)
    n_master = contract.write_master(shards, out / "caas_convergence_master.csv", fields)
    logger.info(f"[contract] {len(disc)} discovery files, {n_disc} discovery rows, {n_bg} background genes, "
                f"{len(shards)} master shards, {n_master} master rows; meta_caas: {counts}")


if __name__ == "__main__":
    main()
