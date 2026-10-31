#!/usr/bin/env python3
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: contract_main.py
#

"""
CAAS_OBSERVED / CAAS_CORE_OBSERVED: Writes the observed contract files from the b_0 directories of the permulation batches.

The files are discovery.tab, background.output, background_genes.output, meta_caas/ and
ct_disambiguation/caas_convergence_master.csv (see src/core/contract.py). With --discovery-file (a
discovery.tab that already exists) only meta_caas/ (from that file) and the master (from the
shards) are written.

Called by:  CAAS_CORE_OBSERVED (subworkflows/CT/caas_permulation.nf) and CAAS_OBSERVED
            (subworkflows/CT_DISAMBIGUATION/ct_observed.nf) → contract_main.py
Usage:
    python3 contract_main.py --b0-dirs <dir> [<dir> ...] --design <trait file or traitfile_H*.tab dir> \\
        --output-dir <dir> [--discovery-file <discovery.tab>] [--verbose]
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import json
import tempfile
import shutil
import sys
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.core import contract
from src.core.labelings import design_max_pairs
from src.core.master import master_fields
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    """Write the observed tables of the b_0 shards (whole set, or only meta and master with --discovery-file)."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--b0-dirs", nargs="+", required=True, help="Directories written by the batches (b0_observed/)")
    p.add_argument("--design", required=True, help="Observed design (trait file or traitfile_H*.tab directory): the master columns")
    p.add_argument("--discovery-file", default=None,
                   help="A discovery.tab that already exists: write only the meta tables (from it) and the master (from the shards)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()
    configure_logging(verbose=args.verbose)

    from src.asr.cache import digest
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    dirs = [Path(d) for d in args.b0_dirs]
    if any(not d.is_dir() for d in dirs):
        raise ValueError('Every b0 input must be an existing directory')
    shards = contract.batch_files(dirs, contract.MASTER_SUFFIX)
    disc = contract.batch_files(dirs, contract.DISCOVERY_SUFFIX)
    bg = contract.batch_files(dirs, contract.BACKGROUND_SUFFIX)
    universe = contract.batch_files(dirs, contract.UNIVERSE_SUFFIX)
    if not args.discovery_file and not (disc or bg or shards):
        return
    expected, excluded, ledgers = {}, set(), []
    for directory in dirs:
        stamp = directory / 'observed_scoring.json'
        if not stamp.is_file():
            if list(directory.glob('*' + contract.DISCOVERY_SUFFIX)) or list(directory.glob('*' + contract.MASTER_SUFFIX)) or list(directory.glob('*' + contract.BACKGROUND_SUFFIX)) or list(directory.glob('*' + contract.UNIVERSE_SUFFIX)) or args.discovery_file:
                raise ValueError(f'Missing observed completion record: {directory}')
            continue
        record = json.loads(stamp.read_text())
        from src.core.scores import SCORE_RULE
        if record.get('version') != 1 or record.get('complete') is not True or record.get('score_rule') != SCORE_RULE:
            raise ValueError(f'Invalid observed completion record: {stamp}')
        if set(expected) & set(record['genes']):
            raise ValueError('Genes overlap between observed batches')
        from src.core.labelings import read_trait_pairs
        if record.get('design') != json.loads(json.dumps(read_trait_pairs(Path(args.design)))):
            raise ValueError('Observed design changed after scoring')
        if ledgers and any(record.get(k) != ledgers[0].get(k) for k in
                           ('model','posterior_threshold','tree_sha256','taxid_sha256','weights_sha256','weighted')):
            raise ValueError('Observed batches use inconsistent scientific settings')
        expected.update(record['genes'])
        excluded.update(record['excluded_genes'])
        for name, checksum in record['sources'].items():
            source = Path(args.discovery_file) if args.discovery_file else directory / name
            if source.name != name or digest(source) != checksum:
                raise ValueError(f'Observed discovery changed after scoring: {source}')
        ledgers.append(record)
    if expected.keys() & excluded:
        raise ValueError('A gene was excluded in one batch but scored in another')
    fields = master_fields(design_max_pairs(Path(args.design)))
    targets = ['ct_disambiguation', 'meta_caas'] + ([] if args.discovery_file else ['discovery.tab','background.output','background_genes.output'])
    if any((out / name).exists() for name in targets):
        raise FileExistsError(f'Existing observed contract artifacts in {out}')
    with tempfile.TemporaryDirectory(prefix='.contract-', dir=out) as scratch:
        stage = Path(scratch)
        (stage / 'ct_disambiguation').mkdir()
        audit = {}
        n_master = contract.write_master(shards, stage / 'ct_disambiguation/caas_convergence_master.csv', fields, audit=audit)
        if audit != expected:
            raise ValueError(f'Observed master coverage differs from scored discovery: missing={sorted(set(expected)-set(audit))[:10]}, extra={sorted(set(audit)-set(expected))[:10]}')
        if args.discovery_file:
            discovery = Path(args.discovery_file)
        else:
            discovery = stage / 'discovery.tab'
            contract.write_discovery(disc, discovery)
            if bg and not universe:
                raise ValueError('Tested positions without a positional universe (*.b0.universe): the batches predate it, regenerate them')
            contract.write_background(bg, universe, stage / 'background.output', stage / 'background_genes.output')
        n_meta = contract.write_meta(discovery, stage / 'meta_caas')
        (stage / 'ct_disambiguation/observed_scoring.json').write_text(json.dumps(
            dict(version=1, complete=True, genes=expected, excluded_genes=sorted(excluded), batches=ledgers), indent=2, sort_keys=True))
        for name in targets:
            shutil.move(str(stage / name), str(out / name))
    logger.info(f'[contract] {len(shards)} master shards, {n_master} master rows, {n_meta} meta rows')


if __name__ == "__main__":
    main()
