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
# File: observed_b0_main.py
#

"""
CAAS_OBSERVED / CAAS_CORE_BATCHED: Scores the observed labeling (b_0) as full records, one master shard per gene.

Reads the b_0 discovery rows a perm-replay batch exported (`<alignment id>.b0.discovery.tsv`, the rows of
discovery.tab), or a discovery.tab that already exists (`--discovery`), scores each gene through `core.observed`
against its ASR (cache first, PAML on a miss) and writes `<gene>.master.csv.gz`: the gene's rows of
caas_convergence_master.csv, in master column order. The merge step (contract_main.py) concatenates the shards.
A missing alignment, unavailable ASR or scoring failure stops the command.

With --strict-output-dir the rows without residue-type overlap are scored again on their own (the analysis of the fully
divergent hypotheses) and written there as the same set of files: the filtered `*.b0.discovery.tsv`, the unchanged
`*.b0.background` and the `<gene>.master.csv.gz` shards of the strict rows. This mode needs --b0-dir
to preserve each gene's discovery and tested-position files.

Called by:  CAAS_CORE_BATCHED (subworkflows/CT/caas_permulation.nf, --b0-dir) and CAAS_OBSERVED
            (subworkflows/CT_DISAMBIGUATION/ct_observed.nf, --discovery) → observed_b0_main.py
Usage:
    python3 observed_b0_main.py --alignment-dir <dir> --tree <newick> (--b0-dir <dir> | --discovery <file>) \\
        --design <trait file or traitfile_H*.tab dir> --output-dir <dir> --asr-cache-dir <dir> [OPTIONS]
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import contextlib
import csv
import json
import hashlib
import logging
import multiprocessing as mp
import shutil
import sys
import time
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.core.driver import load_gene_context
from src.core.labelings import design_max_pairs, observed_pss, read_trait_pairs
from src.core.master import master_fields, write_master_csv
from src.core.contract import BACKGROUND_SUFFIX, UNIVERSE_SUFFIX
from src.core.universe import STRICT_SUFFIX
from src.core.observed import is_strict_row, observed_entries, observed_master_rows, score_observed
from src.data.loaders import load_ensembl_genes
from src.utils.concurrency import init_worker, plan_concurrency, safe_imap_unordered
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)

DISCOVERY_SUFFIX = ".b0.discovery.tsv"


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--alignment-dir", required=True, help="Directory with alignment files")
    p.add_argument("--alignment-manifest", default=None)
    p.add_argument("--manifest-alignment-dir", default=None)
    p.add_argument("--ali-format", default="auto")
    p.add_argument("--tree", required=True, help="Phylogenetic tree (Newick)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--b0-dir", help=f"Directory with the *{DISCOVERY_SUFFIX} files of a perm-replay batch")
    src.add_argument("--discovery", help="A discovery.tab: every gene in it is scored")
    p.add_argument("--design", required=True, help="Observed design: the trait file, or the directory of traitfile_H*.tab")
    p.add_argument("--strict-output-dir", default=None,
                   help="Also write the strict analysis (rows without residue-type overlap) to this directory; needs --b0-dir")
    p.add_argument("--fop-pairs", default=None, help="fop_pairs.tsv holding the b_0 PSS weights (omit for a single contrast)")
    p.add_argument("--output-dir", required=True, help="Directory for the <gene>.master.csv.gz shards")
    p.add_argument("--asr-model", default="lg")
    p.add_argument("--asr-cache-dir", required=True, help="ASR cache dir")
    p.add_argument("--posterior-threshold", type=float, default=0.0)
    p.add_argument("--taxid-mapping", default=None)
    p.add_argument("--ensembl-genes-file", default=None)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--max-tasks-per-child", type=int, default=None)
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args()


# ── Scoring ───────────────────────────────────────────────────────────────────


def read_b0_rows(source, scratch_dir=None, design=None):
    from src.core.discovery import DiscoveryRows
    return DiscoveryRows(source, scratch_dir, design)


def write_strict_inputs(b0_dir, strict_dir):
    """The b_0 files of a batch without the rows that have residue-type overlap: each `*.b0.discovery.tsv` filtered with
    its header kept, each `*.b0.background` copied as it is (the tested positions do not depend on conservation) and each `*.b0.universe_strict` copied as the `*.b0.universe` of the strict directory."""
    strict_dir = Path(strict_dir)
    strict_dir.mkdir(parents=True, exist_ok=True)
    n_kept = n_all = 0
    for path in sorted(Path(b0_dir).glob(f"*{DISCOVERY_SUFFIX}")):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
            fields = reader.fieldnames
        kept = [r for r in rows if is_strict_row(r)]
        n_all += len(rows)
        n_kept += len(kept)
        with open(strict_dir / path.name, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(kept)
    for path in sorted(Path(b0_dir).glob(f"*{BACKGROUND_SUFFIX}")):
        shutil.copyfile(path, strict_dir / path.name)
    # the strict analysis has its own universe (no conserved pair allowed); it takes the name the contract reads
    for path in sorted(Path(b0_dir).glob(f"*{STRICT_SUFFIX}")):
        shutil.copyfile(path, strict_dir / (path.name[:-len(STRICT_SUFFIX)] + UNIVERSE_SUFFIX))
    return n_kept, n_all


_OBSERVED = {}


def _init_observed(config):
    init_worker(1, None)
    _OBSERVED.clear()
    _OBSERVED.update(config)


def _score_gene(job):
    """Pool worker: master rows of one gene (job is the tuple built in main), or (gene, None, None) when it cannot be
    scored. The third element holds the master rows of the strict rows, or None without --strict-output-dir."""
    gene, rows, alignment_path = job
    c = _OBSERVED
    alignment_dir, tree, taxid, model, cache = (c[n] for n in ('alignment_dir','tree','taxid','model','cache'))
    threshold, ensembl, trait_pairs, pss, fields, strict, alignment_format = (c[n] for n in ('threshold','ensembl','trait_pairs','pss','fields','strict','alignment_format'))
    try:
        ctx = load_gene_context(gene, alignment_dir, tree, taxid, model, cache, threshold, ensembl,
                                alignment_path=alignment_path, alignment_format=alignment_format)
        if ctx is None:
            raise RuntimeError(f'ASR unavailable: {gene}')
        results = score_observed(ctx, gene, observed_entries(gene, rows), trait_pairs, pss, threshold)
        strict_rows = None
        if strict:
            strict_in = [r for r in rows if is_strict_row(r)]
            strict_rows = observed_master_rows(
                gene, score_observed(ctx, gene, observed_entries(gene, strict_in), trait_pairs, pss, threshold), fields
            ) if strict_in else []
        return gene, observed_master_rows(gene, results, fields), strict_rows
    except Exception as exc:  # noqa: BLE001 - return failure to the parent, which rejects the batch
        logger.error(f"[observed] {gene} failed: {exc}", exc_info=True)
        return gene, None, None


def main():
    """Score every gene with b_0 rows in parallel and write one <gene>.master.csv.gz per gene."""
    args = parse_arguments()
    configure_logging(verbose=args.verbose, log_file=args.log_file)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from src.core.discovery import position_digest
    from src.utils.io_utils import read_alignment_manifest
    if not 0 <= args.posterior_threshold <= 1:
        raise ValueError('posterior-threshold must lie in [0,1]')
    if args.workers is not None and args.workers < 1:
        raise ValueError('workers must be positive')
    if args.max_tasks_per_child is not None and args.max_tasks_per_child < 1:
        raise ValueError('max-tasks-per-child must be positive')
    if list(out_dir.glob('*.master.csv.gz')) or (out_dir / 'observed_scoring.json').exists():
        raise FileExistsError(f'Existing observed scoring artifacts in {out_dir}')
    trait_pairs = read_trait_pairs(Path(args.design))
    by_gene = read_b0_rows(args.b0_dir or args.discovery, out_dir, trait_pairs)
    alignments = read_alignment_manifest(args.alignment_manifest, args.manifest_alignment_dir or args.alignment_dir) if args.alignment_manifest else {}

    logger.info(f"[observed] {len(by_gene)} genes with b_0 hits")
    fields = master_fields(len(next(iter(trait_pairs.values()))))
    pss = observed_pss(args.fop_pairs, trait_pairs) if args.fop_pairs else None
    ensembl = (load_ensembl_genes(Path(args.ensembl_genes_file)) or set()) if args.ensembl_genes_file else None
    genes = [g for g in by_gene if ensembl is None or g in ensembl]
    excluded = sorted(set(by_gene) - set(genes))
    if args.alignment_manifest and set(genes) - set(alignments):
        raise ValueError(f'Manifest misses discovery genes: {sorted(set(genes)-set(alignments))}')

    strict_dir = Path(args.strict_output_dir) if args.strict_output_dir else None
    if strict_dir is not None:
        if not args.b0_dir:
            raise SystemExit("--strict-output-dir needs --b0-dir to preserve per-gene discovery and background files")
        if list(strict_dir.glob('*.master.csv.gz')) or (strict_dir / 'observed_scoring.json').exists():
            raise FileExistsError(f'Existing strict scoring artifacts in {strict_dir}')
        n_kept, n_all = write_strict_inputs(args.b0_dir, strict_dir)
        logger.info(f"[observed] strict rows (zero residue-type overlap): {n_kept} of {n_all} discovery rows -> {strict_dir}")
    config = dict(alignment_dir=args.alignment_dir, tree=args.tree, taxid=args.taxid_mapping,
                  model=args.asr_model, cache=args.asr_cache_dir, threshold=args.posterior_threshold,
                  ensembl=ensembl, trait_pairs=trait_pairs, pss=pss, fields=fields,
                  strict=strict_dir is not None, alignment_format=args.ali_format)
    jobs = ((g, by_gene[g], str(alignments[g]) if g in alignments else None) for g in genes)
    workers, _ = plan_concurrency(args.workers, 1, logger)
    workers = min(workers, max(1, len(genes)))
    t0 = time.time()
    written, skipped = 0, []
    manager = mp.Pool(processes=workers, maxtasksperchild=args.max_tasks_per_child, initializer=_init_observed, initargs=(config,)) if genes else contextlib.nullcontext(None)
    with manager as pool:
        for gene, rows, strict_rows in (safe_imap_unordered(pool, _score_gene, jobs, chunksize=1) if pool else ()):

            if not rows:
                skipped.append(gene)
                break
            stem = gene
            write_master_csv(rows, out_dir / f"{stem}.master.csv.gz", fields)
            if strict_rows:
                write_master_csv(strict_rows, strict_dir / f"{stem}.master.csv.gz", fields)
            written += 1
    if skipped:
        (out_dir / 'observed_failures.json').write_text(json.dumps({'failed_genes': sorted(skipped)}, indent=2))
        raise RuntimeError(f'{len(skipped)} observed genes failed: {sorted(skipped)[:10]}')

    def ledger(directory, source, strict=False):
        expected = {}
        for gene in genes:
            rows = by_gene[gene]
            if strict:
                rows = [row for row in rows if is_strict_row(row)]
            positions = {int(row['position']) for row in rows}
            if positions:
                expected[gene] = dict(positions=len(positions), sha256=position_digest(positions))
        from src.asr.cache import digest
        src = Path(source)
        sources = [src] if src.is_file() else sorted([*src.glob('*.b0.discovery.tsv'), *src.glob('*.b0.background'), *src.glob('*.b0.universe'), *src.glob('*.b0.universe_strict')])
        from src.core.scores import SCORE_RULE
        record = dict(version=1, strict=strict, score_rule=SCORE_RULE, complete=True, genes=expected, excluded_genes=excluded,
                      tree_sha256=digest(args.tree), model=args.asr_model, posterior_threshold=args.posterior_threshold,
                      design=trait_pairs, taxid_sha256=digest(args.taxid_mapping) if args.taxid_mapping else None,
                      weights_sha256=digest(args.fop_pairs) if args.fop_pairs else None,
                      sources={f.name: digest(f) for f in sources}, weighted=pss is not None,
                      alignments={g:digest(alignments[g]) for g in genes if g in alignments})
        (directory / 'observed_scoring.json').write_text(json.dumps(record, indent=2, sort_keys=True))
    ledger(out_dir, args.b0_dir or args.discovery)
    if strict_dir is not None:
        ledger(strict_dir, strict_dir, strict=True)
    by_gene.close()
    logger.info(f"[observed] {written} genes scored in {time.time() - t0:.1f}s; {len(excluded)} excluded by whitelist")


if __name__ == "__main__":
    main()
