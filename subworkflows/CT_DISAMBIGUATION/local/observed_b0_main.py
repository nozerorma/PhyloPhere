#!/usr/bin/env python3
"""The observed labeling (b_0) as full records: one master shard per gene.

Reads the b_0 discovery rows a perm-replay batch exported (`<alignment id>.b0.discovery.tsv`, the rows of
discovery.tab), or a discovery.tab that already exists (`--discovery`), scores each gene through `core.observed` against its ASR (cache first, PAML on a miss) and writes
`<gene>.master.csv.gz`: the gene's rows of caas_convergence_master.csv, in master column order. The merge step
concatenates the shards. A gene with no alignment or ASR is left out with a warning.
"""
import argparse
import csv
import logging
import multiprocessing as mp
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.core.driver import load_gene_context
from src.core.labelings import design_max_pairs, observed_pss, read_trait_pairs
from src.core.master import master_fields, write_master_csv
from src.core.observed import observed_entries, observed_master_rows, score_observed
from src.data.loaders import load_ensembl_genes
from src.utils.concurrency import init_worker, plan_concurrency
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)

DISCOVERY_SUFFIX = ".b0.discovery.tsv"


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--alignment-dir", required=True, help="Directory with alignment files")
    p.add_argument("--tree", required=True, help="Phylogenetic tree (Newick)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--b0-dir", help=f"Directory with the *{DISCOVERY_SUFFIX} files of a perm-replay batch")
    src.add_argument("--discovery", help="A discovery.tab: every gene in it is scored")
    p.add_argument("--design", required=True, help="Observed design: the trait file, or the directory of traitfile_H*.tab")
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


def read_b0_rows(source):
    """{gene: [discovery rows]} of every b_0 discovery file in a directory, or of one discovery.tab; rows in file order."""
    by_gene = {}
    source = Path(source)
    paths = [source] if source.is_file() else sorted(source.glob(f"*{DISCOVERY_SUFFIX}"))
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                by_gene.setdefault(row["gene"], []).append(row)
    return by_gene


def _score_gene(job):
    (gene, rows, alignment_dir, tree, taxid, model, cache, threshold, ensembl, trait_pairs, pss, fields) = job
    try:
        ctx = load_gene_context(gene, alignment_dir, tree, taxid, model, cache, threshold, ensembl)
        if ctx is None:
            return gene, None
        results = score_observed(ctx, gene, observed_entries(gene, rows), trait_pairs, pss, threshold)
        return gene, observed_master_rows(gene, results, fields)
    except Exception as exc:  # noqa: BLE001 - one gene must not end the batch; the summary reports it
        logger.error(f"[observed] {gene} failed: {exc}", exc_info=True)
        return gene, None


def main():
    args = parse_arguments()
    configure_logging(verbose=args.verbose, log_file=args.log_file)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    by_gene = read_b0_rows(args.b0_dir or args.discovery)
    logger.info(f"[observed] {len(by_gene)} genes with b_0 hits")
    if not by_gene:
        return

    fields = master_fields(design_max_pairs(Path(args.design)))
    trait_pairs = read_trait_pairs(Path(args.design))
    pss = observed_pss(args.fop_pairs) if args.fop_pairs else None
    ensembl = (load_ensembl_genes(Path(args.ensembl_genes_file)) or set()) if args.ensembl_genes_file else None
    if ensembl is not None:
        by_gene = {g: r for g, r in by_gene.items() if g in ensembl}

    jobs = [(g, by_gene[g], args.alignment_dir, args.tree, args.taxid_mapping, args.asr_model, args.asr_cache_dir,
             args.posterior_threshold, ensembl, trait_pairs, pss, fields) for g in sorted(by_gene)]
    workers, _ = plan_concurrency(args.workers, 1, logger)
    t0 = time.time()
    written, skipped = 0, []
    with mp.Pool(processes=workers, maxtasksperchild=args.max_tasks_per_child, initializer=init_worker, initargs=(1, None)) as pool:
        for gene, rows in pool.imap_unordered(_score_gene, jobs, chunksize=1):
            if not rows:
                skipped.append(gene)
                continue
            stem = gene.replace("/", "__").replace("\\", "__").strip() or "_"
            write_master_csv(rows, out_dir / f"{stem}.master.csv.gz", fields)
            written += 1
    logger.info(f"[observed] {written} genes scored in {time.time() - t0:.1f}s; {len(skipped)} left out (no alignment, ASR or results)"
                + (f": {sorted(skipped)[:10]}" if skipped else ""))


if __name__ == "__main__":
    main()
