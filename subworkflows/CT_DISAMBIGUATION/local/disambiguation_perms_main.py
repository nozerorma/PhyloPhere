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
# File: disambiguation_perms_main.py
#

"""
CAAS_CORE_BATCHED: Replays permuted phenotype labelings over cached ASR posteriors (the CAAS permulation null).

Loads each gene's precomputed ASR posteriors once and replays N permuted phenotype labelings (the
full-pool perm-discovery export of a perm-replay run, with the matching resample_*.tab labelings)
over them, scoring each through analyze_gene_disambiguation as the observed pipeline does
(src/utils/gene_wrapper.py, process_all_genes_perms). The real labeling (b_0) is
replayed on its own into <output-dir>/b0, so it never enters the null.

Pass A writes one per-(gene, cycle, position, scheme) detail shard per gene to
<output-dir>/perm_pos_detail/. With --detail-only the run stops there, because pass B (scoring
against genome-wide pools) runs once over the union of the batches (reaggregate_perm_scores.py).
Without it the same process also writes gene_cycle_scores.tsv, perm_pos_cycle_caas.tsv.gz,
perm_pos_quantiles.tsv and perm_pos_sample.tsv; scoring_caas_perms.R turns gene_cycle_scores.tsv
into the genes x N null matrix behind the FCS p.perm.

Called by:  CAAS_CORE_BATCHED (subworkflows/CT/caas_permulation.nf) → disambiguation_perms_main.py --detail-only
Usage:
    python3 disambiguation_perms_main.py --alignment-dir <dir> --tree <newick> --perm-discovery <file or dir> \\
        --resample-dir <dir> --output-dir <dir> --asr-cache-dir <dir> [--detail-only] [OPTIONS]
"""

# ── Standard library ──────────────────────────────────────────────────────────
import sys
import csv
import argparse
import logging
import time
from pathlib import Path
from typing import Tuple

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.utils.gene_wrapper import process_all_genes_perms, _read_resample_labelings
from src.convergence.fop_pool import base_cycle
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_arguments():
    p = argparse.ArgumentParser(
        description="CAAS permulation-excess null (load-once ASR, replay-N labelings)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--alignment-dir", required=True, help="Directory with alignment files")
    p.add_argument("--tree", required=True, help="Phylogenetic tree (Newick)")
    p.add_argument(
        "--perm-discovery", required=True,
        help="Full-pool export_perm_discovery file or directory of files (canonical headers, Cycle column)",
    )
    p.add_argument(
        "--resample-dir", required=True,
        help="Directory with resample_*.tab labelings (cycle, fg_csv, bg_csv)",
    )
    p.add_argument("--output-dir", required=True, help="Output directory")
    p.add_argument("--asr-model", default="lg", help="ASR substitution model (default: lg)")
    p.add_argument("--asr-cache-dir", required=True, help="Precomputed ASR cache dir")
    p.add_argument("--posterior-threshold", type=float, default=0.0)
    p.add_argument("--taxid-mapping", default=None)
    p.add_argument("--ensembl-genes-file", default=None)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--max-tasks-per-child", type=int, default=None)
    p.add_argument(
        "--cycles", default=None,
        help="Comma-separated cycle tags to process (default: all cycles in the export)",
    )
    p.add_argument(
        "--fop-pairs", default=None,
        help="resample_fop_pairs.tsv (FOP mirror): per-(cycle, hypothesis, domain) PSS "
             "weights. When given, the resample dir is expected to hold resample_fop.tab "
             "('<base>~H<m>' labelings) and each base cycle's hypotheses are domain-pooled.",
    )
    # ── CT_POSTPROC filters applied to the null candidate pool ────────────────
    p.add_argument("--postproc-filter", action="store_true",
                   help="Apply the observed CT_POSTPROC cluster + gene filters to "
                        "the per-cycle null CAAS pool before scoring (distribution-"
                        "matches caas_perms.rds to the observed filtered_discovery).")
    p.add_argument("--gene-lengths", default=None,
                   help="Gene annotation TSV (gene, length ...) for the extreme-gene "
                        "density filter. Same file as CT_POSTPROC's gene_ensembl_file.")
    p.add_argument("--clust-minlen", type=int, default=3,
                   help="CT_FILTER minlen (params.filter_minlen; default 3)")
    p.add_argument("--clust-maxcaas", type=float, default=0.7,
                   help="CT_FILTER maxcaas density (params.filter_maxcaas; default 0.7)")
    p.add_argument("--keep-clusters", action="store_true",
                   help="Keep cluster-train positions in the scored pool (params.remove_caas_clusters "
                        "false); they still count for the dubious-gene test.")
    p.add_argument("--train-map-dir", default=None,
                   help="Directory of the trimmer's per-gene MAP tables. Cluster trains then measure their span "
                        "in untrimmed alignment columns (a gene without a MAP keeps trimmed coordinates). "
                        "Needs --postproc-filter; the observed filter must use the same directory.")
    p.add_argument("--train-map-suffix", default=".map.tsv",
                   help="File-name tail of the MAP tables (default .map.tsv)")
    p.add_argument("--gene-filter-mode", default="none",
                   choices=["none", "extreme", "dubious", "both"],
                   help="CAAS_FILTER_GENES mode (params.gene_filter_mode)")
    p.add_argument("--iqr-multiplier", type=float, default=3.0)
    p.add_argument("--extreme-percentile", type=float, default=0.99)
    p.add_argument("--detail-only", action="store_true",
                   help="Stop after pass A: write only perm_pos_detail/ (pass B scores against genome-wide pools, "
                        "so it runs once over the union of the shards).")
    p.add_argument("--seed", type=int, default=1998,
                   help="Pipeline seed (params.seed): perm_pos_sample.tsv reservoir sampling")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args()


def _genes_from_export(perm_discovery_path: Path) -> Tuple[list, dict]:
    """Genes with at least one CAAS in any cycle (the only genes worth replaying), largest workload first.

    Largest first is longest-processing-time scheduling: starting the biggest gene first keeps
    it from starting late and leaving workers idle behind it (see
    docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md). The row count (file mode) or file size
    (directory mode) is a proxy for a gene's cycle workload, available from the same pass that
    discovers the gene names.

    Returns (genes sorted largest first, {gene: size proxy}); the sizes are passed on as
    gene_sizes of process_all_genes_perms to decide which genes are large enough to split their
    replay across several workers.
    """
    sizes: dict = {}
    if perm_discovery_path.is_file():
        with open(perm_discovery_path, "r") as f:
            header = f.readline()
            if not header:
                return [], {}
            cols = header.rstrip("\n").split("\t")
            try:
                gene_idx = cols.index("gene")
            except ValueError:
                return [], {}
            for line in f:
                # maxsplit=2 keeps the split short: the gene is the second column
                parts = line.split("\t", maxsplit=2)
                if len(parts) > gene_idx:
                    g = parts[gene_idx].strip()
                    if g:
                        sizes[g] = sizes.get(g, 0) + 1
    else:
        # Directory mode: file names are "<alignmentID>.perm_replay.discovery.output", where
        # alignmentID is the base name of the alignment file (for example "GENE.Homo_sapiens"
        # from "GENE.Homo_sapiens.fa"). The gene is the first dot-delimited segment, because
        # the bare symbol is what ensembl_genes, find_gene_alignment (which splits on the
        # first dot) and the meta_caas Gene column use.
        for p in perm_discovery_path.iterdir():
            if p.is_file() and not p.name.startswith("."):
                name = p.name.split(".", 1)[0]
                sizes[name] = sizes.get(name, 0) + p.stat().st_size
    genes = [g for g, _ in sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))]
    return genes, sizes


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    """Replay the permuted labelings (the null) and then b_0 on its own, writing under args.output_dir."""
    args = parse_arguments()
    configure_logging(verbose=args.verbose, log_file=args.log_file)

    logger.info("=" * 80)
    logger.info("CAAS Permulation-Excess Null")
    logger.info("=" * 80)
    for k, v in vars(args).items():
        logger.info(f"  {k}: {v}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genes, gene_sizes = _genes_from_export(Path(args.perm_discovery))
    logger.info(f"Genes with CAAS hits across cycles: {len(genes)}")
    if not genes:
        logger.warning("No genes in export_perm_discovery; writing empty null table")

    cycles = None
    if args.cycles:
        cycles = [c.strip() for c in args.cycles.split(",") if c.strip()]

    # The real labeling (b_0) is replayed as its own single-labeling run into <output-dir>/b0: it goes
    # through the identical code but must never enter the null (n_detected, percent-rank pools and gene
    # removal are all computed over the cycles of one call). With no permuted labeling (N = 0) only b_0 is
    # replayed, which --detail-only allows: the pass that scores the null has nothing to do.
    all_tags = cycles if cycles else sorted(_read_resample_labelings(args.resample_dir))
    b0_tags = [c for c in all_tags if base_cycle(c) == "b_0"]
    null_tags = [c for c in all_tags if base_cycle(c) != "b_0"]

    def run(run_cycles, run_dir):
        return process_all_genes_perms(
            genes=genes,
            alignment_dir=args.alignment_dir,
            tree_file=args.tree,
            perm_discovery_file=args.perm_discovery,
            resample_dir=args.resample_dir,
            taxid_mapping_path=args.taxid_mapping,
            asr_model=args.asr_model,
            asr_cache_dir=args.asr_cache_dir,
            posterior_threshold=args.posterior_threshold,
            workers=args.workers,
            output_dir=run_dir,
            ensembl_genes_file=args.ensembl_genes_file,
            cycles=run_cycles,
            max_tasks_per_child=args.max_tasks_per_child,
            fop_pairs_file=args.fop_pairs,
            gene_lengths_file=args.gene_lengths,
            clust_minlen=args.clust_minlen,
            clust_maxcaas=args.clust_maxcaas,
            gene_filter_mode=args.gene_filter_mode,
            iqr_multiplier=args.iqr_multiplier,
            extreme_percentile=args.extreme_percentile,
            postproc_filter=args.postproc_filter,
            remove_clusters=not args.keep_clusters,
            train_map_dir=args.train_map_dir,
            train_map_suffix=args.train_map_suffix,
            gene_sizes=gene_sizes,
            seed=args.seed,
            detail_only=args.detail_only,
        )

    t0 = time.time()
    if not null_tags and not (b0_tags and args.detail_only):
        raise RuntimeError("[perms] no permuted labelings to replay (only b_0 or nothing was found)")
    if null_tags:
        out_path = run(null_tags, output_dir)
    else:
        logger.info("[perms] no permuted labelings (N = 0): the null has no shards")
        (output_dir / "perm_pos_detail").mkdir(parents=True, exist_ok=True)
        out_path = output_dir
    if b0_tags:
        logger.info(f"b_0: replaying {len(b0_tags)} b_0 labeling(s) -> {output_dir / 'b0'}")
        run(b0_tags, output_dir / "b0")
    logger.info(f"Done in {time.time() - t0:.1f}s → {out_path}")


if __name__ == "__main__":
    main()
