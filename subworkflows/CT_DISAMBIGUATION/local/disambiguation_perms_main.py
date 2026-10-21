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

For each gene of a manifest it replays the labelings of a resample directory over the gene's alignment in process (the
CAAS kernel of `ct perm-replay`, src/utils/perm_stream.py), loads the gene's precomputed ASR posteriors once and scores
each labeling that has a detection through analyze_gene_disambiguation as the observed pipeline does
(src/utils/gene_wrapper.py, process_all_genes_perms). No discovery export is written. The real labeling (b_0) is
replayed on its own into <output-dir>/b0, so it never enters the null.

Pass A writes one per-(gene, cycle, position) detail shard per gene to
<output-dir>/perm_pos_detail/. With --detail-only the run stops there, because pass B (scoring
against genome-wide pools) runs once over the union of the batches (reaggregate_perm_scores.py).
Without it the same process also writes gene_cycle_scores.tsv, perm_pos_cycle_caas.tsv.gz,
perm_pos_quantiles.tsv and perm_pos_sample.tsv; scoring_caas_perms.R turns gene_cycle_scores.tsv
into the genes x N null matrix behind the FCS p.perm.
With --strict-output pass A also writes <output-dir>/perm_pos_detail_strict/ (and b0/perm_pos_detail_strict/): the
same replay pooled without the hypotheses whose discovery row has conserved pairs, which is the null of an analysis that
keeps only fully divergent hypotheses. Replay, ASR and scoring are shared; only the pooling is repeated.

Called by:  CAAS_CORE_BATCHED (subworkflows/CT/caas_permulation.nf) → disambiguation_perms_main.py --detail-only
Usage:
    python3 disambiguation_perms_main.py --alignment-dir <dir> --tree <newick> --manifest <file> --design <path> \\
        --ct-dir <dir> --replay-args-file <file> --resample-dir <dir> --output-dir <dir> --asr-cache-dir <dir> \\
        [--detail-only] [OPTIONS]
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
from src.utils.perm_stream import ReplaySettings
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
        "--manifest", required=True,
        help="Tab-separated '<gene id>\\t<alignment file name>' per line: the genes to replay; the gene key is the "
             "id up to its first dot",
    )
    p.add_argument("--manifest-alignment-dir", default="alignments",
                   help="Directory of the alignment files named in --manifest (default: alignments)")
    p.add_argument("--design", required=True,
                   help="Trait file, or directory of traitfile_H*.tab, of the run: the alignment is sliced against it")
    p.add_argument("--ct-dir", required=True, help="Directory of `ct` (the `modules` package with the CAAS kernel)")
    p.add_argument("--replay-args-file", required=True,
                   help="File with the `ct perm-replay` thresholds (--patterns, --max_*_gaps, --max_*_miss, "
                        "--miss_pair, --max_conserved, --min_non_overlapping)")
    p.add_argument("--ali-format", default="clustal", help="Alignment format of the manifest files (default: clustal)")
    p.add_argument("--chunk-labelings", type=int, default=None,
                   help="Labelings per replay range of the kernel (default: the kernel's memory budget)")
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
        help="Comma-separated cycle tags to process (default: all cycles of the resample directory)",
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
    p.add_argument("--strict-output", action="store_true",
                   help="Also write perm_pos_detail_strict/ (and b0/perm_pos_detail_strict/): the shards of the same replay "
                        "pooled without the hypotheses whose discovery row has conserved pairs.")
    p.add_argument("--seed", type=int, default=1998,
                   help="Pipeline seed (params.seed): perm_pos_sample.tsv reservoir sampling")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args()


def _read_manifest(path: Path, alignment_dir: Path) -> Tuple[list, dict]:
    """Genes of the manifest, the largest alignment first, and {gene key: alignment path}.

    Largest first is longest-processing-time scheduling: starting the biggest gene first keeps it from starting
    late and leaving workers idle behind it (see docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md). The gene key is
    the manifest id up to its first dot, the symbol that find_gene_alignment and the master Gene column use.
    """
    alignments: dict = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        gene_id, name = line.split("\t")[:2]
        alignments[gene_id.split(".", 1)[0]] = alignment_dir / name
    sizes = {g: a.stat().st_size for g, a in alignments.items()}
    genes = sorted(alignments, key=lambda g: (-sizes[g], g))
    return genes, alignments


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

    genes, gene_alignments = _read_manifest(Path(args.manifest), Path(args.manifest_alignment_dir))
    logger.info(f"Genes to replay: {len(genes)}")
    if not genes:
        logger.warning("No genes in the manifest; writing empty null table")
    replay = ReplaySettings.from_args_file(
        args.replay_args_file, ct_dir=args.ct_dir, design=args.design, ali_format=args.ali_format,
        chunk_labelings=args.chunk_labelings)

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
            gene_alignments=gene_alignments,
            alignment_dir=args.alignment_dir,
            tree_file=args.tree,
            resample_dir=args.resample_dir,
            replay=replay,
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
            seed=args.seed,
            detail_only=args.detail_only,
            strict_output=args.strict_output,
        )

    t0 = time.time()
    if not null_tags and not (b0_tags and args.detail_only):
        raise RuntimeError("[perms] no permuted labelings to replay (only b_0 or nothing was found)")
    if null_tags:
        out_path = run(null_tags, output_dir)
    else:
        logger.info("[perms] no permuted labelings (N = 0): the null has no shards")
        (output_dir / "perm_pos_detail").mkdir(parents=True, exist_ok=True)
        if args.strict_output:
            (output_dir / "perm_pos_detail_strict").mkdir(parents=True, exist_ok=True)
        out_path = output_dir
    if b0_tags:
        logger.info(f"b_0: replaying {len(b0_tags)} b_0 labeling(s) -> {output_dir / 'b0'}")
        run(b0_tags, output_dir / "b0")
    logger.info(f"Done in {time.time() - t0:.1f}s → {out_path}")


if __name__ == "__main__":
    main()
