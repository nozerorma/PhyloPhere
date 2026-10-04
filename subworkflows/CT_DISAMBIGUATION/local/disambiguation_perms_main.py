#!/usr/bin/env python3
"""CLI for the CAAS permulation-excess null (genome-wide *excess* null for FCS).

Loads each gene's precomputed ASR posteriors ONCE and replays N permuted phenotype
labelings (the full-pool `export_perm_discovery` from a perm-replay run + the matching
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
from typing import Tuple

sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.utils.gene_wrapper import process_all_genes_perms, _read_resample_labelings
from src.convergence.fop_pool import base_cycle
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
    # ── Gap B: CT_POSTPROC filtering of the null candidate pool ───────────────
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
    """Genes that produced ≥1 CAAS in any cycle (the only genes worth replaying),
    ordered largest-workload-first (LPT scheduling: dispatching the biggest gene
    first keeps it from starting late and stranding idle workers behind it — see
    docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md). Row count (file mode) / file
    size (directory mode) is a free-to-compute proxy for a gene's cycle workload,
    already available from the same iteration that discovers the gene names.
    Returns (genes_sorted_largest_first, {gene: size_proxy}) -- the sizes dict is
    reused downstream (process_all_genes_perms' gene_sizes) to decide which genes
    are big enough to split their replay across multiple workers (Stage 2)."""
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
                # gene is column 2 (index 1), splitting maxsplit=2 avoids splitting remaining columns
                parts = line.split("\t", maxsplit=2)
                if len(parts) > gene_idx:
                    g = parts[gene_idx].strip()
                    if g:
                        sizes[g] = sizes.get(g, 0) + 1
    else:
        # Directory mode: filenames are "<alignmentID>.perm_replay.discovery.output",
        # where alignmentID = Nextflow's f.baseName on the alignment file (e.g.
        # "GENE.Homo_sapiens.fa" -> "GENE.Homo_sapiens"). Splitting on ".perm_replay"
        # left the species suffix attached ("GENE.Homo_sapiens"), which then never
        # matches the bare "GENE" symbols used everywhere else (ensembl_genes,
        # find_gene_alignment's own prefix match, meta_caas "Gene" column) --
        # silently zeroing every gene. Take the first dot-delimited segment instead,
        # mirroring find_gene_alignment's own `path.name.split(".", 1)[0]` convention.
        for p in perm_discovery_path.iterdir():
            if p.is_file() and not p.name.startswith("."):
                name = p.name.split(".", 1)[0]
                sizes[name] = sizes.get(name, 0) + p.stat().st_size
    genes = [g for g, _ in sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))]
    return genes, sizes



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
