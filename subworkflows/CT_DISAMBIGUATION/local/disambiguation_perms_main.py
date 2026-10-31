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

For each gene of a manifest it replays the labelings of a resample directory (the cycles of one block, or the real
labeling b_0) over the gene's alignment in process (the CAAS kernel of `ct perm-replay`, src/utils/perm_stream.py), loads
the gene's precomputed ASR posteriors once and scores each labeling that has a detection through analyze_gene_disambiguation
as the observed pipeline does. No discovery export is written. Each scored cycle is reduced as it is produced to a
(gene, cycle) unit and the scored positions of the cycle (src/utils/null_tables.py, src/utils/gene_wrapper.py,
process_all_genes_perms), and the tables of the genes of the batch are written to <output-dir>/null_tables/.
With --strict-output the tables of the same replay pooled without the hypotheses whose discovery row has residue-type overlap
are also written to <output-dir>/null_tables_strict/, which is the null of an analysis that keeps only fully divergent
hypotheses. Replay, ASR and scoring are shared; only the pooling is repeated.

The statistics that need every gene of a cycle (gene removal, the reference pools of the gene score) are computed once
per block of cycles from the tables of all the batches (merge_null_block.py).

Called by:  CAAS_CORE_BATCHED (subworkflows/CT/caas_permulation.nf)
Usage:
    python3 disambiguation_perms_main.py --alignment-dir <dir> --tree <newick> --manifest <file> --design <path> \\
        --ct-dir <dir> --replay-args-file <file> --resample-dir <dir> --output-dir <dir> --asr-cache-dir <dir> [OPTIONS]
"""

# ── Standard library ──────────────────────────────────────────────────────────
import sys
import argparse
import logging
import time
from pathlib import Path
from typing import Tuple

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.utils.gene_wrapper import process_all_genes_perms, _read_resample_labelings
from src.utils.perm_stream import ReplaySettings
from src.core.universe import read_dir as read_universe_dir, STRICT_SUFFIX
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
             "complete id, including embedded dots",
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
    p.add_argument("--universe-dir", required=True, type=Path,
                   help="Directory of the *.b0.universe files (ct perm-replay --export_b0_universe) of the manifest genes: the "
                        "positional universe, the only positions the null scores; with --strict-output also the *.b0.universe_strict files")
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
                   help="Flag the cluster trains of every cycle (CT_POSTPROC filter) so the null follows the same "
                        "filter as the observed filtered_discovery.tsv.")
    p.add_argument("--clust-minlen", type=int, default=3,
                   help="CT_FILTER minlen (params.filter_minlen; default 3)")
    p.add_argument("--clust-min-positions", type=int, default=None, help="Minimum candidates per train; default follows clust-minlen")
    p.add_argument("--clust-maxcaas", type=float, default=0.7,
                   help="CT_FILTER maxcaas density (params.filter_maxcaas; default 0.7)")
    p.add_argument("--keep-clusters", action="store_true",
                   help="Keep cluster-train positions in the scores (params.remove_caas_clusters "
                        "false); they still count for the dubious-gene test.")
    p.add_argument("--train-map-dir", default=None,
                   help="Directory of the trimmer's per-gene MAP tables. Cluster trains then measure their span "
                        "in untrimmed alignment columns (a gene without a MAP keeps trimmed coordinates). "
                        "Needs --postproc-filter; the observed filter must use the same directory.")
    p.add_argument("--train-map-suffix", default=".map.tsv",
                   help="File-name tail of the MAP tables (default .map.tsv)")
    p.add_argument("--strict-output", action="store_true",
                   help="Also write null_tables_strict/: the tables of the same replay pooled without the hypotheses "
                        "whose discovery row has residue-type overlap.")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    args = p.parse_args()
    import math
    if (args.clust_minlen < 1 or (args.clust_min_positions is not None and args.clust_min_positions < 1)
            or not math.isfinite(args.clust_maxcaas) or not 0 <= args.clust_maxcaas <= 1):
        p.error('Invalid cluster span, minimum candidate count or density')
    return args


def _read_manifest(path: Path, alignment_dir: Path) -> Tuple[list, dict]:
    """Genes of the manifest, the largest alignment first, and {gene key: alignment path}.

    Largest first is longest-processing-time scheduling: starting the biggest gene first keeps it from starting
    late and leaving workers idle behind it (see docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md). The gene key is
    the complete manifest identifier, including embedded dots.
    """
    alignments: dict = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        gene_id, name = line.split("\t")[:2]
        if gene_id in alignments:
            raise ValueError(f'Duplicate gene in manifest: {gene_id}')
        alignments[gene_id] = alignment_dir / name
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

    t0 = time.time()
    tags = cycles if cycles else sorted(_read_resample_labelings(args.resample_dir))
    if not tags:
        raise RuntimeError("[perms] no labelings to replay in the resample directory")
    out_path = process_all_genes_perms(
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
        output_dir=output_dir,
        ensembl_genes_file=args.ensembl_genes_file,
        cycles=cycles,
        max_tasks_per_child=args.max_tasks_per_child,
        fop_pairs_file=args.fop_pairs,
        clust_minlen=args.clust_minlen,
        clust_min_positions=args.clust_min_positions,
        clust_maxcaas=args.clust_maxcaas,
        postproc_filter=args.postproc_filter,
        remove_clusters=not args.keep_clusters,
        train_map_dir=args.train_map_dir,
        train_map_suffix=args.train_map_suffix,
        strict_output=args.strict_output,
        universe=read_universe_dir(args.universe_dir),
        universe_strict=read_universe_dir(args.universe_dir, STRICT_SUFFIX) if args.strict_output else None,
    )
    logger.info(f"Done in {time.time() - t0:.1f}s → {out_path}")


if __name__ == "__main__":
    main()
