#!/usr/bin/env python3
"""Rebuild gene_cycle_scores.tsv from an existing perm_pos_detail.tsv.gz.

Why this exists
---------------
The gene-level and position-level null statistics must be computed with exactly
the same formula as the observed side (scoring_compute.R), or the FCS p.perm in
fcs_enrich.R compares two different quantities and silently goes wrong. So any
change to the scoring formula obliges a rebuild of caas_perms.rds.

Rebuilding it from scratch means re-running the ASR replay (CAAS_CORE_BATCHED), whose cost
is the ASR replay across every gene x labeling -- hours. But perm_pos_detail/
(one gz shard per gene; a legacy run may instead have a single concatenated
perm_pos_detail.tsv.gz) already holds every (Gene, cycle, Position, caap_group,
asr_path_score, n_detected, side) row the aggregation needs, so re-scoring
needs no ASR at all (see caas_permulation.nf, which publishes the detail
shards for exactly this). This script does that re-scoring in minutes.

Pipe the output back through scoring_caas_perms.R to regenerate caas_perms.rds:

    python3 reaggregate_perm_scores.py \\
        --detail  <run>/caas_permulation/perm_pos_detail \\
        --output-dir <run>/caas_permulation
    Rscript subworkflows/SCORING/local/src/scoring_caas_perms.R \\
        --gene-cycle-scores <run>/caas_permulation/gene_cycle_scores.tsv \\
        --output            <run>/caas_permulation/caas_perms.rds

IMPORTANT: this deliberately calls gene_wrapper's own _finalize_perm_scores rather than reimplementing the aggregation. The gene
score is F(max)^n over a pool of heavily tied values, so a difference of 1e-16 in
how the per-position sum is accumulated can flip a tie boundary, and the ^n then
amplifies it -- a pandas reimplementation was measured drifting up to 5.9e-3 from
the streaming one. Same code path, or the null stops matching the observed side.

Run from the CT_DISAMBIGUATION/local directory (as the Nextflow process does).
"""
import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.utils.gene_wrapper import (  # noqa: E402
    _cycle_gene_removal_from_detail,
    _finalize_perm_scores,
    iter_detail_rows,
    write_removed_units,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def scan_detail(detail_path: Path):
    """One pass over the detail file (or shard directory): the cycles present and the row count.

    Returns (cycle_tags, n_rows); cycle_tags is sorted so the emitted gene x cycle rows keep a deterministic
    column order.
    """
    cycles = set()
    n_rows = 0
    for row in iter_detail_rows(detail_path):
        cycles.add(row["cycle"])
        n_rows += 1
    return sorted(cycles), n_rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detail", required=True, type=Path,
                    help="perm_pos_detail/ shard directory (current layout) or a legacy "
                         "concatenated perm_pos_detail.tsv.gz, from a previous "
                         "CAAS_CORE_BATCHED run")
    ap.add_argument("--output-dir", required=True, type=Path,
                    help="directory to write gene_cycle_scores.tsv (and the sample/quantile files)")
    ap.add_argument("--seed", type=int, default=1998,
                    help="pipeline seed (params.seed): perm_pos_sample.tsv reservoir sampling")
    # Gene removal is genome-wide (per-cycle IQR / density thresholds over all genes), so it
    # can only be computed here, on the merged detail, never inside a per-batch worker.
    ap.add_argument("--gene-lengths", default=None,
                    help="gene annotation TSV (gene, length ...); enables the dubious/extreme gene removal")
    ap.add_argument("--gene-filter-mode", default="none", choices=["none", "extreme", "dubious", "both"])
    ap.add_argument("--keep-clusters", action="store_true",
                    help="keep cluster-train positions in the scored pool (params.remove_caas_clusters false)")
    ap.add_argument("--iqr-multiplier", type=float, default=3.0)
    ap.add_argument("--extreme-percentile", type=float, default=0.99)
    args = ap.parse_args()

    if not (args.detail.is_dir() or args.detail.is_file()):
        logger.error("detail path not found: %s", args.detail)
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[reaggregate] scanning %s", args.detail)
    cycle_tags, n_rows = scan_detail(args.detail)
    logger.info("[reaggregate] %d rows across %d cycles", n_rows, len(cycle_tags))
    if not cycle_tags:
        logger.error("no cycles found in detail file; nothing to do")
        return 1

    removed = set()
    if args.gene_lengths and args.gene_filter_mode != "none":
        from src.core.postproc import load_gene_lengths
        removed = _cycle_gene_removal_from_detail(
            args.detail, load_gene_lengths(args.gene_lengths), args.gene_filter_mode,
            args.iqr_multiplier, args.extreme_percentile)
        write_removed_units(args.output_dir / "removed_units.tsv", removed)
        logger.info("[reaggregate] gene removal (%s): %d (cycle, group, gene) units",
                    args.gene_filter_mode, len(removed))

    _finalize_perm_scores(
        detail_path=args.detail,
        output_dir=args.output_dir,
        cycle_tags=cycle_tags,
        removed=removed,
        seed=args.seed,
        remove_clusters=not args.keep_clusters,
    )
    logger.info("[reaggregate] wrote %s", args.output_dir / "gene_cycle_scores.tsv")
    # V3-4a: _finalize_perm_scores also emits perm_pos_cycle_caas.tsv.gz (the
    # per-cycle CAAS numerator/denominator behind scoring_compute.R's p.emp) --
    # a rebuild from an existing run's detail shards regenerates it for free.
    logger.info("[reaggregate] wrote %s", args.output_dir / "perm_pos_cycle_caas.tsv.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
