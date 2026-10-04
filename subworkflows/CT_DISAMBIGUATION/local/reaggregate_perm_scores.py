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
import gzip
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.convergence.fop_pool import base_cycle  # noqa: E402
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


def read_roster(labelings_path: Path, by_labeling: bool = False):
    """Cycles of a labelings file (one row per labeling, tag in the first column), the real labeling b_0 left out.

    A cycle is a base cycle ("b_5") unless `by_labeling`: then it is the full tag ("b_5~H3"), the grain of a detail
    whose hypotheses were not pooled.
    """
    tags = set()
    with open(labelings_path) as fh:
        for line in fh:
            tag = line.split("\t", 1)[0].strip() if line.strip() else ""
            if tag and base_cycle(tag) != "b_0":
                tags.add(tag if by_labeling else base_cycle(tag))
    return sorted(tags)


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
    ap.add_argument("--cycles-from", type=Path, default=None,
                    help="labelings file the cycles were replayed from (resample_perms.tab). N of the null is the number "
                         "of its cycles, not the number of cycles that left a row: writes cycle_roster.txt for "
                         "scoring_caas_perms.R, warns about the cycles without a row, and refuses rows of a cycle it "
                         "does not list")
    ap.add_argument("--empty-null", action="store_true",
                    help="the null has no permuted labeling (N = 0): write the null tables empty. Refused when --detail "
                         "holds shards, and without it a null with no shard is an error (a replay that lost every hit "
                         "must not pass for an empty null)")
    ap.add_argument("--iqr-multiplier", type=float, default=3.0)
    ap.add_argument("--extreme-percentile", type=float, default=0.99)
    args = ap.parse_args()

    if not (args.detail.is_dir() or args.detail.is_file()):
        logger.error("detail path not found: %s", args.detail)
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.empty_null:
        shards = sorted(args.detail.glob("*.tsv.gz")) if args.detail.is_dir() else [args.detail]
        if shards:
            logger.error("--empty-null given but %s holds %d shard(s): the null is not empty", args.detail, len(shards))
            return 1
        # The tables go through the same writer as a normal run, over a detail with a header and no rows, so their
        # columns cannot drift from it.
        empty_detail = args.output_dir / ".empty_perm_pos_detail.tsv.gz"
        with gzip.open(empty_detail, "wt") as fh:
            fh.write("")
        try:
            _finalize_perm_scores(detail_path=empty_detail, output_dir=args.output_dir, cycle_tags=[], removed=set(),
                                  seed=args.seed, remove_clusters=not args.keep_clusters)
        finally:
            empty_detail.unlink()
        logger.warning("[reaggregate] empty null (no permuted labeling): no gene x cycle scores were written")
        return 0

    logger.info("[reaggregate] scanning %s", args.detail)
    cycle_tags, n_rows = scan_detail(args.detail)
    logger.info("[reaggregate] %d rows across %d cycles", n_rows, len(cycle_tags))
    if not cycle_tags:
        logger.error("no cycles found in detail file; nothing to do")
        return 1

    if args.cycles_from is not None:
        # The roster is kept at the grain of the detail: base cycles when the hypotheses were pooled (the production
        # case), the labeling tags otherwise.
        roster = read_roster(args.cycles_from, by_labeling=any("~" in c for c in cycle_tags))
        present = set(cycle_tags)
        stray = sorted(present - set(roster))
        if stray:
            logger.error("detail rows of %d cycle(s) that are not in %s: %s", len(stray), args.cycles_from, ", ".join(stray[:5]))
            return 1
        silent = sorted(set(roster) - present)
        if silent:
            logger.warning("[reaggregate] %d of %d replayed cycles left no row (%s%s): they count in N as cycles with no signal",
                           len(silent), len(roster), ", ".join(silent[:5]), ", ..." if len(silent) > 5 else "")
        (args.output_dir / "cycle_roster.txt").write_text("".join(f"{c}\n" for c in roster))

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
