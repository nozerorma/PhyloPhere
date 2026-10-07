#!/usr/bin/env python3
# reaggregate_perm_scores.py — Pass B of the permulation null: gene x cycle scores from the perm_pos_detail shards.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/

"""
CAAS_CORE_MERGE: Scores the permulation null (pass B) from the per-gene perm_pos_detail shards, with no ASR replay.

Why it exists
-------------
The gene-level and position-level null statistics must be computed with exactly the same formula as
the observed side (scoring_compute.R), or the FCS p.perm of fcs_enrich.R compares two different
quantities. The scoring pass is therefore the aggregation of gene_wrapper.py (_finalize_perm_scores),
called here and not reimplemented: the gene score is F(max)^n over a pool of heavily tied values, so a
difference of 1e-16 in how the per-position sum is accumulated can flip a tie boundary, and the ^n then
amplifies it.

The detail is the output of pass A (disambiguation_perms_main.py --detail-only): one gz shard per gene in
perm_pos_detail/ (or a single concatenated perm_pos_detail.tsv.gz), holding every (Gene, cycle, Position,
caap_group, asr_path_score, n_detected, clust, side) row that pass B needs. Pass B takes minutes, whereas the
ASR replay that produces the detail takes hours. It can also be run alone on the detail of an earlier run to
rebuild its null tables, then followed by scoring_caas_perms.R for caas_perms.rds:

    python3 reaggregate_perm_scores.py \\
        --detail  <run>/caas_permulation/perm_pos_detail \\
        --output-dir <run>/caas_permulation
    Rscript subworkflows/SCORING/local/src/scoring_caas_perms.R \\
        --gene-cycle-scores <run>/caas_permulation/gene_cycle_scores.tsv \\
        --output            <run>/caas_permulation/caas_perms.rds

Run from the CT_DISAMBIGUATION/local directory (the process copies it to the work directory first).

Called by:  CAAS_CORE_MERGE (subworkflows/CT/caas_permulation.nf) → reaggregate_perm_scores.py
Inputs:     --detail  perm_pos_detail/ shard directory or a concatenated perm_pos_detail.tsv.gz
            --cycles-from  optional labelings file (resample_perms.tab) that fixes the cycle roster N
            --gene-lengths, --gene-filter-mode  optional dubious/extreme gene removal
Outputs:    <output-dir>/gene_cycle_scores.tsv, perm_pos_cycle_caas.tsv.gz, perm_pos_quantiles.tsv,
            perm_pos_sample.tsv; cycle_roster.txt with --cycles-from; removed_units.tsv with gene removal
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import gzip
import logging
import sys
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.convergence.fop_pool import base_cycle  # noqa: E402
from src.core.scores import AGGREGATIONS, DEFAULT_AGGREGATION  # noqa: E402
from src.utils.gene_wrapper import (  # noqa: E402
    _cycle_gene_removal_from_detail,
    _finalize_perm_scores,
    iter_detail_rows,
    write_removed_units,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


# ── Detail scan and roster ────────────────────────────────────────────────────


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


# ── CLI ───────────────────────────────────────────────────────────────────────


def main() -> int:
    """Run pass B over the detail and write the null tables; returns the process exit status."""
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
    ap.add_argument("--score-aggregation", default=DEFAULT_AGGREGATION, choices=AGGREGATIONS,
                    help="scheme aggregation of the position score (params.caas_score_aggregation); "
                         "the observed score must use the same")
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
                                  seed=args.seed, remove_clusters=not args.keep_clusters,
                                  aggregation=args.score_aggregation)
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
        aggregation=args.score_aggregation,
    )
    logger.info("[reaggregate] wrote %s", args.output_dir / "gene_cycle_scores.tsv")
    # _finalize_perm_scores also writes perm_pos_cycle_caas.tsv.gz: the per-cycle position
    # CAAS scores behind the position-level p.emp of scoring_compute.R.
    logger.info("[reaggregate] wrote %s", args.output_dir / "perm_pos_cycle_caas.tsv.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
