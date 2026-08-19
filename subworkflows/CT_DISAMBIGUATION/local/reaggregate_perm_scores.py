#!/usr/bin/env python3
"""Rebuild gene_cycle_scores.tsv from an existing perm_pos_detail.tsv.gz.

Why this exists
---------------
The gene-level and position-level null statistics must be computed with exactly
the same formula as the observed side (scoring_compute.R), or the FCS p.perm in
fcs_enrich.R compares two different quantities and silently goes wrong. So any
change to the scoring formula obliges a rebuild of caas_perms.rds.

Rebuilding it from scratch means re-running CAAS_PERMS_DISAMBIGUATE, whose cost
is the ASR replay across every gene x labeling -- hours. But perm_pos_detail.tsv.gz
already holds every (Gene, cycle, Position, caap_group, asr_path_score,
n_detected, ct, cb) row the aggregation needs, so re-scoring needs no ASR at all
(see caas_permulation.nf, which publishes the detail file for exactly this).
This script does that re-scoring in minutes.

Pipe the output back through scoring_caas_perms.R to regenerate caas_perms.rds:

    python3 reaggregate_perm_scores.py \\
        --detail  <run>/caas_permulation/perm_pos_detail.tsv.gz \\
        --output-dir <run>/caas_permulation
    Rscript subworkflows/SCORING/local/src/scoring_caas_perms.R \\
        --gene-cycle-scores <run>/caas_permulation/gene_cycle_scores.tsv \\
        --output            <run>/caas_permulation/caas_perms.rds

IMPORTANT: this deliberately calls gene_wrapper's own build_percent_rank_lookup
and _finalize_perm_scores rather than reimplementing the aggregation. The gene
score is F(max)^n over a pool of heavily tied values, so a difference of 1e-16 in
how the per-position sum is accumulated can flip a tie boundary, and the ^n then
amplifies it -- a pandas reimplementation was measured drifting up to 5.9e-3 from
the streaming one. Same code path, or the null stops matching the observed side.

Run from the CT_DISAMBIGUATION/local directory (as the Nextflow process does).
"""
import argparse
import csv
import gzip
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.utils.gene_wrapper import (  # noqa: E402
    build_percent_rank_lookup,
    _finalize_perm_scores,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def scan_detail(detail_path: Path):
    """One pass over the detail file for the two things pass A would have produced.

    Returns (cycle_tags, hist_by_cycle) where hist_by_cycle[cycle][n_detected] is
    the count of (gene, position, scheme) candidates that cycle discovered at that
    replication level -- the input build_percent_rank_lookup expects.
    """
    hist_by_cycle = {}
    n_rows = 0
    with gzip.open(detail_path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cyc = row["cycle"]
            d = int(row["n_detected"])
            per_cycle = hist_by_cycle.setdefault(cyc, {})
            per_cycle[d] = per_cycle.get(d, 0) + 1
            n_rows += 1
    # Sorted so the emitted gene x cycle rows keep a deterministic column order.
    return sorted(hist_by_cycle), hist_by_cycle, n_rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detail", required=True, type=Path,
                    help="perm_pos_detail.tsv.gz from a previous CAAS_PERMS_DISAMBIGUATE run")
    ap.add_argument("--output-dir", required=True, type=Path,
                    help="directory to write gene_cycle_scores.tsv (and the sample/quantile files)")
    args = ap.parse_args()

    if not args.detail.is_file():
        logger.error("detail file not found: %s", args.detail)
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[reaggregate] scanning %s", args.detail)
    cycle_tags, hist_by_cycle, n_rows = scan_detail(args.detail)
    logger.info("[reaggregate] %d rows across %d cycles", n_rows, len(cycle_tags))
    if not cycle_tags:
        logger.error("no cycles found in detail file; nothing to do")
        return 1

    rank_lookup = build_percent_rank_lookup(hist_by_cycle)
    _finalize_perm_scores(
        detail_path=args.detail,
        output_dir=args.output_dir,
        cycle_tags=cycle_tags,
        rank_lookup=rank_lookup,
    )
    logger.info("[reaggregate] wrote %s", args.output_dir / "gene_cycle_scores.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
