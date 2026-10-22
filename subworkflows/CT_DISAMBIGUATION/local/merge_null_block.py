#!/usr/bin/env python3
# merge_null_block.py — Merge of the tile tables of one block of cycles of the permutation null.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/

"""
CAAS_CORE_MERGE: merges the tile tables of one block of cycles (src/utils/null_merge.py) into the gene x cycle
statistics of the null, with no ASR replay.

The gene-level and position-level null statistics must be computed with exactly the same formula as the observed side
(scoring_compute.R, observed_core_scores.py), or the FCS p.perm of fcs_enrich.R compares two different quantities. Both
sides use core.scores; the gene statistic of a direction is chosen with --gene-stat.

A cycle belongs to one block and every gene of a cycle is in the block's tiles, so gene removal and the reference pools
of the gene score, which are decided within a cycle, are exact here. The real labeling b_0 is merged the same way as a
block of its own (--include-b0).

Called by:  CAAS_CORE_MERGE (subworkflows/CT/caas_permulation.nf) → merge_null_block.py
Inputs:     --tiles       the null_tables/ directories (units.parquet, scores.parquet) of every tile of the block
            --labelings   the block's labelings file (resample_perms.tab): its cycles are the cycles of the block, the
                          ones with no detection included
            --universe    gene universe (cleaned_background_main.txt), optional
            --gene-lengths / --map-dir, --gene-filter-mode  gene removal, optional
Outputs:    <out>/matrices/   raw matrices and manifest.json (assemble_block_matrices.R turns them into an RDS)
            <out>/positions.parquet, cycle_summary.tsv, cycle_hist.tsv, perm_pos_sample.tsv, perm_pos_quantiles.tsv,
            cycles.txt, counts.json, removed_units.tsv with gene removal

Run from the CT_DISAMBIGUATION/local directory (the process copies it to the work directory first).
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.convergence.fop_pool import base_cycle  # noqa: E402
from src.utils.null_merge import GENE_STATS, merge_block  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def read_roster(labelings_path: Path, by_labeling: bool, include_b0: bool):
    """Cycles of a labelings file (one row per labeling, tag in the first column).

    A cycle is a base cycle ("b_5") unless `by_labeling`: then it is the full tag ("b_5~H3"), the grain of tables whose
    hypotheses were not pooled. The real labeling b_0 is left out of a block of permuted cycles and is the only cycle
    of the block of the real labeling (`include_b0`).
    """
    tags = set()
    with open(labelings_path) as fh:
        for line in fh:
            tag = line.split("\t", 1)[0].strip() if line.strip() else ""
            if not tag:
                continue
            if (base_cycle(tag) == "b_0") != include_b0:
                continue
            tags.add(tag if by_labeling else base_cycle(tag))
    return sorted(tags)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", required=True, nargs="+", type=Path, help="null_tables/ directories of the tiles of the block")
    ap.add_argument("--labelings", required=True, type=Path, help="labelings file of the block (resample_perms.tab)")
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--include-b0", action="store_true", help="the block of the real labeling b_0")
    ap.add_argument("--universe", type=Path, default=None, help="gene universe, one gene per line (a 'Gene' header is ignored)")
    ap.add_argument("--gene-stat", default="size_adj_max", choices=GENE_STATS,
                    help="gene statistic of the matrices named `stat` (the others stay in the parts)")
    ap.add_argument("--seed", type=int, default=1998, help="pipeline seed (params.seed): the per-cycle sample")
    ap.add_argument("--gene-lengths", default=None, help="gene annotation TSV (gene, length ...); enables the dubious/extreme removal")
    ap.add_argument("--gene-filter-mode", default="none", choices=["none", "extreme", "dubious", "both"])
    ap.add_argument("--iqr-multiplier", type=float, default=3.0)
    ap.add_argument("--extreme-percentile", type=float, default=0.99)
    ap.add_argument("--map-dir", "--train-map-dir", type=Path, default=None, help="directory of the trimmer's MAP tables (untrimmed lengths)")
    ap.add_argument("--map-suffix", type=str, default=".map.tsv")
    args = ap.parse_args()

    for d in args.tiles:
        if not (Path(d) / "units.parquet").is_file():
            logger.error("no tile tables under %s", d)
            return 1
    metas = [json.loads((Path(d) / "meta.json").read_text()) for d in args.tiles]
    by_labeling = {bool(m.get("by_labeling")) for m in metas}
    if len(by_labeling) != 1:
        logger.error("the tiles of the block disagree on the grain of their cycles")
        return 1
    cycles = read_roster(args.labelings, by_labeling.pop(), args.include_b0)
    if not cycles:
        logger.error("no cycle in %s", args.labelings)
        return 1

    universe = None
    if args.universe is not None and args.universe.is_file():
        universe = [g.strip() for g in args.universe.read_text().splitlines() if g.strip() and g.strip() != "Gene"]

    removal_on = (bool(args.gene_lengths) or bool(args.map_dir)) and args.gene_filter_mode != "none"
    lengths = {}
    if removal_on:
        from src.core.postproc import load_gene_lengths
        lengths = load_gene_lengths(path=str(args.gene_lengths) if args.gene_lengths else None,
                                    map_dir=str(args.map_dir) if args.map_dir else None, map_suffix=args.map_suffix)
    per_cycle = int(os.environ.get("CAAS_PERMS_SAMPLE_PER_CYCLE", "200"))
    counts = merge_block(args.tiles, args.output_dir, cycles, universe=universe, lengths=lengths,
                         gene_filter_mode=args.gene_filter_mode if removal_on else "none",
                         iqr_multiplier=args.iqr_multiplier, extreme_percentile=args.extreme_percentile,
                         gene_stat=args.gene_stat, seed=args.seed, sample_per_cycle=per_cycle)
    silent = len(cycles) - len(counts["cycles_present"])
    if silent:
        logger.warning("[merge] %d of %d cycles of the block left no row: they count in N as cycles with no signal",
                       silent, len(cycles))
    logger.info("[merge] %d cycles, %d units (%d removed), %d scored positions",
                counts["cycles"], counts["units"], counts["removed_units"], counts["scored_positions"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
