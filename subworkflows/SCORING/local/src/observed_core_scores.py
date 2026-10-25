#!/usr/bin/env python3
# observed_core_scores.py — Observed position and gene CAAS scores (the b_0 slice) from filtered_discovery.tsv.
# PhyloPhere | subworkflows/SCORING/local/src/

"""
ObservedCoreScores: scores the observed labeling with the same functions that score the
permulation null (core.scores), so observed and null statistics are built by one formula.

Called by:  SCORING_COMPUTE Nextflow process (scoring_compute.nf → observed_core_scores.py);
            its two outputs are read by scoring_compute.R (--core_positions, --core_genes)
Inputs:     --input  filtered_discovery.tsv, one row per (Gene, Position, side) with an asr_path_score
                     column; an asr_path_score that is not numeric counts as missing
Outputs:    --positions-out  TSV [Gene, Position, side, CAAS_score]; one row per detected position side,
                             CAAS_score being its asr_path_score (core.scores)
            --genes-out      TSV [Gene, gene_caas_score, gene_caas_score_top,
                             gene_caas_score_bottom]; the gene statistic chosen with --gene-stat
                             (max only), NA when the gene has no scored position
                             in that direction

Needs only the standard library plus src/core/scores.py, which the process stages next to it.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

# core.scores lives with the disambiguation core. Staged in a work dir, the script finds the
# copy of src/core/scores.py that the process places next to it; run from the repository, it
# falls back to CT_DISAMBIGUATION/local.
_here = Path(__file__).resolve()
for _cand in (_here.parent, _here.parents[3] / "CT_DISAMBIGUATION" / "local"):
    if (_cand / "src" / "core" / "scores.py").exists():
        sys.path.insert(0, str(_cand))
        break
from src.core.scores import GENE_STATS, gene_scores, position_score  # noqa: E402


def _num(x):
    """float(x), or None when x is not numeric or is NaN."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def read_position_scores(path):
    """{(gene, position, side): asr_path_score or None}.

    A position side is scored once; a repeated (gene, position, side) row raises ValueError
    instead of being averaged silently.
    """
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["Gene"], int(row["Position"]), row["side"])
            if key in out:
                raise ValueError(f"duplicate row for {key}: a position side is scored once")
            out[key] = _num(row["asr_path_score"])
    return out


def _fmt(x):
    """Text of a score: "NA" for None, else repr() (shortest string that round-trips the float)."""
    return "NA" if x is None else repr(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="filtered_discovery.tsv")
    ap.add_argument("--positions-out", required=True)
    ap.add_argument("--genes-out", required=True)
    ap.add_argument("--gene-stat", default="max", choices=GENE_STATS,
                    help="gene statistic of the gene scores; the permulation null must hold the same one")
    a = ap.parse_args()

    try:
        rows = read_position_scores(a.input)
    except ValueError as e:
        print(f"Error: duplicate rows: {e}", file=sys.stderr)
        return 1

    scores = {k: position_score(v) for k, v in rows.items()}
    with open(a.positions_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "Position", "side", "CAAS_score"])
        for (g, p, side), v in sorted(scores.items()):
            w.writerow([g, p, side, _fmt(v)])

    # Group positions by gene; global takes the best side of each position.
    by_pos = {}
    for (g, p, side), v in scores.items():
        by_pos.setdefault((g, p), {})[side] = v
    pools = {d: () for d in ("all", "top", "bottom")}  # ignored by the raw maximum
    per_gene = {}
    for (g, _p), sides in by_pos.items():
        per_gene.setdefault(g, []).append(sides)
    with open(a.genes_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "gene_caas_score", "gene_caas_score_top", "gene_caas_score_bottom"])
        for g in sorted(per_gene):
            s = gene_scores(per_gene[g], pools, a.gene_stat)
            w.writerow([g, _fmt(s["all"]), _fmt(s["top"]), _fmt(s["bottom"])])
    print(f"observed core scores: {len(scores)} position rows, {len(per_gene)} genes; "
          "gene_stat=max", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
