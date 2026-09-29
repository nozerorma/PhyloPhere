#!/usr/bin/env python3
"""Observed position and gene CAAS scores (the b_0 slice), computed with core.scores.

Reads filtered_discovery.tsv and writes the scored quantities scoring_compute.R integrates:

* positions table: Gene, Position, side, CAAS_score (mean of asr_path_score over the schemes
  that scored the position, per side)
* genes table: Gene, gene_caas_score, gene_caas_score_top_all, gene_caas_score_bottom_all
  (size_adj_max against the pool of the same direction; NA when the gene has no scored
  position in it)

Row scope: the five scoring schemes; an asr_path_score that is not numeric counts as missing.
Stdlib only (runs in the same environment as the R step).
"""

import argparse
import csv
import math
import sys
from pathlib import Path

# core.scores lives with the disambiguation core. When this script is staged into a work dir,
# the process copies src/core/scores.py next to it.
_here = Path(__file__).resolve()
for _cand in (_here.parent, _here.parents[3] / "CT_DISAMBIGUATION" / "local"):
    if (_cand / "src" / "core" / "scores.py").exists():
        sys.path.insert(0, str(_cand))
        break
from src.core.scores import direction_values, gene_scores, position_score  # noqa: E402

SCHEMES = ("US", "GS4", "GS3", "GS2", "GS1")


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def read_position_scores(path):
    """{(gene, position, side): {scheme: asr_path_score or None}} for the scoring schemes."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            scheme = row["caap_group"]
            if scheme not in SCHEMES:
                continue
            key = (row["Gene"], int(row["Position"]), row["side"])
            schemes = out.setdefault(key, {})
            if scheme in schemes:
                raise ValueError(f"duplicate scheme {scheme} for {key}: a position is scored once per scheme")
            schemes[scheme] = _num(row["asr_path_score"])
    return out


def _fmt(x):
    return "NA" if x is None else repr(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="filtered_discovery.tsv")
    ap.add_argument("--positions-out", required=True)
    ap.add_argument("--genes-out", required=True)
    a = ap.parse_args()

    try:
        rows = read_position_scores(a.input)
    except ValueError as e:
        print(f"Error: duplicate scheme rows: {e}", file=sys.stderr)
        return 1

    scores = {k: position_score(s) for k, s in rows.items()}
    with open(a.positions_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "Position", "side", "CAAS_score"])
        for (g, p, side), v in sorted(scores.items()):
            w.writerow([g, p, side, _fmt(v)])

    # Reference pools: every scored position of the run, per direction.
    by_pos = {}
    for (g, p, side), v in scores.items():
        by_pos.setdefault((g, p), {})[side] = v
    pools = {d: sorted(v) for d, v in direction_values(by_pos.values()).items()}
    per_gene = {}
    for (g, _p), sides in by_pos.items():
        per_gene.setdefault(g, []).append(sides)
    with open(a.genes_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "gene_caas_score", "gene_caas_score_top_all", "gene_caas_score_bottom_all"])
        for g in sorted(per_gene):
            s = gene_scores(per_gene[g], pools)
            w.writerow([g, _fmt(s["all"]), _fmt(s["top"]), _fmt(s["bottom"])])
    print(f"observed core scores: {len(scores)} position rows, {len(per_gene)} genes; "
          f"pools all={len(pools['all'])} top={len(pools['top'])} bottom={len(pools['bottom'])}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
