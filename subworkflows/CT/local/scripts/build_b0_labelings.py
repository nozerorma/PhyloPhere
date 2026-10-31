#!/usr/bin/env python3
"""Labelings of the real (observed) design, tagged ``b_0``, in the null's own formats.

Emits the rows the permulation null consumes, so the observed labeling can be
replayed through exactly the same path as the permuted ones:

  plain mode : one row  ``b_0 <TAB> fg_csv <TAB> bg_csv``           (resample_*.tab shape)
  FOP mode   : one row per hypothesis ``b_0~H<m> <TAB> fg <TAB> bg`` (fop_labelings.tab shape)
               plus the matching ``fop_pairs.tsv`` rows (cycle = b_0)

Inputs come from the observed contrast selection: a trait file (species, 0/1, pair id;
1 = foreground) or the multi-hypothesis directory holding ``traitfile_H*.tab`` and
``contrast_hypotheses_pairs.tsv``. Pairing is by index: fg[k] <-> bg[k] is pair k+1.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.design_io import read_design, pair_metadata
from prepare_labelings import labelings_of_design


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", required=True, help="observed trait file, or the multi-hypothesis directory")
    ap.add_argument("--fop", action="store_true", help="emit the b_0~H<m> hypothesis fan-out + fop_pairs rows")
    ap.add_argument("--labelings-out", required=True, help="labeling rows (append-ready, no header)")
    ap.add_argument("--pairs-out", help="fop_pairs rows (no header); FOP mode only")
    a = ap.parse_args()

    design = read_design(a.config, a.fop)
    labelings = labelings_of_design(design, a.fop)
    with open(a.labelings_out, "w") as out:
        for tag, (fg, bg) in labelings.items():
            out.write(f"{tag}\t{','.join(fg)}\t{','.join(bg)}\n")
    if a.fop and a.pairs_out:
        p = Path(a.config)
        metadata = (p if p.is_dir() else p.parent) / "contrast_hypotheses_pairs.tsv"
        weights = pair_metadata(metadata, labelings, observed=True)
        with open(a.pairs_out, "w") as out:
            for (cycle, hyp, pair), weight in weights.items():
                sides = design[hyp][pair]
                out.write(f"{cycle}\t{hyp}\t{pair}\t{sides['1']}\t{sides['0']}\t{weight}\n")
    print(f"[b_0] {'FOP: ' + str(len(labelings)) + ' hypotheses' if a.fop else 'plain'}")


if __name__ == "__main__":
    main()
