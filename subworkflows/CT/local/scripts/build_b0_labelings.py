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

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION" / "local"))
from src.core.labelings import read_design  # noqa: E402  (the one reader of the observed design)


def pair_rows(pairs_path):
    """contrast_hypotheses_pairs.tsv rows, as written to fop_pairs.tsv (they carry the species and PSS)."""
    return list(csv.DictReader(open(pairs_path), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", required=True, help="observed trait file, or the multi-hypothesis directory")
    ap.add_argument("--fop", action="store_true", help="emit the b_0~H<m> hypothesis fan-out + fop_pairs rows")
    ap.add_argument("--labelings-out", required=True, help="labeling rows (append-ready, no header)")
    ap.add_argument("--pairs-out", help="fop_pairs rows (no header); FOP mode only")
    a = ap.parse_args()

    labelings = read_design(a.config)
    if not labelings:
        sys.exit("ERROR: no fg/bg pairs read from the observed design")
    with open(a.labelings_out, "w") as out:
        if a.fop:
            for tag, lab in labelings.items():
                out.write(f"{tag}\t{','.join(lab.fg)}\t{','.join(lab.bg)}\n")
        else:
            # plain mode is the canonical contrast: the single trait file, or H1 of a directory
            lab = labelings.get("b_0") or labelings.get("b_0~H1")
            if lab is None:
                sys.exit("ERROR: the observed design has no canonical (H1) contrast")
            out.write(f"b_0\t{','.join(lab.fg)}\t{','.join(lab.bg)}\n")
    if a.fop and a.pairs_out:
        rows = pair_rows(Path(a.config) / "contrast_hypotheses_pairs.tsv")
        with open(a.pairs_out, "w") as out:
            for r in rows:  # columns as in fop_pairs.tsv: cycle, hypothesis_id, pair, species1, species2, pss_score
                out.write("\t".join(["b_0", r["hypothesis_id"], r["pair"], r["species1"], r["species2"], r["pss_score"]]) + "\n")
    print(f"[b_0] {'FOP: ' + str(len(labelings)) + ' hypotheses' if a.fop else 'plain'}")


if __name__ == "__main__":
    main()
