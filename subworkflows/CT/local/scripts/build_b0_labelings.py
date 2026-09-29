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


def read_traitfile(path):
    """traitfile*.tab -> (fg, bg) ordered by pair id."""
    fg, bg = {}, {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3 or not row[2].strip().lstrip("-").isdigit():
                continue
            (fg if row[1].strip() == "1" else bg)[int(row[2])] = row[0].strip()
    pairs = sorted(set(fg) & set(bg))
    return [fg[p] for p in pairs], [bg[p] for p in pairs]


def canonical_traitfile(cfg):
    """A directory of traitfile_H*.tab resolves to H1 (the canonical contrast)."""
    cfg = Path(cfg)
    if cfg.is_dir():
        h1 = cfg / "traitfile_H1.tab"
        return h1 if h1.is_file() else sorted(cfg.glob("*.tab"))[0]
    return cfg


def read_hypotheses(pairs_path):
    """contrast_hypotheses_pairs.tsv -> ({hyp: (fg, bg)} in file order, raw rows)."""
    rows = list(csv.DictReader(open(pairs_path), delimiter="\t"))
    hyps = {}
    for r in rows:
        hyps.setdefault(r["hypothesis_id"], []).append(r)
    labelings = {}
    for h, rs in hyps.items():
        rs = sorted(rs, key=lambda r: int(float(r["pair"])))
        labelings[h] = ([r["species1"] for r in rs], [r["species2"] for r in rs])
    return labelings, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", required=True, help="observed trait file, or the multi-hypothesis directory")
    ap.add_argument("--fop", action="store_true", help="emit the b_0~H<m> hypothesis fan-out + fop_pairs rows")
    ap.add_argument("--labelings-out", required=True, help="labeling rows (append-ready, no header)")
    ap.add_argument("--pairs-out", help="fop_pairs rows (no header); FOP mode only")
    a = ap.parse_args()

    if a.fop:
        pairs_path = Path(a.config) / "contrast_hypotheses_pairs.tsv"
        if not pairs_path.is_file():
            sys.exit(f"ERROR: FOP b_0 needs {pairs_path}")
        labelings, rows = read_hypotheses(pairs_path)
        with open(a.labelings_out, "w") as out:
            for h, (fg, bg) in labelings.items():
                out.write(f"b_0~{h}\t{','.join(fg)}\t{','.join(bg)}\n")
        if a.pairs_out:
            with open(a.pairs_out, "w") as out:
                for r in rows:  # columns as in fop_pairs.tsv: cycle, hypothesis_id, pair, species1, species2, pss_score
                    out.write("\t".join(["b_0", r["hypothesis_id"], r["pair"], r["species1"],
                                         r["species2"], r["pss_score"]]) + "\n")
        print(f"[b_0] FOP: {len(labelings)} hypotheses")
    else:
        fg, bg = read_traitfile(canonical_traitfile(a.config))
        if not fg:
            sys.exit("ERROR: no fg/bg pairs read from the observed trait file")
        with open(a.labelings_out, "w") as out:
            out.write(f"b_0\t{','.join(fg)}\t{','.join(bg)}\n")
        print(f"[b_0] plain: {len(fg)} pairs")


if __name__ == "__main__":
    main()
