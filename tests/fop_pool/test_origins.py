#!/usr/bin/env python3
# test_origins.py: independent origins of a residue on a species tree.
# PhyloPhere | tests/fop_pool/
"""OriginCounter: carriers that form a clade of carriers are one origin, isolated carriers are one origin each,
tips without a residue do not decide, and the count does not depend on how many carriers a clade has."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"))
from src.asr.tree_parser import parse_newick  # noqa: E402
from src.convergence.origins import OriginCounter  # noqa: E402

TREE = "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,G:1):1);"


def counter(states):
    """OriginCounter on TREE with the residues of `states` ({tip: residue}); a tip not in `states` has no residue."""
    return OriginCounter(parse_newick(TREE), lambda tip, col: states.get(tip))


def main():
    base = {t: "V" for t in "ABCDEFG"}
    # A, B (sisters) and C (its sister D holds another residue): two origins
    c = counter({**base, "A": "I", "B": "I", "C": "I"})
    assert c.count(10, "I", ["A", "B", "C"]) == 2
    # the whole left clade carries I: one origin, however many carriers it has
    c = counter({**base, "A": "I", "B": "I", "C": "I", "D": "I"})
    assert c.count(10, "I", ["A", "B", "C", "D"]) == 1
    assert c.count(10, "I", ["A"]) == 1
    # carriers on both sides of the root: one origin each
    c = counter({**base, "A": "I", "E": "I"})
    assert c.count(10, "I", ["A", "E"]) == 2
    # a tip without a residue does not break a clade: A carries I and its sister B has a gap
    c = counter({"A": "I", "C": "V", "D": "V", "E": "V", "F": "V", "G": "V"})
    assert c.count(10, "I", ["A"]) == 1
    # but it does not make one either: the clade (A,B) has a residue only in A, still one origin; C is not a carrier
    assert c.count(10, "V", ["C", "D"]) == 1
    # a carrier whose clade is the tip itself
    c = counter({**base, "G": "I"})
    assert c.count(10, "I", ["G"]) == 1
    # residues that are not carried give no origin
    assert c.count(10, "L", ["G"]) == 0
    # the same counter serves another column
    c2 = OriginCounter(parse_newick(TREE), lambda tip, col: ("I" if (col == 1 and tip in "AB") else "V"))
    assert c2.count(1, "I", ["A", "B"]) == 1 and c2.count(2, "I", ["A", "B"]) == 0
    print("test_origins: OK")


if __name__ == "__main__":
    main()
