#!/usr/bin/env python3
# =============================================================================
# boot_fop_collapse_equivtest.py — FOP (Gap A) base-cycle collapse
# =============================================================================
# Under params.caas_perms_fop the OBSERVED bootstrap resamples over
# fop_labelings.tab, whose cycle tags are "<base>~H<m>" (one row per null cycle
# and fanned Dunn-independent alternative hypothesis). recovery_boot must be
# reported in BASE-CYCLE units: a base cycle HITS an observed (Gene@Position,
# scheme) iff ANY of its ~H<m> labelings calls a CAAS there.
#
# Part 1: unit test of collapse_fop_hits_by_base() in isolation.
# Part 2: end-to-end ct-bootstrap-equivalent invocation over a hand-built
#         fop_labelings.tab + toy alignment; asserts the collapsed count is
#         3/4 (base cycles b_1,b_2,b_4 hit; b_3 does not) — NOT 4/12.
#
# Run:  micromamba run -n phylophere python \
#         subworkflows/CT/local/modules/boot_fop_collapse_equivtest.py
# =============================================================================

import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_LOCAL = os.path.dirname(_HERE)
if _LOCAL not in sys.path:
    sys.path.insert(0, _LOCAL)

from modules.boot import collapse_fop_hits_by_base, boot_on_single_alignment

FAILS = 0


def _check(label, cond, detail=""):
    global FAILS
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    if not cond:
        FAILS += 1


class Sliced:
    def __init__(self, genename, species, d):
        self.genename = genename
        self.species = species
        self.d = d


# --------------------------------------------------------------------------
# Part 1 — pure collapse function
# --------------------------------------------------------------------------
def test_collapse_unit():
    print("Part 1: collapse_fop_hits_by_base (isolated)")
    all_labelings = []
    for b in (1, 2, 3, 4):
        for h in (1, 2, 3):
            all_labelings.append(f"b_{b}~H{h}")

    # b_1 hits only via ~H2 ; b_2 via ~H1 and ~H3 ; b_3 no hit ; b_4 via ~H2
    per_key_hits = {
        "GENEX@10": ["b_1~H2", "b_2~H1", "b_2~H3", "b_4~H2"],
        "GENEX@99": [],  # never hit
    }
    collapsed, base_total = collapse_fop_hits_by_base(per_key_hits, all_labelings)
    _check("base_total is 4 distinct base cycles", base_total == 4, f"got {base_total}")
    _check("P collapses 4 labeling-hits -> 3 base cycles",
           collapsed["GENEX@10"] == 3, f"got {collapsed['GENEX@10']}")
    _check("empty hit list -> 0", collapsed["GENEX@99"] == 0)

    # tag with no ~H suffix is its own base cycle
    c2, bt2 = collapse_fop_hits_by_base({"k": ["b_7", "b_7~H2"]}, ["b_7", "b_7~H2", "b_8~H1"])
    _check("bare tag folds with its ~H siblings", c2["k"] == 1 and bt2 == 2,
           f"got {c2['k']}, {bt2}")


# --------------------------------------------------------------------------
# Part 2 — end-to-end through boot_on_single_alignment(fop_mode=True)
# --------------------------------------------------------------------------
def test_end_to_end(vectorize):
    os.environ["CT_BOOTSTRAP_VECTORIZE"] = "1" if vectorize else "0"
    import importlib
    import modules.boot_vec as bv
    import modules.boot as boot
    importlib.reload(bv)
    importlib.reload(boot)

    tmp = tempfile.mkdtemp(prefix="boot_fop_")
    # s1..s4 carry 'A' at position 10 ; s5..s8 carry 'W'.
    species = [f"s{i}" for i in range(1, 9)]
    pos = {sp: ("A@10" if int(sp[1:]) <= 4 else "W@10") for sp in species}
    sliced = Sliced("GENEX", species, [pos])

    HIT = ("s1,s2", "s5,s6")      # fg all A, bg all W  -> CAAS (pattern 1)
    NOHIT = ("s1,s2", "s3,s4")    # fg all A, bg all A  -> not a CAAS

    layout = {
        "b_1": [NOHIT, HIT, NOHIT],
        "b_2": [HIT, NOHIT, HIT],
        "b_3": [NOHIT, NOHIT, NOHIT],
        "b_4": [NOHIT, HIT, NOHIT],
    }
    resample_dir = os.path.join(tmp, "nw_tree.resampled.output")
    os.makedirs(resample_dir)
    lines = []
    for b, hyps in layout.items():
        for m, (fg, bg) in enumerate(hyps, 1):
            lines.append(f"{b}~H{m}\t{fg}\t{bg}")
    with open(os.path.join(resample_dir, "fop_labelings.tab"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    # a decoy resample_*.tab must be IGNORED under --fop
    with open(os.path.join(resample_dir, "resample_000.tab"), "w") as fh:
        fh.write("b_1\ts1,s2\ts3,s4\n")

    out = os.path.join(tmp, "fop.out")
    boot.boot_on_single_alignment(
        trait_config_file="DUMMY",
        resampled_traits=resample_dir,           # directory -> redirected to fop_labelings.tab
        sliced_object=sliced,
        max_fg_gaps="NO", max_bg_gaps="NO", max_overall_gaps="NO",
        max_fg_miss="NO", max_bg_miss="NO", max_overall_miss="NO",
        the_admitted_patterns="1,2,3",
        output_file=out,
        max_conserved=0,
        caap_mode=False,
        fop_mode=True,
    )
    rows = [l.split("\t") for l in open(out).read().splitlines() if l.strip()]
    tag = "vectorized" if vectorize else "scalar"
    print(f"Part 2 ({tag}): end-to-end fop_mode")
    _check(f"[{tag}] one output row", len(rows) == 1, f"got {rows}")
    r = rows[0]
    _check(f"[{tag}] position/scheme", r[0] == "GENEX@10" and r[1] == "US", f"got {r[:2]}")
    _check(f"[{tag}] hits are base-cycle units (3, not 4)", r[2] == "3", f"got {r[2]}")
    _check(f"[{tag}] total is base cycles (4, not 12)", r[3] == "4", f"got {r[3]}")
    _check(f"[{tag}] proportion 0.75", abs(float(r[4]) - 0.75) < 1e-9, f"got {r[4]}")


def main():
    test_collapse_unit()
    test_end_to_end(vectorize=True)
    test_end_to_end(vectorize=False)
    print(f"\n{'ALL PASS' if FAILS == 0 else 'FAILURES'} — {FAILS} failure(s)")
    sys.exit(0 if FAILS == 0 else 1)


if __name__ == "__main__":
    main()
