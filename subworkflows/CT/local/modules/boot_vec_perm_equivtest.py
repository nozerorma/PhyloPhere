#!/usr/bin/env python3
# =============================================================================
# boot_vec_perm_equivtest.py — perm_discovery row parity (BOOTSTRAP_PERMS path)
# =============================================================================
# Proves the vectorized perm_discovery emitter reproduces the scalar caasboot
# perm_discovery output EXACTLY (as a row multiset; row ORDER legitimately
# differs because scalar iterates hash-ordered set intersections, and the
# downstream concat→disambiguate treats rows independently). Also checks the
# count output (.bootstraped.output equivalent) stays byte-identical.
#
# Covers: classical + CAAP, max_conserved 0 and >0, gaps/ambiguity/missing,
# directory mode (what BOOTSTRAP_PERMS uses) and single-file mode.
#
# Run:  micromamba run -n phylophere python \
#         subworkflows/CT/local/modules/boot_vec_perm_equivtest.py
# =============================================================================

import os
import sys
import random
import tempfile
import importlib

_HERE = os.path.dirname(os.path.abspath(__file__))
_LOCAL = os.path.dirname(_HERE)
if _LOCAL not in sys.path:
    sys.path.insert(0, _LOCAL)

STD = "ACDEFGHIKLMNPQRSTVWY"
FAILS = 0


class Sliced:
    def __init__(self, genename, species, d):
        self.genename = genename
        self.species = species
        self.d = d


def make_positions(rng, species, n_pos):
    d = []
    for p in range(n_pos):
        pos = {}
        for sp in species:
            r = rng.random()
            aa = "-" if r < 0.12 else (rng.choice("XBZJU") if r < 0.17 else rng.choice(STD))
            pos[sp] = f"{aa}@{p}"
        d.append(pos)
    return d


def write_resample_dir(rng, allsp, n_files, cyc, path):
    os.makedirs(path, exist_ok=True)
    for f in range(n_files):
        lines = []
        for c in range(cyc):
            pool = list(allsp); rng.shuffle(pool)
            lines.append("\t".join([f"b_{f}_{c}", ",".join(pool[:6]), ",".join(pool[6:12])]))
        with open(os.path.join(path, f"resample_{f:03d}.tab"), "w") as fh:
            fh.write("\n".join(lines) + "\n")


def write_resample_file(rng, allsp, cyc, path):
    lines = []
    for c in range(cyc):
        pool = list(allsp); rng.shuffle(pool)
        lines.append("\t".join([f"b_{c}", ",".join(pool[:6]), ",".join(pool[6:12])]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def run(vectorize, sliced, resampled, out, perm_out, caap, max_conserved, single):
    os.environ["CT_BOOTSTRAP_VECTORIZE"] = "1" if vectorize else "0"
    import modules.boot_vec as bv; importlib.reload(bv)
    import modules.boot as boot; importlib.reload(boot)
    from modules.init_bootstrap import simtrait_revive
    res = simtrait_revive(resampled) if single else resampled
    boot.boot_on_single_alignment(
        trait_config_file="DUMMY",
        resampled_traits=res,
        sliced_object=sliced,
        max_fg_gaps="NO", max_bg_gaps="NO", max_overall_gaps="3",
        max_fg_miss="NO", max_bg_miss="NO", max_overall_miss="3",
        the_admitted_patterns="1,2,3",
        output_file=out,
        max_conserved=max_conserved,
        caap_mode=caap,
        export_perm_discovery=perm_out,
    )


def rows_multiset(path):
    """Return (header, sorted body lines) for a perm_discovery file."""
    if not os.path.exists(path):
        return None, []
    lines = open(path).read().splitlines()
    if not lines:
        return None, []
    return lines[0], sorted(lines[1:])


def check(label, **kw):
    global FAILS
    tmp = kw.pop("tmp")
    caap = kw["caap"]; mc = kw["max_conserved"]; single = kw["single"]
    species = [f"sp{i}" for i in range(14)]
    extra = [f"ex{i}" for i in range(4)]
    allsp = species + extra
    rng = random.Random(hash(label) & 0xffff)
    sliced = Sliced("GENEX", species, make_positions(rng, species, 26))

    tag = label.replace(" ", "_")
    if single:
        res = os.path.join(tmp, f"{tag}.tab"); write_resample_file(random.Random(5), allsp, 500, res)
    else:
        res = os.path.join(tmp, f"{tag}_dir"); write_resample_dir(random.Random(5), allsp, 3, 200, res)

    ov, op_v = os.path.join(tmp, f"{tag}_v.out"), os.path.join(tmp, f"{tag}_v.disc")
    os_, op_s = os.path.join(tmp, f"{tag}_s.out"), os.path.join(tmp, f"{tag}_s.disc")
    run(True, sliced, res, ov, op_v, caap, mc, single)
    run(False, sliced, res, os_, op_s, caap, mc, single)

    # 1. count output byte-identical
    cv = open(ov).read() if os.path.exists(ov) else "<none>"
    cs = open(os_).read() if os.path.exists(os_) else "<none>"
    count_ok = cv == cs
    # 2. perm_discovery: header identical + body multiset identical
    hv, bv_ = rows_multiset(op_v)
    hs, bs_ = rows_multiset(op_s)
    head_ok = (hv == hs)
    body_ok = (bv_ == bs_)
    ok = count_ok and head_ok and body_ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:<34} rows(s/v)={len(bs_)}/{len(bv_)} "
          f"count={'ok' if count_ok else 'DIFF'} header={'ok' if head_ok else 'DIFF'} "
          f"body={'ok' if body_ok else 'DIFF'}")
    if not ok:
        FAILS += 1
        if not body_ok:
            sset, vset = set(bs_), set(bv_)
            only_s = list(sset - vset)[:2]
            only_v = list(vset - sset)[:2]
            for r in only_s: print(f"        only in scalar: {r}")
            for r in only_v: print(f"        only in vector: {r}")


def main():
    tmp = tempfile.mkdtemp(prefix="boot_vec_perm_")
    print(f"boot_vec_perm_equivtest — scratch={tmp}\n")
    for single in (False, True):
        mode = "single" if single else "dir"
        check(f"classical mc0 {mode}", caap=False, max_conserved=0, single=single, tmp=tmp)
        check(f"classical mc2 {mode}", caap=False, max_conserved=2, single=single, tmp=tmp)
        check(f"caap mc0 {mode}",      caap=True,  max_conserved=0, single=single, tmp=tmp)
        check(f"caap mc2 {mode}",      caap=True,  max_conserved=2, single=single, tmp=tmp)
    print(f"\n{'ALL PASS' if FAILS == 0 else 'FAILURES'} — {FAILS} failure(s)")
    sys.exit(0 if FAILS == 0 else 1)


if __name__ == "__main__":
    main()
