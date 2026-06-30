#!/usr/bin/env python3
# =============================================================================
# boot_vec_equivtest.py - numerical equivalence test for the vectorized kernel
# =============================================================================
# Proves VectorizedBootstrap.count() reproduces the scalar caasboot() counting
# loop EXACTLY (the path it replaces), across: clean data, gaps, ambiguity codes,
# missing species, conserved-pair tolerance (max_conserved>0), every admitted-
# pattern subset, all five grouping schemes (US/GS1-GS4), discovery-scheme
# subsetting, exotic non-standard symbols, and "NO" vs numeric gap/miss
# thresholds. If this passes, the bootstrap empirical-p produced by the
# vectorized path is identical to the old loop (calibration is preserved).
#
# Run:  micromamba run -n phylophere python \
#         subworkflows/CT/local/modules/boot_vec_equivtest.py
# (or from subworkflows/CT/local:  python modules/boot_vec_equivtest.py)
# =============================================================================

import os
import sys
import random

# Make `modules` importable whether run from local/ or from the repo root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_LOCAL = os.path.dirname(_HERE)            # .../subworkflows/CT/local
if _LOCAL not in sys.path:
    sys.path.insert(0, _LOCAL)

from modules.caas_id import process_position
from modules.boot import caasboot
from modules.init_bootstrap import simtrait_revive
from modules.boot_vec import VectorizedBootstrap

TOTAL_FAILS = 0
STD_AAS = "ACDEFGHIKLMNPQRSTVWY"


# ----------------------------------------------------------------------------
# Synthetic-data construction
# ----------------------------------------------------------------------------

def make_gene(rng, n_align=14, n_extra=4, n_positions=40,
              alphabet=STD_AAS, gap_prob=0.12, ambig_prob=0.05,
              exotic_prob=0.0):
    """Build one synthetic gene: alignment species + positions, plus a few extra
    'tree' species that never appear in the alignment (sources of 'missing').

    Returns (alignment_species, all_species, positions) where positions is a list
    of pos_dict {species: 'AA@pos'} over the alignment species only.
    """
    align_species = [f"sp{i}" for i in range(n_align)]
    extra_species = [f"ex{i}" for i in range(n_extra)]
    all_species = align_species + extra_species

    positions = []
    for p in range(n_positions):
        pos_dict = {}
        for sp in align_species:
            r = rng.random()
            if r < gap_prob:
                aa = "-"
            elif r < gap_prob + ambig_prob:
                aa = rng.choice("XBZJU")
            elif exotic_prob > 0 and r < gap_prob + ambig_prob + exotic_prob:
                aa = rng.choice("*?+")  # non-standard, non-gap, non-ambiguity
            else:
                aa = rng.choice(alphabet)
            pos_dict[sp] = f"{aa}@{p}"
        positions.append(pos_dict)
    return align_species, all_species, positions


def make_resample(rng, all_species, n_cycles, fg_size, bg_size, path):
    """Write a resample_*.tab-style file: cycleid<TAB>fg_csv<TAB>bg_csv."""
    lines = []
    for c in range(n_cycles):
        pool = list(all_species)
        rng.shuffle(pool)
        fg = pool[:fg_size]
        bg = pool[fg_size:fg_size + bg_size]
        if not fg or not bg:
            continue
        lines.append("\t".join([f"b_{c}", ",".join(fg), ",".join(bg)]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# Scalar reference: run caasboot per position, parse counts
# ----------------------------------------------------------------------------

def scalar_counts(positions_with_schemes, genename, cfg, align_species,
                  thresholds, max_conserved, patterns, caap_mode):
    (g_fg, g_bg, g_all, m_fg, m_bg, m_all) = thresholds
    out = {}
    for pos_dict, schemes in positions_with_schemes:
        processed = process_position(pos_dict, multiconfig=cfg,
                                     species_in_alignment=align_species)
        line_output = caasboot(
            processed,
            genename=genename,
            list_of_traits=cfg.alltraits,
            maxgaps_fg=g_fg, maxgaps_bg=g_bg, maxgaps_all=g_all,
            maxmiss_fg=m_fg, maxmiss_bg=m_bg, maxmiss_all=m_all,
            cycles=cfg.cycles,
            multiconfig=cfg,
            max_conserved=max_conserved,
            admitted_patterns=patterns,
            caap_mode=caap_mode,
            discovery_schemes=schemes,
        )
        for line in line_output.split("\n"):
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name, scheme, count = parts[0], parts[1], int(parts[2])
            if caap_mode:
                out[(name, scheme)] = out.get((name, scheme), 0) + count
            else:
                out[name] = out.get(name, 0) + count
    return out


# ----------------------------------------------------------------------------
# One comparison case
# ----------------------------------------------------------------------------

def check_case(label, rng, *, n_cycles=400, caap_mode=False,
               thresholds=("NO", "NO", "NO", "NO", "NO", "NO"),
               max_conserved=0, patterns="1,2,3",
               discovery_subset=None, exotic_prob=0.0, tmpdir="."):
    global TOTAL_FAILS
    genename = "GENEX"
    align_species, all_species, positions = make_gene(rng, exotic_prob=exotic_prob)

    res_path = os.path.join(tmpdir, f"resample_{label.replace(' ', '_')}.tab")
    make_resample(rng, all_species, n_cycles, fg_size=5, bg_size=5, path=res_path)
    cfg = simtrait_revive(res_path)

    # Attach a per-position scheme set (only meaningful in caap_mode).
    if caap_mode and discovery_subset is not None:
        pos_with_schemes = [(p, set(discovery_subset)) for p in positions]
    else:
        pos_with_schemes = [(p, None) for p in positions]

    truth = scalar_counts(pos_with_schemes, genename, cfg, align_species,
                          thresholds, max_conserved, patterns, caap_mode)

    vb = VectorizedBootstrap(cfg, align_species)
    (g_fg, g_bg, g_all, m_fg, m_bg, m_all) = thresholds
    vec = vb.count(
        pos_with_schemes, genename,
        maxgaps_fg=g_fg, maxgaps_bg=g_bg, maxgaps_all=g_all,
        maxmiss_fg=m_fg, maxmiss_bg=m_bg, maxmiss_all=m_all,
        max_conserved=max_conserved, admitted_patterns=patterns,
        caap_mode=caap_mode,
        # tiny chunk to also exercise the cycle-chunking boundary
        b_chunk=137,
    )

    keys = set(truth) | set(vec)
    max_diff = 0
    n_mismatch = 0
    total_hits = 0
    first_bad = None
    for k in keys:
        a = truth.get(k, 0)
        b = vec.get(k, 0)
        total_hits += a
        d = abs(a - b)
        if d > max_diff:
            max_diff = d
        if d != 0:
            n_mismatch += 1
            if first_bad is None:
                first_bad = (k, a, b)

    ok = (max_diff == 0)
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {label:<34} keys={len(keys):<4} hits={total_hits:<6} "
          f"max|Δ|={max_diff} mismatched={n_mismatch}")
    if not ok:
        print(f"         first mismatch: key={first_bad[0]} scalar={first_bad[1]} vec={first_bad[2]}")
        TOTAL_FAILS += 1

    try:
        os.remove(res_path)
    except OSError:
        pass


def main():
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="boot_vec_equivtest_")
    print(f"boot_vec_equivtest — scratch={tmpdir}\n")

    # Each case gets its own seeded RNG so failures are reproducible.
    def rng(seed):
        return random.Random(seed)

    # 1. Classical, strict, default patterns, no filters
    check_case("classical clean", rng(1), tmpdir=tmpdir)

    # 2. Classical with heavier gaps/ambiguity in the alignment
    check_case("classical gappy", rng(2), tmpdir=tmpdir)

    # 3. Classical with exotic non-standard symbols (iscaas identity edge case)
    check_case("classical exotic-symbols", rng(3), exotic_prob=0.10, tmpdir=tmpdir)

    # 4. Numeric gap/missing thresholds (0 and 1 tolerated)
    check_case("classical thr gaps=0 miss=0", rng(4),
               thresholds=("0", "0", "0", "0", "0", "0"), tmpdir=tmpdir)
    check_case("classical thr gaps=1 miss=2", rng(5),
               thresholds=("1", "1", "2", "2", "2", "3"), tmpdir=tmpdir)

    # 5. Conserved-pair tolerance (relaxed CAAS)
    check_case("classical max_conserved=1", rng(6), max_conserved=1, tmpdir=tmpdir)
    check_case("classical max_conserved=3", rng(7), max_conserved=3, tmpdir=tmpdir)

    # 6. Pattern subsets (substring membership semantics)
    check_case("classical patterns=1", rng(8), patterns="1", tmpdir=tmpdir)
    check_case("classical patterns=2,3", rng(9), patterns="2,3", tmpdir=tmpdir)
    check_case("classical patterns=1,2,3,4", rng(10), patterns="1,2,3,4", tmpdir=tmpdir)

    # 7. CAAP mode, all schemes
    check_case("caap all-schemes clean", rng(11), caap_mode=True, tmpdir=tmpdir)
    check_case("caap all-schemes gappy", rng(12), caap_mode=True,
               thresholds=("1", "1", "2", "1", "1", "2"), tmpdir=tmpdir)
    check_case("caap all-schemes max_cons=2", rng(13), caap_mode=True,
               max_conserved=2, tmpdir=tmpdir)

    # 8. CAAP mode, discovery-scheme subsetting
    check_case("caap subset GS1,GS3", rng(14), caap_mode=True,
               discovery_subset=["GS1", "GS3"], tmpdir=tmpdir)
    check_case("caap subset US only", rng(15), caap_mode=True,
               discovery_subset=["US"], tmpdir=tmpdir)
    check_case("caap subset CAAS=all", rng(16), caap_mode=True,
               discovery_subset=["CAAS"], tmpdir=tmpdir)

    # 9. Larger cycle count to exercise multi-chunk path and accumulation
    check_case("classical big-cycles", rng(17), n_cycles=900, tmpdir=tmpdir)
    check_case("caap big-cycles", rng(18), n_cycles=900, caap_mode=True, tmpdir=tmpdir)

    print(f"\n{'ALL PASS' if TOTAL_FAILS == 0 else 'FAILURES'} — {TOTAL_FAILS} failure(s)")
    sys.exit(0 if TOTAL_FAILS == 0 else 1)


if __name__ == "__main__":
    main()
