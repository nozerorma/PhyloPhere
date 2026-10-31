'''
Positional universe of the observed design.

MODULE NAME: universe.py
DESCRIPTION: positions that can score under the observed hypotheses H1..Hn, with the CT kernel itself
DEPENDENCIES: perm_replay.py, perm_replay_io.py
CALLED BY: perm_replay.py (export_b0_universe), tests

A position belongs to the universe of a gene when
  1. b_0 examined it (the `background_sink` rule of `VectorizedPermReplay.count`: some b_0 hypothesis passes the gap,
     missing-species and miss_pair limits and has a residue on both sides), and
  2. some orientation of the pairs of some observed hypothesis passes the whole CT test: the same limits, the same
     patterns, `max_conserved` and `min_non_overlapping`.

An orientation swaps the two members of some pairs (foreground <-> background inside a pair). No pair is created,
removed or re-matched, so the hypotheses keep the pairs the contrast selection fixed; the Dunn index of a hypothesis
depends on the pair partition and not on which member is foreground, so it needs no re-evaluation (the equality is
checked against lean_contrast_selector.R in tests/test_positional_universe.py).

The orientations are replayed by `VectorizedPermReplay`, the kernel of the permulations and of b_0, so there is no
second implementation of the CAAS rule. Every orientation is evaluated (2**K per hypothesis): the per-side gap and
missing limits and miss_pair are not symmetric under a swap, so no orientation can be dropped as redundant in general.
The observed orientation is one of them, which makes every observed detection part of the universe.
'''

from modules.perm_replay_io import labelings_to_multicfg

# Safety bound on the pairs per hypothesis (2**K orientations are replayed); not a scientific parameter.
DEFAULT_MAX_PAIRS = 12


def observed_labelings(cfg):
    """{tag: (fg, bg)} of the b_0 labelings (every hypothesis of the observed design) of a labelings multicfg."""
    from modules.perm_replay import _fop_base_cycle
    labelings = {t: (tuple(cfg.trait2fg[t]), tuple(cfg.trait2bg[t])) for t in cfg.alltraits
                 if _fop_base_cycle(t) == "b_0"}
    if not labelings:
        raise ValueError("the positional universe needs the b_0 labelings of the observed design")
    return labelings


def orientation_labelings(labelings, max_pairs=DEFAULT_MAX_PAIRS):
    """Every orientation of the pairs of every labeling: {"<tag>#<mask>": (fg, bg)}.

    Pair k is fg[k] with bg[k]; bit k of `mask` swaps its members. mask 0 is the labeling itself."""
    out = {}
    for tag, (fg, bg) in labelings.items():
        k = len(fg)
        if len(bg) != k:
            raise ValueError(f"unpaired labeling {tag}")
        if k > max_pairs:
            raise ValueError(f"{tag} has {k} pairs: its 2**{k} orientations exceed the bound of {max_pairs} pairs")
        for mask in range(1 << k):
            swap = [(mask >> j) & 1 for j in range(k)]
            out[f"{tag}#{mask}"] = (tuple(bg[j] if swap[j] else fg[j] for j in range(k)),
                                    tuple(fg[j] if swap[j] else bg[j] for j in range(k)))
    return out


def eligible_positions(cfg, sliced_object, positions, genename, thresholds, max_pairs=DEFAULT_MAX_PAIRS):
    """Position numbers (strings, as in background.output) that pass the CT test under some orientation.

    `cfg` is the labelings multicfg of the b_0 hypotheses, `thresholds` the keyword arguments `count` takes for the
    CT limits: maxgaps_fg/bg/all, maxmiss_fg/bg/all, max_conserved, admitted_patterns, miss_pair, min_non_overlapping."""
    from modules.perm_replay import VectorizedPermReplay
    oriented = labelings_to_multicfg(orientation_labelings(observed_labelings(cfg), max_pairs))
    counts = VectorizedPermReplay(oriented, sliced_object.species).count(positions, genename, **thresholds)
    return {name.split("@", 1)[1] for name, n in counts.items() if n > 0}
