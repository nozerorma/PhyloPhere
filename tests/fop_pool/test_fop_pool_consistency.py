#!/usr/bin/env python3
# test_fop_pool_consistency.py: the pooled score is strength (mean over the domains of the noisy-OR of the pair independences) times breadth (n_sp / S).
# PhyloPhere | tests/fop_pool/
"""Checks pool_domains: the worked example, species counted once, the mean (not the maximum) over the species
combinations of two domains, a missing S, and the observed / null parity of the reducer."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"))
from src.convergence.fop_pool import pool_domains  # noqa: E402
from src.core.driver import _expand_pooled  # noqa: E402


def rec(hyp, species, pairs, scores=None, domains=(1, 2, 3)):
    """One hypothesis record on the top side. `species` {domain: species carrying the derived residue},
    `pairs` [(a, b, c)] the independence of two converging domains; `scores` defaults to the noisy-OR of `pairs`."""
    if scores is None:
        scores = {d: 1 - math.prod(1 - c for a, b, c in pairs if d in (a, b)) for d in species
                  if any(d in (a, b) for a, b, c in pairs)}
    side = {"domain_scores": scores, "domain_der": {d: "A" for d in species}, "domain_der_enc": {d: "A" for d in species},
            "domain_anc": {d: "V" for d in species}, "domain_der_species": dict(species),
            "pair_lca": [[a, b, 99, c] for a, b, c in pairs]}
    return {"hyp": hyp, "sides": {"top": side, "bottom": {}, "domain_meta": {d: {} for d in domains}}}


def top(pooled):
    return pooled["top"]["asr_path_score"]


def main():
    # Worked example: K = 3, S = 20; A and B converge with c = 0.9, the residue sits in 2 species of A and 1 of B,
    # C does not change. s = (0.9, 0.9, 0): strength 0.6, breadth 3 / 20, score 0.09. The two species of A are two
    # combinations of the same pair of domains (same c), so their mean is 0.9.
    pair = [(1, 2, 0.9)]
    h1 = rec("H1", {1: "sA1", 2: "sB1"}, pair)
    h2 = rec("H2", {1: "sA2", 2: "sB1"}, pair)
    pooled = pool_domains([h1, h2], None, 20)
    assert math.isclose(top(pooled), 0.09), top(pooled)
    assert pooled["top"]["n_sp"] == 3 and pooled["n_hypotheses"] == 2

    # Without S the breadth is 1 (also S = 0): the score is the strength.
    assert math.isclose(top(pool_domains([h1, h2], None, None)), 0.6)
    assert math.isclose(top(pool_domains([h1, h2], None, 0)), 0.6)

    # Species are counted once whatever the number of hypotheses that use them; hypotheses repeating the same
    # combination add nothing to the pair either.
    same = pool_domains([rec(f"H{i}", {1: "sA", 2: "sB"}, pair) for i in (1, 2, 3)], None, 20)
    assert same["top"]["n_sp"] == 2 and math.isclose(top(same), 0.6 * 2 / 20)

    # A species of a domain with score 0 does not count.
    zero = pool_domains([rec("H1", {1: "sA", 2: "sB"}, pair, scores={1: 0.9, 2: 0.0})], None, 20)
    assert zero["top"]["n_sp"] == 1

    # Without pair_lca the score is 0 and the pooled record still carries M.
    none = pool_domains([rec("H1", {1: "sA", 2: "sB"}, [], scores={})], None, 20)
    assert top(none) == 0.0 and none["n_hypotheses"] == 1
    rows = _expand_pooled(7, "US", None, none)
    assert rows[0].side == "none" and rows[0].n_hypotheses == 1

    # Two combinations of the same two domains take the MEAN of their independences, not the maximum.
    lo, hi = rec("H1", {1: "sA1", 2: "sB1"}, [(1, 2, 0.5)]), rec("H2", {1: "sA2", 2: "sB1"}, [(1, 2, 0.9)])
    mixed = pool_domains([lo, hi], None, None)
    assert math.isclose(top(mixed), (0.7 + 0.7 + 0.0) / 3), top(mixed)

    # A domain with two partners takes the noisy-OR of their independences.
    three = pool_domains([rec("H1", {1: "a", 2: "b", 3: "c"}, [(1, 2, 0.5), (1, 3, 0.5), (2, 3, 0.0)])], None, None)
    assert math.isclose(top(three), ((1 - 0.25) + 0.5 + 0.5) / 3), top(three)

    # The score is clamped to [0, 1]: more species than S (a design that lacks one) caps the breadth at 1.
    assert top(pool_domains([h1, h2], None, 2)) <= 1.0

    # Observed and null parity: the same records, design keys and S give the same pooled result whatever the
    # hypothesis labels are called, and the design keys complete K with the domains no record holds.
    keys = {(h, d): 1.0 for h in ("H1", "H2") for d in (1, 2, 3, 4)}
    a = pool_domains([rec("H1", {1: "sA", 2: "sB"}, pair)], keys, 20)
    b = pool_domains([rec("H7", {1: "sA", 2: "sB"}, pair)], keys, 20)
    assert top(a) == top(b) and math.isclose(top(a), 0.45 * 2 / 20), top(a)    # K = 4: strength 1.8 / 4

    print("test_fop_pool_consistency: OK")


if __name__ == "__main__":
    main()
