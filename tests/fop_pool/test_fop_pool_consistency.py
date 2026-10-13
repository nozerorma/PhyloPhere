#!/usr/bin/env python3
# test_fop_pool_consistency.py: the pooled score of the observed labeling and of a null cycle goes through one rule.
# PhyloPhere | tests/fop_pool/
"""Checks pool_domains: the support shrinkage M / (M + 2) is a property of the design (the number H of hypotheses of
the PSS map), not of whether a map happens to be present or empty."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"))
from src.convergence.fop_pool import pool_domains, effective_number  # noqa: E402
from src.core.driver import _expand_pooled  # noqa: E402


def record(hyp, scores):
    """One hypothesis record: domain scores on the top side, three domains in the design."""
    side = {"domain_scores": scores, "domain_der": {d: "A" for d in scores}, "domain_der_enc": {d: "A" for d in scores},
            "domain_anc": {d: "V" for d in scores}}
    return {"hyp": hyp, "sides": {"top": side, "bottom": {}, "domain_meta": {1: {}, 2: {}, 3: {}}}}


def pss(hyps):
    return {(h, d): 2.0 for h in hyps for d in (1, 2, 3)}


def top_score(pooled):
    return pooled["top"]["asr_path_score"]


def main():
    one = [record("H1", {1: 1.0, 2: 1.0})]
    base = top_score(pool_domains(one, None))                       # no PSS: weights 1, no shrinkage
    assert math.isclose(base, 2 / 3), base

    # Single-hypothesis design with a PSS map: not shrunk, equal weights (constant PSS) -> same as no map.
    single = top_score(pool_domains(one, pss(["H1"])))
    assert math.isclose(single, base), (single, base)

    # Two-hypothesis design: one detecting hypothesis (M = 1) is shrunk by 1 / 3.
    multi = top_score(pool_domains(one, pss(["H1", "H2"])))
    assert math.isclose(multi, base / 3), (multi, base)

    # M = 2 of the same design: both detect, factor 2 / 4.
    two = [record("H1", {1: 1.0, 2: 1.0}), record("H2", {1: 1.0, 2: 1.0})]
    both = top_score(pool_domains(two, pss(["H1", "H2"])))
    assert math.isclose(both, base / 2), (both, base)

    # An empty map and None are the same input.
    assert top_score(pool_domains(one, {})) == top_score(pool_domains(one, None))

    # Observed and null parity: the same records and the same map give the same pooled result, whatever the
    # hypothesis labels of the cycle are called, because only the map decides the design size.
    obs = pool_domains(one, pss(["H1", "H2", "H3"]))
    nul = pool_domains([record("H1", {1: 1.0, 2: 1.0})], pss(["H1", "H2", "H3"]))
    assert obs == nul
    # The pooled record of a null position carries M, the number of hypotheses that detected it (also for a
    # position no side of which converges), so the null and the observed support can be compared.
    two_pooled = pool_domains(two, pss(["H1", "H2", "H3"]))
    rows = _expand_pooled(7, "US", None, two_pooled)
    assert rows and all(r.n_hypotheses == 2 for r in rows), rows
    none = _expand_pooled(7, "US", None, pool_domains([record("H1", {})], pss(["H1", "H2", "H3"])))
    assert none[0].side == "none" and none[0].n_hypotheses == 1
    # Effective number of independent hypotheses: identical supports count once, disjoint ones in full.
    a, b, c = frozenset({("1", 10), ("2", 20)}), frozenset({("1", 11), ("2", 21)}), frozenset({("1", 10), ("2", 21)})
    assert effective_number([]) == 0.0 and effective_number([a]) == 1.0
    assert effective_number([a, a, a, a]) == 1.0
    assert effective_number([a, b]) == 2.0
    assert 1.0 < effective_number([a, a, b]) < 3.0
    assert 1.0 < effective_number([a, b, c]) < 3.0          # c shares one pair with each of a and b
    # pool_domains reports it per side: hypotheses with the same supporting pairs are one effective hypothesis.
    def rec_pair(hyp, node):
        side = {"domain_scores": {1: 1.0, 2: 1.0}, "domain_der": {1: "A", 2: "A"}, "domain_der_enc": {1: "A", 2: "A"},
                "domain_anc": {1: "V", 2: "V"}}
        return {"hyp": hyp, "sides": {"top": side, "bottom": {}, "domain_meta": {1: {"mrca_id": node}, 2: {"mrca_id": node + 1}, 3: {}}}}
    same = pool_domains([rec_pair("H1", 5), rec_pair("H2", 5), rec_pair("H3", 5)], pss(["H1", "H2", "H3"]))
    diff = pool_domains([rec_pair("H1", 5), rec_pair("H2", 7), rec_pair("H3", 9)], pss(["H1", "H2", "H3"]))
    assert same["top"]["n_eff"] == 1.0 and diff["top"]["n_eff"] == 3.0, (same["top"]["n_eff"], diff["top"]["n_eff"])
    assert same["n_hypotheses"] == diff["n_hypotheses"] == 3
    print("test_fop_pool_consistency: OK")


if __name__ == "__main__":
    main()
