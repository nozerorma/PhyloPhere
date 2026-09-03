#!/usr/bin/env python3
"""Tests for fop_pool.py — the permulation-null twin of fop_pool.R.

Run: python -m pytest test_fop_pool.py   (or: python test_fop_pool.py)
The FOP-case numbers are the SAME scenario as test_fop_pool.R so the two
implementations are checked against one shared oracle.
"""
import importlib.util
import sys
from pathlib import Path

# .../local (has src/) so fop_pool's `from src.biochem.grouping import ...` resolves
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

_spec = importlib.util.spec_from_file_location(
    "fop_pool", Path(__file__).with_name("fop_pool.py"))
_fp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fp)
p_at_least_2, pool_hypotheses, base_cycle = _fp.p_at_least_2, _fp.pool_hypotheses, _fp.base_cycle


def approx(a, b, tol=1e-6):
    return a is not None and b is not None and abs(a - b) < tol


def test_p_at_least_2():
    assert approx(p_at_least_2([1, 1]), 1.0)
    assert approx(p_at_least_2([0.5, 0.5]), 0.25)
    assert approx(p_at_least_2([0.9, 0.1]), 0.09)
    assert approx(p_at_least_2([0.8, 0.8]), 0.64)
    assert approx(p_at_least_2([0.8, 0.8, 0.8]), 0.896)
    assert p_at_least_2([1.0]) == 0.0


def test_base_cycle():
    assert base_cycle("b_5~H3") == "b_5"
    assert base_cycle("b_12") == "b_12"


def test_single_hypothesis_passthrough():
    rec = [{"hyp": "H1", "asr_path_score": 0.42, "core": 0.6,
            "independence": 0.9, "mrca_diversity": 0.5,
            "derived_agreement": 1.0, "conservation_gate": 1.0,
            "pair_scores": {1: 0.8, 2: 0.7}}]
    out = pool_hypotheses(rec, {})
    assert approx(out["asr_path_score"], 0.42)
    assert out["n_hypotheses"] == 1


def test_fop_two_hypotheses():
    # Same scenario as test_fop_pool.R's FOP case — two-job PSS weighting.
    # H1: strong domain-1 pair (s=0.9, own pss 10) but weak domain 2 (own pss 1).
    # H2: domain-1 pair s=0.2 (own pss 2). Domain 2 s=0.8 for both.
    #   Job A (per-domain c_d, pair-own PSS):
    #     c_1 = wmean(s=(0.9,0.2), w=(10,2)) = (9 + 0.4)/12 = 0.7833333
    #     c_2 = wmean(s=(0.8,0.8), w=(1,2))  = 0.8
    #   Job B (axis pools, per-hypothesis MEAN PSS):
    #     w(H1) = mean(10,1) = 5.5 ; w(H2) = mean(2,2) = 2
    #     diversity = wmean((1.0,0.5), w=(5.5,2)) = (5.5+1)/7.5 = 0.8666667
    recs = [
        {"hyp": "H1", "asr_path_score": 0.70, "independence": 1.0,
         "mrca_diversity": 1.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.72, "pair_scores": {1: 0.9, 2: 0.8}},
        {"hyp": "H2", "asr_path_score": 0.20, "independence": 1.0,
         "mrca_diversity": 0.5, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.24, "pair_scores": {1: 0.2, 2: 0.8}},
    ]
    pss = {("H1", 1): 10.0, ("H1", 2): 1.0, ("H2", 1): 2.0, ("H2", 2): 2.0}
    c1, c2, div = 9.4 / 12, 0.8, 6.5 / 7.5
    out = pool_hypotheses(recs, pss)
    assert out["n_hypotheses"] == 2
    assert approx(out["core"], p_at_least_2([c1, c2])), out["core"]
    assert approx(out["mrca_diversity"], div), out["mrca_diversity"]
    exp = 1.0 * out["core"] * (0.75 + 0.25 * div) * 1.0 * 1.0
    assert approx(out["asr_path_score"], exp), (out["asr_path_score"], exp)
    assert 0.20 < out["asr_path_score"] < 0.70
    # Job A pair-PSS-driven: the s=0.9 pair (own pss 10) dominates c_1 even though
    # its hypothesis H1 has a weak sibling domain (pss 1); a hyp-min weight -> ~0.43.
    assert p_at_least_2([c1, c2]) and c1 > 0.75
    # Job B mean-PSS-driven: H1 (mean 5.5) outweighs H2 (2), so pooled diversity
    # leans to H1's 1.0, above the min-weighted 0.667.
    assert out["mrca_diversity"] > 0.8


def test_fop_equal_weights_no_pss():
    recs = [
        {"hyp": "H1", "asr_path_score": 0.70, "independence": 1.0,
         "mrca_diversity": 1.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.72, "pair_scores": {1: 0.9, 2: 0.8}},
        {"hyp": "H2", "asr_path_score": 0.20, "independence": 1.0,
         "mrca_diversity": 0.5, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.24, "pair_scores": {1: 0.2, 2: 0.8}},
    ]
    out = pool_hypotheses(recs, {})   # no weights -> equal
    # c_1 = mean(0.9, 0.2) = 0.55 ; c_2 = 0.8
    assert approx(out["core"], p_at_least_2([0.55, 0.8]))
    assert approx(out["mrca_diversity"], 0.75)


def test_point2_conserved_gate_dedup():
    # Same scenario as test_fop_pool.R's POINT 2 case. H1 and H2 SHARE one
    # conserved pair (same node "cn_shared", cons 0.6); H3 has a unique conserved
    # pair ("cn_uniq", cons 0.4). Rebuilt gate = 0.5 + 0.5*mean(0.6,0.4) = 0.75,
    # counting the shared pair ONCE — NOT the 3-way per-hypothesis gate mean
    # (0.8, 0.8, 0.7) = 0.76667.
    recs = [
        {"hyp": "H1", "asr_path_score": 0.5, "independence": 1.0,
         "mrca_diversity": 0.0, "derived_agreement": 1.0, "conservation_gate": 0.8,
         "core": 0.5, "pair_scores": {1: 0.7, 2: 0.6},
         "conserved_pair_scores": {5: 0.6}, "conserved_pair_nodes": {5: "cn_shared"}},
        {"hyp": "H2", "asr_path_score": 0.5, "independence": 1.0,
         "mrca_diversity": 0.0, "derived_agreement": 1.0, "conservation_gate": 0.8,
         "core": 0.5, "pair_scores": {1: 0.7, 2: 0.6},
         "conserved_pair_scores": {5: 0.6}, "conserved_pair_nodes": {5: "cn_shared"}},
        {"hyp": "H3", "asr_path_score": 0.4, "independence": 1.0,
         "mrca_diversity": 0.0, "derived_agreement": 1.0, "conservation_gate": 0.7,
         "core": 0.4, "pair_scores": {1: 0.6, 2: 0.5},
         "conserved_pair_scores": {5: 0.4}, "conserved_pair_nodes": {5: "cn_uniq"}},
    ]
    out = pool_hypotheses(recs, {})
    assert approx(out["conservation_gate"], 0.75), out["conservation_gate"]
    assert not approx(out["conservation_gate"], (0.8 + 0.8 + 0.7) / 3)

    # No conserved-pair identity -> fallback to pooling per-hypothesis gates.
    plain = [{k: v for k, v in r.items()
              if k not in ("conserved_pair_scores", "conserved_pair_nodes")}
             for r in recs]
    out2 = pool_hypotheses(plain, {})
    assert approx(out2["conservation_gate"], (0.8 + 0.8 + 0.7) / 3), out2["conservation_gate"]

    # Identity present but no conserved pairs -> novel position, gate 1.0.
    empty = [{**r, "conserved_pair_scores": {}, "conserved_pair_nodes": {}}
             for r in recs]
    out3 = pool_hypotheses(empty, {})
    assert approx(out3["conservation_gate"], 1.0), out3["conservation_gate"]


def test_point3_harvest_wide_da():
    # Twin of test_fop_pool.R's POINT 3 case. Two null hypotheses, K=2 domains,
    # bottom-side changed pairs landing on raw {I, V} in domain 1 and {I, V} in
    # domain 2 (a between-hypothesis split). No node identity on the null, so the
    # dedup is by (domain, side, raw_aa): 4 distinct changed residues remain.
    #   US : bot I,V,I,V -> plurality 2/4 = 0.5
    #   GS4: I,V both 'h' -> 4/4 = 1.0
    recs = [
        {"hyp": "H1", "asr_path_score": 0.5, "independence": 1.0,
         "mrca_diversity": 0.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.5, "pair_scores": {1: 0.7, 2: 0.6},
         "pair_derived_top": {}, "pair_derived_bot": {1: "I", 2: "V"}},
        {"hyp": "H2", "asr_path_score": 0.4, "independence": 1.0,
         "mrca_diversity": 0.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.4, "pair_scores": {1: 0.7, 2: 0.6},
         "pair_derived_top": {}, "pair_derived_bot": {1: "V", 2: "I"}},
    ]
    out_us = pool_hypotheses(recs, {}, scheme="US")
    assert approx(out_us["derived_agreement"], 0.5), out_us["derived_agreement"]
    out_gs4 = pool_hypotheses(recs, {}, scheme="GS4")
    assert approx(out_gs4["derived_agreement"], 1.0), out_gs4["derived_agreement"]

    # No raw residues -> fall back to per-hypothesis da wmean (all 1.0).
    plain = [{k: v for k, v in r.items()
              if k not in ("pair_derived_top", "pair_derived_bot")} for r in recs]
    assert approx(pool_hypotheses(plain, {}, scheme="US")["derived_agreement"], 1.0)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    ok = True
    for fn in fns:
        try:
            fn()
            print("PASS ", fn.__name__)
        except AssertionError as e:
            ok = False
            print("FAIL ", fn.__name__, "--", e)
    print("\n" + ("ALL TESTS PASSED" if ok else "SOME TESTS FAILED"))
    sys.exit(0 if ok else 1)
