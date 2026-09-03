#!/usr/bin/env python3
"""Tests for fop_pool.py — the permulation-null twin of fop_pool.R.

Run: python -m pytest test_fop_pool.py   (or: python test_fop_pool.py)
The FOP-case numbers are the SAME scenario as test_fop_pool.R so the two
implementations are checked against one shared oracle.
"""
import importlib.util
import sys
from pathlib import Path

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
    # Same scenario as test_fop_pool.R's FOP case.
    # Domain 1: H1 s=0.9 (w=8), H2 s=0.3 (w=2)  -> c_1 = (7.2+0.6)/10 = 0.78
    #   wait: weights are min-pair-PSS per hypothesis. H1 min(9,8)=8, H2 min(2,8)=2.
    #   c_1 = (0.9*8 + 0.3*2) / (8+2) = 7.8/10 = 0.78
    # Domain 2: both s=0.8 -> c_2 = 0.8
    # core = p_at_least_2(0.78, 0.8) = 0.624
    recs = [
        {"hyp": "H1", "asr_path_score": 0.70, "independence": 1.0,
         "mrca_diversity": 1.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.72, "pair_scores": {1: 0.9, 2: 0.8}},
        {"hyp": "H2", "asr_path_score": 0.20, "independence": 1.0,
         "mrca_diversity": 0.5, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.24, "pair_scores": {1: 0.3, 2: 0.8}},
    ]
    pss = {("H1", 1): 9.0, ("H1", 2): 8.0, ("H2", 1): 2.0, ("H2", 2): 8.0}
    out = pool_hypotheses(recs, pss)
    assert out["n_hypotheses"] == 2
    assert approx(out["core"], p_at_least_2([0.78, 0.8])), out["core"]
    # diversity pooled = wmean(1.0, 0.5; w=8,2) = (8+1)/10 = 0.9
    assert approx(out["mrca_diversity"], 0.9), out["mrca_diversity"]
    exp = 1.0 * out["core"] * (0.75 + 0.25 * 0.9) * 1.0 * 1.0
    assert approx(out["asr_path_score"], exp), (out["asr_path_score"], exp)
    # sits between the two hypotheses, not at the max (weighted-mean semantics)
    assert 0.20 < out["asr_path_score"] < 0.70


def test_fop_equal_weights_no_pss():
    recs = [
        {"hyp": "H1", "asr_path_score": 0.70, "independence": 1.0,
         "mrca_diversity": 1.0, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.72, "pair_scores": {1: 0.9, 2: 0.8}},
        {"hyp": "H2", "asr_path_score": 0.20, "independence": 1.0,
         "mrca_diversity": 0.5, "derived_agreement": 1.0, "conservation_gate": 1.0,
         "core": 0.24, "pair_scores": {1: 0.3, 2: 0.8}},
    ]
    out = pool_hypotheses(recs, {})   # no weights -> equal
    # c_1 = mean(0.9, 0.3) = 0.6 ; c_2 = 0.8
    assert approx(out["core"], p_at_least_2([0.6, 0.8]))


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
