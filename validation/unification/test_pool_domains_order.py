"""pool_domains gives the same bits whatever order the hypotheses arrive in."""
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.convergence.fop_pool import pool_domains  # noqa: E402


def _records(rng, m, k=4):
    recs, pss = [], {}
    for i in range(1, m + 1):
        h = f"H{i}"
        side = lambda: {"domain_scores": {d: rng.random() for d in range(1, k + 1)}, "domain_der_enc": {}, "domain_der": {}, "domain_anc": {}}
        recs.append({"hyp": h, "sides": {"top": side(), "bottom": side(), "domain_meta": {d: {} for d in range(1, k + 1)}}})
        for d in range(1, k + 1):
            pss[(h, d)] = rng.random()
    return recs, pss


@pytest.mark.parametrize("m", [3, 12, 100])
def test_pooled_scores_do_not_depend_on_the_order_of_the_hypotheses(m):
    rng = random.Random(m)
    for _ in range(30):
        recs, pss = _records(rng, m)
        ref = pool_domains(recs, pss)
        shuffled = recs[:]
        rng.shuffle(shuffled)
        got = pool_domains(shuffled, pss)
        for side in ("top", "bottom"):
            assert got[side]["asr_path_score"] == ref[side]["asr_path_score"]
            assert got[side]["domain_scores"] == ref[side]["domain_scores"]
            assert got[side]["domain_weights"] == ref[side]["domain_weights"]


def test_two_hypotheses_and_the_unweighted_mean_are_unchanged():
    rng = random.Random(1)
    recs, _ = _records(rng, 2)
    out = pool_domains(recs, None)
    top = [r["sides"]["top"]["domain_scores"] for r in recs]
    for d in (1, 2, 3, 4):
        assert out["top"]["domain_scores"][d] == pytest.approx((top[0][d] + top[1][d]) / 2, abs=1e-15)
    assert out["top"]["asr_path_score"] == pytest.approx(sum(out["top"]["domain_scores"].values()) / 4, abs=1e-15)
