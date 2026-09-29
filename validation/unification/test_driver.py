"""core.driver.pool_labelings: hypotheses of a cycle collapse into one domain-pooled record set."""
import sys
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.driver import pool_labelings  # noqa: E402


def _rec(top_scores, bottom_changed=False, position=5, group="US"):
    """One axes-only record: `sides` as compute_domain_scores returns it (two Voronoi domains)."""
    top = {"domain_scores": top_scores, "domain_der_enc": {d: "D" for d in top_scores}}
    bottom = {"domain_scores": {1: 0.0, 2: 0.0}, "domain_der_enc": {1: "K"} if bottom_changed else {}}
    return NS(position=position, caap_group=group, hypothesis=None,
              sides={"top": top, "bottom": bottom, "domain_meta": {1: {}, 2: {}}})


def test_fop_pools_hypotheses_with_pss_weights():
    # domain scores per hypothesis and their PSS weights; expected core = sum(w*s)/sum(w) over domains
    h1, h2 = _rec({1: 1.0, 2: 0.0}), _rec({1: 1.0, 2: 0.5})
    pss = {"b_1": {("H1", 1): 3.0, ("H1", 2): 1.0, ("H2", 1): 1.0, ("H2", 2): 1.0}}
    s_bar = {1: (1.0 + 1.0) / 2, 2: (0.0 + 0.5) / 2}
    w_bar = {1: (3.0 + 1.0) / 2, 2: (1.0 + 1.0) / 2}
    expected = sum(w_bar[d] * s_bar[d] for d in s_bar) / sum(w_bar.values())
    out = pool_labelings([("b_1~H1", [h1]), ("b_1~H2", [h2])], pss)
    assert [c for c, _ in out] == ["b_1"]  # the ~H<m> hypotheses collapse into their base cycle
    rows = out[0][1]
    assert [(r.position, r.caap_group, r.side) for r in rows] == [(5, "US", "top")]  # bottom has no change
    assert rows[0].asr_path_score == pytest.approx(expected, abs=1e-15)
    assert rows[0].asr_path_score == pytest.approx(0.75, abs=1e-15)


def test_fop_without_pss_weights_every_domain_equally():
    out = pool_labelings([("b_2~H1", [_rec({1: 1.0, 2: 0.0})]), ("b_2~H2", [_rec({1: 0.0, 2: 0.0})])], {})
    assert out[0][1][0].asr_path_score == pytest.approx(0.25)  # mean over two hypotheses and two domains


def test_non_fop_keeps_each_labeling_a_cycle_and_reports_both_sides():
    out = pool_labelings([("b_3", [_rec({1: 1.0, 2: 1.0}, bottom_changed=True)]), ("b_4", [_rec({1: 0.0, 2: 0.0})])], None)
    assert [c for c, _ in out] == ["b_3", "b_4"]
    assert [r.side for r in out[0][1]] == ["top", "bottom"]
    assert out[0][1][0].asr_path_score == pytest.approx(1.0)


def test_position_without_any_changed_domain_is_one_none_row():
    quiet = NS(position=7, caap_group="GS1", hypothesis=None,
               sides={"top": {"domain_scores": {}}, "bottom": {"domain_scores": {}}, "domain_meta": {1: {}}})
    out = pool_labelings([("b_5", [quiet])], None)
    assert [(r.position, r.caap_group, r.side, r.asr_path_score) for r in out[0][1]] == [(7, "GS1", "none", 0.0)]
