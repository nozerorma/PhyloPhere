"""core.pooling.pooled_sides, and the parity of the observed and null emitters that use it (M3)."""
import sys
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.convergence.disambiguate_single import _emit_pooled_side_rows  # noqa: E402
from src.convergence.fop_pool import pool_domains  # noqa: E402
from src.core.driver import pool_labelings  # noqa: E402
from src.core.pooling import pooled_sides  # noqa: E402
from src.data.models import ConvergenceResult  # noqa: E402


def _sides(top, bottom=None):
    """compute_domain_scores-shaped record: scores per domain and the residue each changed domain ended on."""
    bottom = bottom or {}
    return {
        "top": {"domain_scores": top, "domain_der_enc": {d: "D" for d in top}, "domain_der": {d: "D" for d in top},
                "domain_anc": {d: "A" for d in top}},
        "bottom": {"domain_scores": bottom or {1: 0.0, 2: 0.0}, "domain_der_enc": {d: "K" for d in bottom}},
        "domain_meta": {1: {}, 2: {}},
    }


def test_pooled_sides_reports_only_sides_with_a_changed_domain():
    pooled = pool_domains([{"hyp": "H1", "sides": _sides({1: 1.0, 2: 0.5})}], None)
    got = pooled_sides(pooled)
    assert [g["side"] for g in got] == ["top"]
    assert got[0]["asr_path_score"] == pytest.approx(0.75) and got[0]["derived_agreement"] == 1.0
    assert got[0]["participating_hyps"] == "H1" and got[0]["domain_scores"] == {1: 1.0, 2: 0.5}
    assert "convergence_type" in got[0]
    assert pooled_sides(pool_domains([{"hyp": "H1", "sides": _sides({})}], None)) == []
    assert pooled_sides(None) == []


def test_observed_and_null_emitters_agree_per_side():
    """Same hypotheses, same weights: the observed rows and the null rows carry the same numbers."""
    hyps = {"H1": _sides({1: 1.0, 2: 0.0}, {1: 0.5}), "H2": _sides({1: 1.0, 2: 0.5}, {1: 0.5})}
    pss = {("H1", 1): 3.0, ("H1", 2): 1.0, ("H2", 1): 1.0, ("H2", 2): 1.0}
    base = ConvergenceResult(gene="G", position=5, tag="POS5", caas="D/K", ancestral="A", derived="D",
                             convergence_type="divergent")
    hyp_rows = [{"hyp": h, "sides": s} for h, s in hyps.items()]
    observed = _emit_pooled_side_rows(base, hyp_rows, pss, all_rows=[base, base])

    records = [NS(position=5, caap_group="US", sides=s, hypothesis=None) for s in hyps.values()]
    null = pool_labelings([(f"b_1~{h}", [r]) for h, r in zip(hyps, records)], {"b_1": pss})
    assert [c for c, _ in null] == ["b_1"]

    assert [r.side for r in observed] == [r.side for r in null[0][1]] == ["top", "bottom"]
    for o, n in zip(observed, null[0][1]):
        assert o.asr_path_score == n.asr_path_score
        assert o.derived_agreement == n.derived_agreement
        assert o.domain_scores == n.domain_scores
