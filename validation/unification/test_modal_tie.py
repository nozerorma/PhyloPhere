"""A tie between residues of a domain is settled by a rule that does not depend on the order of the hypotheses, and it is reported."""
import itertools
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.convergence.fop_pool import _modal_str, pool_domains  # noqa: E402
from src.core.master import master_fields  # noqa: E402
from src.core.pooling import pooled_sides  # noqa: E402


def _rec(hyp, top):
    """One hypothesis whose top side changed the domains of `top` = {domain: (encoded residue, raw derived residue)}."""
    return {"hyp": hyp, "sides": {
        "top": {"domain_scores": {d: 0.5 for d in top}, "domain_der_enc": {d: e for d, (e, _) in top.items()},
                "domain_der": {d: r for d, (_, r) in top.items()}, "domain_anc": {d: "A" for d in top}},
        "bottom": {"domain_scores": {}, "domain_der_enc": {}, "domain_der": {}, "domain_anc": {}},
        "domain_meta": {d: {} for d in top}}}


# domain 1 is tied (Y in H1, X in H2); domain 2 holds X everywhere: the tie decides whether the two domains agree
TIED = [_rec("H1", {1: ("Y", "y"), 2: ("X", "x")}), _rec("H2", {1: ("X", "x"), 2: ("X", "x")})]
UNTIED = [_rec("H1", {1: ("X", "x"), 2: ("X", "x")}), _rec("H2", {1: ("X", "x"), 2: ("X", "x")})]


def test_the_modal_residue_is_the_most_frequent_and_a_tie_goes_to_the_smallest_whatever_the_order():
    assert _modal_str(["Y", "Y", "X"]) == "Y"
    assert _modal_str(["B", "A"]) == _modal_str(["A", "B"]) == "A"
    assert _modal_str(["", None, "Q"]) == "Q" and _modal_str([]) is None and _modal_str(["", None]) is None


@pytest.mark.parametrize("order", list(itertools.permutations(range(3))))
def test_a_tied_domain_gives_the_same_pooled_summary_in_every_order_of_the_hypotheses(order):
    recs = TIED + [_rec("H3", {1: ("Y", "y"), 2: ("X", "x")})]          # H1 Y, H2 X, H3 Y: not a tie once H3 votes
    recs = [recs[i] for i in order]
    got = pooled_sides(pool_domains(recs, None))
    ref = pooled_sides(pool_domains(sorted(recs, key=lambda r: r["hyp"]), None))
    assert got == ref


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_a_two_way_tie_does_not_change_the_label_with_the_order(order):
    out = pooled_sides(pool_domains([TIED[i] for i in order], None))[0]
    assert out["convergence_type"] == "convergent" and out["derived_agreement"] == 1.0   # X wins the tie in domain 1: both domains hold X
    assert out["domain_der"] == {1: "x", 2: "x"}


def test_the_tie_is_reported_and_an_unanimous_pool_is_not():
    assert pooled_sides(pool_domains(TIED, None))[0]["agreement_ambiguous"] is True
    assert pooled_sides(pool_domains(UNTIED, None))[0]["agreement_ambiguous"] is False


def test_the_master_has_the_ambiguity_column_next_to_the_agreement():
    fields = master_fields(2)
    assert fields.index("agreement_ambiguous") == fields.index("derived_agreement") + 1
