#!/usr/bin/env python3
"""Regression tests for the hypothesis-pooling representative-collapse fixes.

Covers:
  1. ``_emit_pooled_side_rows``: the ``trait`` "H1"-collision fix + the new
     ``tag_support``/``caas_support``/``participating_hypotheses`` columns.
  2. ``pool_domains``: the new ``domain_der_support``/``domain_anc_support``
     support-string columns (siblings of the unchanged modal ``domain_der``).
  3. ``score_domains_side``: the new ``pair_lca`` capture.
  4. ``build_result_from_row`` (CT_DISAMBIGUATION debug-tree plotting): the
     stale ``mrca_{idx}_*`` -> ``domain_{idx}_*`` column-name fix.

Run:  python -m pytest test_pooling_representative_collapse.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.data.models import ConvergenceResult  # noqa: E402
from src.convergence.disambiguate_single import _emit_pooled_side_rows  # noqa: E402
from src.convergence.fop_pool import pool_domains  # noqa: E402
from src.convergence.path_scores import (  # noqa: E402
    build_node_index, score_domains_side,
)
from src.plots.plot_utils import build_result_from_row  # noqa: E402


def _mk(**kw) -> ConvergenceResult:
    base = dict(gene="G", position=10, tag="t", caas="A/B", ancestral="A",
                derived="B", convergence_type="convergent")
    base.update(kw)
    return ConvergenceResult(**base)


# ── 1. _emit_pooled_side_rows ────────────────────────────────────────────────

def test_emit_pooled_side_rows_trait_collision_and_support_columns():
    """3-row pool: hypotheses 'H1', 'H2', and an unresolved (None) row, with
    differing tag/caas. Before the fix, the unresolved row's `hyp` fell back to
    the literal string "H1" and could collide with the genuine "H1" label,
    falsely collapsing `len(hyp_labels)` to 1 and letting `trait` pass through
    an arbitrary hypothesis instead of nulling on real multi-hypothesis pools.
    """
    base = _mk(hypothesis="H1", tag="tagA", caas="A/B")
    rows = [
        _mk(hypothesis="H1", tag="tagA", caas="A/B"),
        _mk(hypothesis="H2", tag="tagB", caas="A/C"),
        _mk(hypothesis=None, tag="tagC", caas="A/D"),
    ]
    # A domain changed under H1 and H2 only (the unresolved row contributes no
    # changed domain -> it must not appear in `participating_hypotheses`).
    sides_h1 = {"top": {"domain_scores": {"1": 0.9}, "domain_der": {"1": "V"},
                         "domain_der_enc": {"1": "V"}, "domain_anc": {"1": "A"},
                         "agree_num": 1, "agree_den": 1, "n_changed": 1,
                         "convergence_type": "single", "pair_lca": []},
                "bottom": {"domain_scores": {}, "domain_der": {}, "domain_der_enc": {},
                           "domain_anc": {}, "agree_num": 0, "agree_den": 0,
                           "n_changed": 0, "convergence_type": "no_change", "pair_lca": []},
                "domain_meta": {"1": {"mrca_id": 5, "state": "A", "posterior": 0.9}}}
    sides_h2 = {"top": {"domain_scores": {"1": 0.8}, "domain_der": {"1": "V"},
                         "domain_der_enc": {"1": "V"}, "domain_anc": {"1": "A"},
                         "agree_num": 1, "agree_den": 1, "n_changed": 1,
                         "convergence_type": "single", "pair_lca": []},
                "bottom": {"domain_scores": {}, "domain_der": {}, "domain_der_enc": {},
                           "domain_anc": {}, "agree_num": 0, "agree_den": 0,
                           "n_changed": 0, "convergence_type": "no_change", "pair_lca": []},
                "domain_meta": {"1": {"mrca_id": 5, "state": "A", "posterior": 0.9}}}
    sides_none = {"top": {"domain_scores": {}, "domain_der": {}, "domain_der_enc": {},
                           "domain_anc": {}, "agree_num": 0, "agree_den": 0,
                           "n_changed": 0, "convergence_type": "no_change", "pair_lca": []},
                  "bottom": {"domain_scores": {}, "domain_der": {}, "domain_der_enc": {},
                             "domain_anc": {}, "agree_num": 0, "agree_den": 0,
                             "n_changed": 0, "convergence_type": "no_change", "pair_lca": []},
                  "domain_meta": {"1": {"mrca_id": 5, "state": "A", "posterior": 0.9}}}
    rows[0].sides = sides_h1
    rows[1].sides = sides_h2
    rows[2].sides = sides_none

    hyp_rows = [
        {"hyp": (getattr(r, "hypothesis", None) or f"_H_UNRESOLVED_{i}"),
         "sides": getattr(r, "sides", None) or {}}
        for i, r in enumerate(rows)
    ]
    out = _emit_pooled_side_rows(base, hyp_rows, all_rows=rows)

    # Genuine multi-hypothesis pool (3 distinct labels, no collision) -> trait nulled.
    assert all(r.hypothesis is None for r in out)

    top_row = next(r for r in out if r.side == "top")
    assert top_row.tag_support == "tagA:1,tagB:1,tagC:1"
    assert top_row.caas_support == "A/B:1,A/C:1,A/D:1"
    # Only H1/H2 drove a changed domain on top -> the unresolved row is excluded.
    assert top_row.participating_hypotheses == "H1,H2"


def test_emit_pooled_side_rows_single_hypothesis_keeps_label():
    """A lone non-FOP row (hypothesis=None) still reduces to one label; `trait`
    stays empty as before (unchanged regression check)."""
    base = _mk(hypothesis=None)
    rows = [_mk(hypothesis=None)]
    rows[0].sides = {
        "top": {"domain_scores": {"1": 0.5}, "domain_der": {"1": "V"},
                "domain_der_enc": {"1": "V"}, "domain_anc": {"1": "A"},
                "agree_num": 1, "agree_den": 1, "n_changed": 1,
                "convergence_type": "single", "pair_lca": []},
        "bottom": {"domain_scores": {}, "domain_der": {}, "domain_der_enc": {},
                   "domain_anc": {}, "agree_num": 0, "agree_den": 0,
                   "n_changed": 0, "convergence_type": "no_change", "pair_lca": []},
        "domain_meta": {"1": {"mrca_id": 5, "state": "A", "posterior": 0.9}},
    }
    hyp_rows = [{"hyp": "_H_UNRESOLVED_0", "sides": rows[0].sides}]
    out = _emit_pooled_side_rows(base, hyp_rows, all_rows=rows)
    top_row = next(r for r in out if r.side == "top")
    assert top_row.hypothesis is None  # single label -> base.hypothesis (None) passes through


# ── 2. pool_domains domain support columns ───────────────────────────────────

def test_pool_domains_support_columns_additive():
    """2-hypothesis case where domain 3 reconstructs 'V' in one hypothesis and
    'I' in another: the new support column must show both, while the existing
    modal `domain_der` stays first-seen-wins (unchanged regression check)."""
    hyp_records = [
        {"hyp": "H1", "sides": {
            "top": {"domain_scores": {"3": 0.9}, "domain_der": {"3": "V"},
                    "domain_der_enc": {"3": "V"}, "domain_anc": {"3": "A"}},
            "bottom": {},
            "domain_meta": {"3": {"mrca_id": 1}},
        }},
        {"hyp": "H2", "sides": {
            "top": {"domain_scores": {"3": 0.7}, "domain_der": {"3": "I"},
                    "domain_der_enc": {"3": "I"}, "domain_anc": {"3": "A"}},
            "bottom": {},
            "domain_meta": {"3": {"mrca_id": 1}},
        }},
    ]
    out = pool_domains(hyp_records)
    assert out["top"]["domain_der_support"]["3"] == "I:1,V:1"  # alphabetical tiebreak
    assert out["top"]["domain_der"]["3"] == "V"  # unchanged modal (first-seen-wins)


# ── 3. score_domains_side pair_lca capture ───────────────────────────────────

def test_score_domains_side_captures_pair_lca():
    """2-domain-pair case: pair_lca must carry exactly the expected
    (domain_a, domain_b, lca_node_id, contrib) tuple, and domain_scores must be
    unchanged from before the pair_lca addition (purely additive)."""
    # Tiny tree: root(0) -> {1, 2}; each domain's MRCA is a leaf under root.
    class _Node:
        def __init__(self, node_id, children=None):
            self.node_id = node_id
            self.children = children or []
            self.parent = None

    n1, n2 = _Node(1), _Node(2)
    root = _Node(0, [n1, n2])
    n1.parent = root
    n2.parent = root
    node_index = build_node_index(root)

    per_node_dist = {0: {"V": 0.05, "A": 0.95}}
    domains = [
        {"d": 1, "mrca_id": 1, "anc_enc": "A", "der_enc": "V"},
        {"d": 2, "mrca_id": 2, "anc_enc": "A", "der_enc": "V"},
    ]
    out = score_domains_side(domains, node_index, per_node_dist, scheme="US")
    assert len(out["pair_lca"]) == 1
    a, b, lca, contrib = out["pair_lca"][0]
    assert (a, b) == (1, 2)
    assert lca == 0
    assert abs(contrib - 0.95) < 1e-9
    # domain_scores unaffected by the addition.
    assert set(out["domain_scores"]) == {1, 2}


# ── 4. build_result_from_row stale-column fix ────────────────────────────────

def test_build_result_from_row_reads_domain_columns():
    """Regression: `mrca_{idx}_*` no longer exists in the schema (renamed to
    `domain_{idx}_*`); before the fix this always returned an empty
    focal_nodes/node_mapping. A row using the current `domain_{idx}_*` schema
    must now populate them correctly."""
    row = pd.Series({
        "msa_pos": 9,
        "all_mrca_node": 0,
        "all_mrca_state": "A",
        "all_mrca_posterior": 0.9,
        "domain_1_node": 1,
        "domain_1_state": "V",
        "domain_1_posterior": 0.8,
    })
    result = build_result_from_row(row)
    assert result is not None
    assert result["node_mapping"]["focal_nodes"] == [1]
    assert result["node_mapping"]["focal_1"] == 1
    assert result["node_state_details"]["focal_states"] == ["V"]


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
