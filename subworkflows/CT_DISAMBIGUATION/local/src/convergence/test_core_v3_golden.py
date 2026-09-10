#!/usr/bin/env python3
"""Golden net for ``path_scores.compute_domain_scores`` (scoring_v2 core v3).

Freezes the FULL return dict of ``compute_domain_scores`` over the hand-worked
scenarios in ``docs/scoring_v3_core.md`` Appendix A. Expected values live in
``golden/core_v3_golden.json`` and were authored **by hand from the doc** (V3-0):
V3-1 lands ``compute_domain_scores`` and switches the fixture to generator-output
— the regen diff against this file must be empty (proof the impl matches the
spec).

Marked ``xfail(strict=True)`` until V3-1: ``compute_domain_scores`` does not exist
yet, so the call raises ``AttributeError`` -> xfail passes. Once V3-1 lands a
correct implementation the test xpasses, ``strict`` turns that into a failure,
and V3-1 drops the marker.

Run:  python -m pytest test_core_v3_golden.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
from _common import _load, build_tree, diff_report, load_json  # noqa: E402

ps = _load("path_scores")
_GOLDEN = load_json("core_v3_golden.json")
TOL = 1e-9


def _norm(d):
    if isinstance(d, dict):
        return {str(k): _norm(v) for k, v in d.items()}
    if isinstance(d, list):
        return [_norm(v) for v in d]
    return d


def _run(scn):
    root, _ = build_tree([tuple(e) for e in scn["edges"]], list(scn["node_ids"]))
    node_index = ps.build_node_index(root)
    pnd = {int(k): v for k, v in scn["posteriors"].items()}
    return ps.compute_domain_scores(
        scn["pair_details"], pnd, node_index, scn["scheme"]
    )


@pytest.mark.xfail(strict=True, reason="compute_domain_scores lands in V3-1")
def test_core_v3_golden():
    failures: list[str] = []
    for entry in _GOLDEN:
        got = _run(entry["scenario"])
        failures += diff_report(
            _norm(entry["expected"]), _norm(got), TOL, entry["scenario"]["name"]
        )
    assert not failures, "\n".join(failures)


def test_golden_covers_the_enumerated_cases():
    """Doc-independent: the fixture carries every Appendix-A scenario. Stays green
    (it only reads the JSON)."""
    names = {e["scenario"]["name"] for e in _GOLDEN}
    required = {
        "one_domain_changed", "two_domains_same_residue_clean",
        "two_domains_diff_residue_US", "two_domains_coencoded_GS3",
        "two_domains_opposite_sides", "three_domains_VVL_majority",
        "shared_lca_contaminated", "no_domain_changed",
        "domain_without_reconstruction",
    }
    assert required <= names, f"missing golden scenarios: {required - names}"
