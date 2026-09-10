#!/usr/bin/env python3
"""Golden net for ``fop_pool.pool_domains`` (scoring_v2 core v3 — treeless pooler).

Freezes the FULL return dict of ``pool_domains`` over the hand-worked scenarios in
``docs/scoring_v3_core.md`` Appendix B (mean-of-means over the K fixed Voronoi
domains; PSS enters only at the pooling step; ``M = 1`` degenerates to the plain
PSS-weighted mean). Expected values authored **by hand from the doc** (V3-0).

Marked ``xfail(strict=True)`` until V3-2: ``pool_domains`` does not exist yet, so
the call raises ``AttributeError`` -> xfail passes. V3-2 lands it, regenerates the
fixture, and drops the marker.

Run:  python -m pytest test_pool_domains_golden.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
from _common import _load, diff_report, load_json  # noqa: E402

fop = _load("fop_pool")
_GOLDEN = load_json("pool_domains_golden.json")
TOL = 1e-9


def _norm(d):
    if isinstance(d, dict):
        return {str(k): _norm(v) for k, v in d.items()}
    if isinstance(d, list):
        return [_norm(v) for v in d]
    return d


def _pss_map(rows):
    return {(h, int(dom)): float(w) for h, dom, w in rows} or None


@pytest.mark.xfail(strict=True, reason="pool_domains lands in V3-2")
def test_pool_domains_golden():
    failures: list[str] = []
    for entry in _GOLDEN:
        got = fop.pool_domains(entry["hyp_records"], _pss_map(entry["pss"]))
        failures += diff_report(
            _norm(entry["expected"]), _norm(got), TOL, entry["name"]
        )
    assert not failures, "\n".join(failures)


def test_golden_covers_the_enumerated_cases():
    """Doc-independent: every Appendix-B scenario is present. Stays green."""
    names = {e["name"] for e in _GOLDEN}
    required = {
        "m1_passthrough", "m1_pss_shifts_to_d1", "m2_domain_changed_only_in_H1",
        "m2_asymmetric_pss_same_domain", "side_none_collapse",
        "no_reconstruction_in_denominator", "m3_split_scheme_resolved_US",
        "m3_split_scheme_resolved_GS3",
    }
    assert required <= names, f"missing golden scenarios: {required - names}"
