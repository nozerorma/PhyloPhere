#!/usr/bin/env python3
"""V3-3: the permulation null and the observed path score one statistic.

``_perms_worker`` (null) and ``_emit_pooled_side_rows`` (observed) are now both
thin wrappers over ``fop_pool.pool_domains``. The FCS ``p.perm`` comparison is
only calibrated if, for the SAME harvest + PSS, both routes emit the same
``core_s`` per side and collapse to ``side="none"`` under the same condition
(risk R-null-1 in the V3-3 brief). This replays every Appendix-B scenario
through both and asserts bit-for-bit parity on the direction key + score.

Run:  python -m pytest test_null_domain_pool_wiring.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))

from _common import _load, load_json  # noqa: E402
from src.data.models import ConvergenceResult  # noqa: E402
from src.convergence.disambiguate_single import _emit_pooled_side_rows  # noqa: E402

fop = _load("fop_pool")
_GOLDEN = load_json("pool_domains_golden.json")
TOL = 1e-12


def _pss_map(rows):
    return {(h, int(dom)): float(w) for h, dom, w in rows} or None


def _null_expand(pooled):
    """Mirror of ``_perms_worker._expand_pooled``: pooled -> {side: core}, or
    {"none": 0.0} when no domain participates on either side."""
    out = {}
    for s in ("top", "bottom"):
        agg = pooled.get(s) or {}
        if int(agg.get("n_participating", 0) or 0) > 0:
            out[s] = float(agg.get("core", agg.get("asr_path_score", 0.0)) or 0.0)
    return out or {"none": 0.0}


def test_null_and_observed_agree_per_scenario():
    base = ConvergenceResult(
        gene="G", position=10, tag="t", caas="A/B", ancestral="A", derived="B",
        convergence_type="no_change", caap_group="US", hypothesis="H1",
    )
    failures = []
    for entry in _GOLDEN:
        name = entry["name"]
        pss = _pss_map(entry["pss"])
        pooled = fop.pool_domains(entry["hyp_records"], pss)
        null_sides = _null_expand(pooled)

        obs_rows = _emit_pooled_side_rows(base, entry["hyp_records"], pss)
        obs_sides = {r.side: float(r.core or 0.0) for r in obs_rows}

        if set(obs_sides) != set(null_sides):
            failures.append(f"{name}: side sets differ obs={set(obs_sides)} "
                            f"null={set(null_sides)}")
            continue
        for s in obs_sides:
            if abs(obs_sides[s] - null_sides[s]) > TOL:
                failures.append(f"{name}/{s}: obs core {obs_sides[s]} != "
                                f"null core {null_sides[s]}")
            # observed rows also carry asr_path_score == core
            for r in obs_rows:
                if r.side == s and abs(float(r.asr_path_score or 0.0)
                                      - null_sides[s]) > TOL:
                    failures.append(f"{name}/{s}: observed asr_path_score "
                                    f"{r.asr_path_score} != core {null_sides[s]}")
    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
