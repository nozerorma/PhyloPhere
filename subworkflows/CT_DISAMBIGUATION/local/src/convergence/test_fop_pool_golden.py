#!/usr/bin/env python3
"""T0 shared-oracle golden for fop_pool.py <-> fop_pool.R (roadmap tier T0, H8).

``golden/fop_pool_fixture.json`` holds a handful of hypothesis-pooling scenarios
plus the frozen ``pool_hypotheses`` output. This test and its R twin
(``subworkflows/SCORING/local/src/test_fop_pool_golden.R``) both read that ONE
file and assert their pooling reproduces ``expected`` to 1e-9 — the shared oracle
that catches the R/Python twins drifting (they already have, twice).

Run:  python -m pytest test_fop_pool_golden.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
from _common import _load, load_json  # noqa: E402

fop = _load("fop_pool")
_FIX = load_json("fop_pool_fixture.json")
TOL = 1e-9
_KEYS = ("asr_path_score", "core", "independence", "mrca_diversity",
         "derived_agreement", "conservation_gate")


def _keys_to_int(d):
    return {int(k): v for k, v in d.items()}


def _prep(records):
    fields = ("pair_scores", "pair_top_scores", "pair_bottom_scores",
              "pair_derived_top", "pair_derived_bot",
              "conserved_pair_scores", "conserved_pair_nodes")
    return [
        {**r, **{k: _keys_to_int(r[k]) for k in fields if k in r}}
        for r in records
    ]


def _pss(rows):
    return {(h, int(d)): float(v) for h, d, v in rows} or None


def test_fop_pool_golden():
    failures = []
    for sc in _FIX:
        got = fop.pool_hypotheses(_prep(sc["hyp_records"]), _pss(sc["pss"]),
                                  scheme=sc["scheme"])
        for k in _KEYS:
            a, b = sc["expected"].get(k), got.get(k)
            if a is None and b is None:
                continue
            if a is None or b is None or abs(float(a) - float(b)) > TOL:
                failures.append(f"{sc['name']}.{k}: {b!r} != {a!r}")
    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    try:
        test_fop_pool_golden()
        print(f"ALL {len(_FIX)} FOP GOLDEN SCENARIOS MATCH")
    except AssertionError as e:
        print(e)
        sys.exit(1)
