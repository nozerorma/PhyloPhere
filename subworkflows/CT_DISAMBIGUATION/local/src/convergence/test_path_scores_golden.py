#!/usr/bin/env python3
"""T0 golden net for compute_asr_path_score (roadmap tier T0).

Freezes the FULL return dict of ``compute_asr_path_score`` over a set of
synthetic scenarios (one changed pair, >=2 same-side, opposite-only, both-sides,
conserved pair, hop+1 contamination, n>2 mixed residues, GS co-encoding, MRCA at
root, sibling merge, no changed pairs, soft posteriors).

The expected values live in ``golden/path_scores_golden.json`` and are produced
by ``golden/gen_golden.py`` — every tier that changes the scoring maths re-runs
the generator and commits the diff instead of hand-editing numbers here.

Run:  python -m pytest test_path_scores_golden.py
  or: python test_path_scores_golden.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
from _common import _load, diff_report, load_json, run_scenario  # noqa: E402

ps = _load("path_scores")
_GOLDEN = load_json("path_scores_golden.json")
TOL = 1e-12


def _check(entry) -> list[str]:
    got = run_scenario(entry["scenario"], ps)
    # json round-trips dict int keys to str; normalise both sides for the diff
    return diff_report(_norm(entry["expected"]), _norm(got), TOL,
                       entry["scenario"]["name"])


def _norm(d):
    if isinstance(d, dict):
        return {str(k): _norm(v) for k, v in d.items()}
    if isinstance(d, list):
        return [_norm(v) for v in d]
    return d


def test_path_scores_golden():
    failures: list[str] = []
    for entry in _GOLDEN:
        failures += _check(entry)
    assert not failures, "\n".join(failures)


def test_golden_covers_the_enumerated_cases():
    names = {e["scenario"]["name"] for e in _GOLDEN}
    required = {
        "single_changed_pair", "two_same_side_converge", "opposite_sides_only",
        "pair_changes_both_sides", "with_conserved_pair", "contaminated_hop1",
        "n_gt_2_mixed_residues", "mrca_at_root", "sibling_merge",
        "no_changed_pairs",
    }
    assert required <= names, f"missing golden scenarios: {required - names}"


def test_conserved_pair_columns_survive():
    """H2 guard: the conserved-pair plumbing (roadmap decision G) must keep
    emitting its maps even after T1 removes conservation_gate from the product."""
    entry = next(e for e in _GOLDEN
                 if e["scenario"]["name"] == "with_conserved_pair")
    got = run_scenario(entry["scenario"], ps)
    for key in ("conserved_pair_scores", "conserved_pair_nodes",
                "pair_ancestral", "pair_derived_top", "pair_derived_bot"):
        assert key in got, f"{key} dropped from compute_asr_path_score return"
    assert got["conserved_pair_scores"], "conserved pair 3 lost its score"
    assert set(got["conserved_pair_scores"]) == set(got["conserved_pair_nodes"])


if __name__ == "__main__":
    fails = []
    for e in _GOLDEN:
        fails += _check(e)
    if fails:
        print("\n".join(fails))
        print("\nSOME GOLDEN CHECKS FAILED")
        sys.exit(1)
    print(f"ALL {len(_GOLDEN)} GOLDEN SCENARIOS MATCH")
