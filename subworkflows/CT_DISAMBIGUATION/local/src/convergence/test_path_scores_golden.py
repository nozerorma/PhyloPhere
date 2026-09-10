#!/usr/bin/env python3
"""Golden net for compute_asr_path_score (scoring_v2 — per-side pairwise core).

Freezes the FULL return dict of ``compute_asr_path_score`` over the synthetic
scenarios in ``golden/gen_golden.py``. Every scenario returns the per-side
``{"top": <row>, "bottom": <row>}`` shape — T4a retired the flat one-row path
and its ``native_side_split`` flag.

Expected values live in ``golden/path_scores_golden.json`` and are produced by
``golden/gen_golden.py`` — every tier that changes the maths re-runs the
generator and commits the diff instead of hand-editing numbers here.

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


def _norm(d):
    if isinstance(d, dict):
        return {str(k): _norm(v) for k, v in d.items()}
    if isinstance(d, list):
        return [_norm(v) for v in d]
    return d


def _check(entry) -> list[str]:
    got = run_scenario(entry["scenario"], ps)
    return diff_report(_norm(entry["expected"]), _norm(got), TOL,
                       entry["scenario"]["name"])


def _by_name(name):
    return next(e for e in _GOLDEN if e["scenario"]["name"] == name)


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
        "no_changed_pairs", "both_sides_two_rows",
    }
    assert required <= names, f"missing golden scenarios: {required - names}"


def test_return_shape_is_always_per_side():
    """Every scenario returns {"top": row, "bottom": row}; the retired global
    axes are gone from each row (T3-doc §13)."""
    got = run_scenario(_by_name("with_conserved_pair")["scenario"], ps)
    assert set(got) == {"top", "bottom"}
    for side in ("top", "bottom"):
        for gone in ("independence", "replication", "mrca_diversity",
                     "conservation_gate", "strength", "core_top", "core_bottom"):
            assert gone not in got[side], f"{gone} still emitted in {side}"


def test_native_side_split_returns_two_rows():
    """{"top": row, "bottom": row}; sides never recombine (no 1-(1-t)(1-b)); the
    conserved pair counts in |D_s| on BOTH sides."""
    got = run_scenario(_by_name("both_sides_two_rows")["scenario"], ps)
    assert set(got) == {"top", "bottom"}
    top, bottom = got["top"], got["bottom"]
    assert abs(top["asr_path_score"] - 0.5158540) < 1e-6
    assert abs(bottom["asr_path_score"] - 0.5430045) < 1e-6
    assert top["n_pairs_side"] == 3 and bottom["n_pairs_side"] == 3
    assert top["n_conserved"] == 1 and bottom["n_conserved"] == 1
    assert top["n_participating"] == 2 and bottom["n_participating"] == 2
    assert abs(top["conserved_pair_scores"][4] - 0.90) < 1e-9
    union = 1.0 - (1.0 - top["asr_path_score"]) * (1.0 - bottom["asr_path_score"])
    assert abs(max(top["asr_path_score"], bottom["asr_path_score"]) - union) > 1e-3


def test_t3_core_pareado():
    """Pin core_top / core_bottom of the two scenarios the T3 rewrite moves
    (T3-doc §10): with_conserved_pair drops (conserved pair now in the
    denominator), both_sides_two_rows splits into two independent rows."""
    wcp = run_scenario(_by_name("with_conserved_pair")["scenario"], ps)
    # 2 converging pairs (0.7737809 each) + 1 conserved (0) over |D_top| = 3
    assert abs(wcp["top"]["asr_path_score"] - (2 * 0.7737809 / 3)) < 1e-6
    assert wcp["bottom"]["asr_path_score"] == 0.0

    bs = run_scenario(_by_name("both_sides_two_rows")["scenario"], ps)
    assert abs(bs["top"]["asr_path_score"] - (2 * 0.7737809 / 3)) < 1e-6
    assert abs(bs["bottom"]["asr_path_score"] - (2 * 0.8145063 / 3)) < 1e-6


def test_n_gt_2_majority_no_longer_penalised():
    """T3-doc §10: the majority (3x->V) convergence rises once the lone ->T pair
    stops dragging a position-wide derived_agreement and instead just adds a 0
    to the denominator."""
    got = run_scenario(_by_name("n_gt_2_mixed_residues")["scenario"], ps)
    assert abs(got["top"]["asr_path_score"] - 0.74135) < 1e-4


def test_conserved_pair_columns_survive():
    """H2 guard (roadmap decision G): the conserved-pair plumbing must keep
    emitting its maps in each per-side row."""
    got = run_scenario(_by_name("with_conserved_pair")["scenario"], ps)
    top = got["top"]
    for key in ("conserved_pair_scores", "conserved_pair_nodes",
                "pair_ancestral", "pair_derived_top", "pair_derived_bot"):
        assert key in top, f"{key} dropped from compute_asr_path_score return"
    assert top["conserved_pair_scores"], "conserved pair 3 lost its score"
    assert set(top["conserved_pair_scores"]) == set(top["conserved_pair_nodes"])


if __name__ == "__main__":
    fails = []
    for e in _GOLDEN:
        fails += _check(e)
    if fails:
        print("\n".join(fails))
        print("\nSOME GOLDEN CHECKS FAILED")
        sys.exit(1)
    print(f"ALL {len(_GOLDEN)} GOLDEN SCENARIOS MATCH")
