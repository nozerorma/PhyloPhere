#!/usr/bin/env python3
"""T3b: two-row observed plumbing for --native_side_split.

`compute_asr_path_score(native_side_split=True)` returns {"top": row, "bottom":
row}; `disambiguate_single._split_result_by_side` expands one position-level
ConvergenceResult into one row per *participating* side (a "both" position -> two
rows keyed by `side`), and `convert_convergence_result_to_dict` carries each
row's `side` + per-side scalar through to the flat CSV dict.

Run:  python -m pytest test_native_side_split.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))

from _common import _load, build_tree  # noqa: E402
from src.data.models import ConvergenceResult  # noqa: E402
from src.convergence.disambiguate_single import _split_result_by_side  # noqa: E402
from src.utils.gene_wrapper import convert_convergence_result_to_dict  # noqa: E402

ps = _load("path_scores")


# §10.1 fixture from docs/scoring_v2_T3_core_pareado.md (same as the T3a golden).
_BS_EDGES = [(0, 1), (0, 2), (1, 10), (10, 3), (10, 4), (2, 20), (20, 5),
             (20, 6), (2, 21), (21, 9), (21, 12), (0, 30), (30, 7), (30, 8)]
_BS_NODES = [0, 1, 2, 10, 20, 21, 30, 3, 4, 5, 6, 9, 12, 7, 8]
_BS_POST = {n: {"A": 0.90, "V": 0.05, "L": 0.05}
            for n in [0, 1, 2, 10, 20, 21, 30]}
_BS_PAIRS = [
    {"pair_id": 1, "node_id": 3, "focal_state": "A",
     "top_tip_mode": "V", "bottom_tip_mode": "L"},   # both sides, diff residues
    {"pair_id": 2, "node_id": 5, "focal_state": "A",
     "top_tip_mode": "V", "bottom_tip_mode": "A"},    # top only
    {"pair_id": 3, "node_id": 7, "focal_state": "A",
     "top_tip_mode": "A", "bottom_tip_mode": "L"},    # bottom only
    {"pair_id": 4, "node_id": 9, "focal_state": "A",
     "top_tip_mode": "A", "bottom_tip_mode": "A"},    # conserved
]


def _split_dict(pairs=_BS_PAIRS, conserved="4", meta=True):
    root, _ = build_tree(_BS_EDGES, _BS_NODES)
    node_index = ps.build_node_index(root)
    return ps.compute_asr_path_score(
        pairs, dict(_BS_POST), node_index, "US", meta, conserved,
        native_side_split=True,
    )


def _base(**kw) -> ConvergenceResult:
    d = dict(gene="G", position=42, tag="t", caas="A/B", ancestral="A",
             derived="V", convergence_type="convergent", change_top="parallel",
             change_bottom="parallel", change_side="both")
    d.update(kw)
    return ConvergenceResult(**d)


def test_both_position_expands_to_two_side_rows():
    rows = _split_result_by_side(_base(), _split_dict())
    assert [r.side for r in rows] == ["top", "bottom"]
    top, bottom = rows
    assert abs(top.asr_path_score - 0.5158540) < 1e-6
    assert abs(bottom.asr_path_score - 0.5430045) < 1e-6
    # sides never recombine
    assert top.asr_path_score != bottom.asr_path_score
    # per-side directional pair-score maps are one-sided
    assert top.pair_bottom_path_scores is None
    assert bottom.pair_top_path_scores is None
    # block independence retired
    assert top.independence is None


def test_one_sided_position_stays_single_row():
    # drop the bottom-only and conserved pairs -> only the top side participates
    pairs = [_BS_PAIRS[0] | {"bottom_tip_mode": "A"}, _BS_PAIRS[1]]
    rows = _split_result_by_side(_base(change_side="top"),
                                 _split_dict(pairs, conserved="", meta=False))
    assert len(rows) == 1 and rows[0].side == "top"
    assert abs(rows[0].asr_path_score - 0.7737809) < 1e-6


def test_no_participating_side_collapses_to_one_none_row():
    pairs = [{"pair_id": 1, "node_id": 3, "focal_state": "A",
              "top_tip_mode": "A", "bottom_tip_mode": "A"}]
    rows = _split_result_by_side(_base(change_side="none"),
                                 _split_dict(pairs, conserved="", meta=False))
    assert len(rows) == 1
    assert rows[0].side == "none" and rows[0].asr_path_score == 0.0


def test_flat_dict_carries_per_side_side_and_score():
    rows = _split_result_by_side(_base(), _split_dict())
    dicts = [convert_convergence_result_to_dict(r, None) for r in rows]
    assert {d["side"] for d in dicts} == {"top", "bottom"}
    by_side = {d["side"]: d for d in dicts}
    assert abs(by_side["top"]["asr_path_score"] - 0.5158540) < 1e-6
    assert abs(by_side["bottom"]["asr_path_score"] - 0.5430045) < 1e-6
    # same msa_pos on both rows -> DB stores two rows per (gene, msa_pos)
    assert by_side["top"]["msa_pos"] == by_side["bottom"]["msa_pos"] == 42


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
