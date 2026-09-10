#!/usr/bin/env python3
"""Tests for fop_pool.pool_hypotheses_pairwise — the FOP null's real per-side
pairwise core (scoring_v2 T3c SC2b).

Oracle: for a harvest of ONE hypothesis, pooling is identity — the pooled
``core_s`` must equal ``compute_asr_path_score(native_side_split=True)`` for that
hypothesis's own pairs. Multi-hypothesis cases are checked against
``aggregate_core_side`` run by hand over the node-deduped union.

Run: python -m pytest test_fop_pool_pairwise.py
"""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # .../local (has src/)

_HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ps = _load("path_scores")
fp = _load("fop_pool")


# ── shared tiny tree (golden both_sides_two_rows, T3-doc §10.1) ───────────────
EDGES = [(0, 1), (0, 2), (1, 10), (10, 3), (10, 4), (2, 20), (20, 5), (20, 6),
         (2, 21), (21, 9), (21, 12), (0, 30), (30, 7), (30, 8)]
NODE_IDS = [0, 1, 2, 10, 20, 21, 30, 3, 4, 5, 6, 9, 12, 7, 8]
POST = {n: {"A": 0.90, "V": 0.05, "L": 0.05}
        for n in (0, 1, 2, 10, 20, 21, 30)}


class _Node:
    __slots__ = ("node_id", "children", "parent")

    def __init__(self, nid):
        self.node_id, self.children, self.parent = nid, [], None


def _tree():
    nodes = {i: _Node(i) for i in NODE_IDS}
    kids = set()
    for p, c in EDGES:
        nodes[c].parent = nodes[p]
        nodes[p].children.append(nodes[c])
        kids.add(c)
    root = next(nodes[i] for i in NODE_IDS if i not in kids)
    return ps.build_node_index(root)


NODE_INDEX = _tree()
PND = {int(k): dict(sorted(v.items())) for k, v in POST.items()}

PAIRS = [
    {"pair_id": 1, "node_id": 3, "focal_state": "A", "top_tip_mode": "V", "bottom_tip_mode": "L"},
    {"pair_id": 2, "node_id": 5, "focal_state": "A", "top_tip_mode": "V", "bottom_tip_mode": "A"},
    {"pair_id": 3, "node_id": 7, "focal_state": "A", "top_tip_mode": "A", "bottom_tip_mode": "L"},
    {"pair_id": 4, "node_id": 9, "focal_state": "A", "top_tip_mode": "A", "bottom_tip_mode": "A"},
]


def approx(a, b, tol=1e-9):
    return a is not None and b is not None and abs(a - b) < tol


def _rec_from_split(hyp, pairs, conserved_pair="", is_conserved_meta=False):
    split = ps.compute_asr_path_score(
        pairs, PND, NODE_INDEX, "US", is_conserved_meta, conserved_pair,
    )
    return {"hyp": hyp, "sides": split}, split


def test_single_hypothesis_is_identity():
    rec, split = _rec_from_split("H1", PAIRS, conserved_pair="4", is_conserved_meta=True)
    out = fp.pool_hypotheses_pairwise([rec], NODE_INDEX, PND, {}, "US")
    assert out["n_hypotheses"] == 1
    # identity: one-hyp harvest == compute_asr_path_score for that hyp's pairs
    assert approx(out["top"]["asr_path_score"], split["top"]["asr_path_score"])
    assert approx(out["bottom"]["asr_path_score"], split["bottom"]["asr_path_score"])
    # golden §10.1: core_top ≈ 0.51585, core_bottom ≈ 0.54301
    assert approx(out["top"]["asr_path_score"], 0.51585, tol=1e-4)
    assert approx(out["bottom"]["asr_path_score"], 0.54301, tol=1e-4)
    assert out["top"]["n_pairs_side"] == 3      # P1, P2, + conserved P4
    assert out["bottom"]["n_pairs_side"] == 3   # P1, P3, + conserved P4


def test_two_hypotheses_shared_node_pool_iso():
    # H1 and H2 both carry P1 (node 3) and P2 (node 5) on top, same residues.
    # No PSS -> equal weight -> pooled iso == the (identical) per-hyp iso, so the
    # pooled top core equals the single-hyp top core exactly.
    rec1, split = _rec_from_split("H1", PAIRS[:2])
    rec2, _ = _rec_from_split("H2", PAIRS[:2])
    out = fp.pool_hypotheses_pairwise([rec1, rec2], NODE_INDEX, PND, {}, "US")
    assert out["n_hypotheses"] == 2
    assert approx(out["top"]["asr_path_score"], split["top"]["asr_path_score"])
    assert out["top"]["n_participating"] == 2


def test_two_hypotheses_disjoint_nodes_union():
    # H1 carries only P1 (node 3, ->V top); H2 only P2 (node 5, ->V top). Neither
    # hypothesis alone has a same-residue partner (core 0). The UNION pairs them
    # -> core > 0. It does NOT match the both-in-one-hyp core: P1's isolation was
    # measured in H1 with no sibling merge point, so its private walk runs one
    # node deeper (0.857 vs 0.902) -> a slightly lower pooled core. That is the
    # intended semantics: the harvest never co-tested P1 and P2.
    rec1, _ = _rec_from_split("H1", [PAIRS[0]])
    rec2, _ = _rec_from_split("H2", [PAIRS[1]])
    solo = fp.pool_hypotheses_pairwise([rec1], NODE_INDEX, PND, {}, "US")
    assert approx(solo["top"]["asr_path_score"], 0.0)
    union = fp.pool_hypotheses_pairwise([rec1, rec2], NODE_INDEX, PND, {}, "US")
    _, split = _rec_from_split("Hx", PAIRS[:2])
    assert 0.0 < union["top"]["asr_path_score"] < split["top"]["asr_path_score"]
    assert union["top"]["n_participating"] == 2


def test_pss_weights_pooled_iso():
    # Same shared node across two hyps but different iso would need different
    # posteriors; here iso is identical so PSS cannot move the mean — assert the
    # weighting path at least runs and stays identity.
    rec1, split = _rec_from_split("H1", PAIRS[:2])
    rec2, _ = _rec_from_split("H2", PAIRS[:2])
    pss = {("H1", 1): 3.0, ("H2", 1): 1.0, ("H1", 2): 2.0, ("H2", 2): 2.0}
    out = fp.pool_hypotheses_pairwise([rec1, rec2], NODE_INDEX, PND, pss, "US")
    assert approx(out["top"]["asr_path_score"], split["top"]["asr_path_score"])


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    ok = True
    for fn in fns:
        try:
            fn()
            print("PASS ", fn.__name__)
        except AssertionError as e:
            ok = False
            print("FAIL ", fn.__name__, "--", e)
    print("\n" + ("ALL TESTS PASSED" if ok else "SOME TESTS FAILED"))
    sys.exit(0 if ok else 1)
