#!/usr/bin/env python3
"""T3c SC3: in-tree FOP pooling of the observed hypothesis harvest.

disambiguate_single._pool_observed_fop groups per-hypothesis ConvergenceResult
rows by (position, scheme) and, for a ≥2-hypothesis group, collapses them with
fop_pool.pool_hypotheses_pairwise (the same routine the null uses). A
1-hypothesis / non-FOP group passes through untouched.

Run: python -m pytest test_pool_observed_fop.py
"""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
_HERE = Path(__file__).resolve().parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ps = _load("path_scores")
ds = _load("disambiguate_single")
from src.data.models import ConvergenceResult  # noqa: E402

EDGES = [(0, 1), (0, 2), (1, 10), (10, 3), (10, 4), (2, 20), (20, 5), (20, 6),
         (2, 21), (21, 9), (21, 12), (0, 30), (30, 7), (30, 8)]
NODE_IDS = [0, 1, 2, 10, 20, 21, 30, 3, 4, 5, 6, 9, 12, 7, 8]
POST = {n: {"A": 0.90, "V": 0.05, "L": 0.05} for n in (0, 1, 2, 10, 20, 21, 30)}


class _N:
    __slots__ = ("node_id", "children", "parent")

    def __init__(self, i):
        self.node_id, self.children, self.parent = i, [], None


class _Tree:
    def __init__(self):
        nodes = {i: _N(i) for i in NODE_IDS}
        kids = set()
        for p, c in EDGES:
            nodes[c].parent = nodes[p]
            nodes[p].children.append(nodes[c])
            kids.add(c)
        self.root = next(nodes[i] for i in NODE_IDS if i not in kids)


TREE = _Tree()
POSTERIOR_DATA = {nid: {6: dict(sorted(d.items()))} for nid, d in POST.items()}  # site 6 (pos0=5)

PAIRS = [
    {"pair_id": 1, "node_id": 3, "focal_state": "A", "top_tip_mode": "V", "bottom_tip_mode": "L"},
    {"pair_id": 2, "node_id": 5, "focal_state": "A", "top_tip_mode": "V", "bottom_tip_mode": "A"},
]


def _split(pairs):
    pnd = {nid: POSTERIOR_DATA[nid][6] for nid in POSTERIOR_DATA}
    ni = ps.build_node_index(TREE.root)
    return ps.compute_asr_path_score(pairs, pnd, ni, "US", False, "")


def _cr(hyp, side, split, pairs):
    sd = split.get(side) or {}
    return ConvergenceResult(
        gene="G", position=5, tag="POS5", caas="V/A", ancestral="A", derived="V",
        convergence_type=sd.get("convergence_type", "single"),
        caap_group="US", hypothesis=hyp, position_one_based=6,
        pair_details=pairs, side=side,
        asr_path_score=sd.get("asr_path_score", 0.0), core=sd.get("core", 0.0),
    )


def _rows_for_hyp(hyp, split, pairs):
    out = []
    for s in ("top", "bottom"):
        if int((split.get(s) or {}).get("n_participating", 0) or 0) > 0:
            r = _cr(hyp, s, split, pairs)
            r.path_split = split
            out.append(r)
    return out


def test_single_hypothesis_passthrough():
    sp = _split(PAIRS)
    rows = _rows_for_hyp("H1", sp, PAIRS)
    out = ds._pool_observed_fop(rows, TREE, POSTERIOR_DATA)
    assert out == rows  # untouched


def test_two_hypotheses_pooled_identity():
    # Both hypotheses carry the identical pair set -> union dedups by node, iso
    # pools over identical values -> pooled core == compute_asr_path_score.
    sp = _split(PAIRS)
    rows = _rows_for_hyp("H1", sp, PAIRS) + _rows_for_hyp("H2", sp, PAIRS)
    assert len(rows) == 4
    out = ds._pool_observed_fop(rows, TREE, POSTERIOR_DATA)
    by_side = {r.side: r for r in out}
    # top: both pairs -> V, converge. bottom: only pair 1 -> L, no partner.
    assert set(by_side) == {"top", "bottom"}
    assert all(r.hypothesis is None for r in out)
    assert abs(by_side["top"].asr_path_score - sp["top"]["asr_path_score"]) < 1e-9
    assert by_side["top"].asr_path_score > 0.5
    assert abs(by_side["bottom"].asr_path_score - sp["bottom"]["asr_path_score"]) < 1e-9
    assert by_side["bottom"].asr_path_score == 0.0


def test_non_fop_untouched():
    sp = _split(PAIRS)
    rows = _rows_for_hyp(None, sp, PAIRS)  # hypothesis None -> single contrast
    out = ds._pool_observed_fop(rows, TREE, POSTERIOR_DATA)
    assert out == rows


def test_pss_weights_threaded():
    # H1 carries pair 1 (node 3), H2 carries pair 2 (node 5). Distinct nodes ->
    # the union has both. PSS is passed straight through to
    # pool_hypotheses_pairwise; with identical iso the weighting cannot move the
    # result, so assert it at least runs and stays a valid pooled row.
    sp1 = _split([PAIRS[0]])
    sp2 = _split([PAIRS[1]])
    rows = _rows_for_hyp("H1", sp1, [PAIRS[0]]) + _rows_for_hyp("H2", sp2, [PAIRS[1]])
    pss = {("H1", 1): 3.0, ("H2", 2): 1.0}
    out = ds._pool_observed_fop(rows, TREE, POSTERIOR_DATA, pss)
    by_side = {r.side: r for r in out}
    assert "top" in by_side and by_side["top"].hypothesis is None
    assert by_side["top"].asr_path_score > 0.0


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
