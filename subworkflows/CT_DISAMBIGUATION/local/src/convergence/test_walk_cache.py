#!/usr/bin/env python3
"""walk_cache equivalence: compute_asr_path_score must be bit-identical with and
without the memoisation, and the cache must actually be reused across repeated
scorings of the same site (the FOP-null / multi-hypothesis case).

Run: python test_walk_cache.py
"""
import importlib.util
import sys
from pathlib import Path

# path_scores imports `from src.biochem.grouping import get_grouping_scheme`
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # .../local (has src/)
_spec = importlib.util.spec_from_file_location(
    "path_scores", Path(__file__).with_name("path_scores.py"))
ps = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ps)


class N:
    __slots__ = ("node_id", "children", "parent")
    def __init__(self, nid, parent=None):
        self.node_id = nid
        self.parent = parent
        self.children = []


def build_tree():
    # root(0) -> a(1) -> b(2) -> {leaf c(3), sub d(4) -> {e(5), f(6)}}
    #                 -> g(7) -> {h(8), i(9)}
    nodes = {i: N(i) for i in range(10)}
    def link(p, c):
        nodes[c].parent = nodes[p]; nodes[p].children.append(nodes[c])
    link(0, 1); link(1, 2); link(1, 7)
    link(2, 3); link(2, 4); link(4, 5); link(4, 6)
    link(7, 8); link(7, 9)
    return nodes[0], nodes


def posteriors():
    # {node_id: {aa: prob}} for one site. US scheme -> raw residues.
    return {
        0: {"A": 0.9, "V": 0.1},
        1: {"A": 0.8, "V": 0.2},
        2: {"A": 0.7, "V": 0.3},
        4: {"A": 0.5, "V": 0.5},
        5: {"V": 0.9, "A": 0.1},
        6: {"V": 0.85, "A": 0.15},
        7: {"A": 0.6, "T": 0.4},
        8: {"T": 0.8, "A": 0.2},
        9: {"T": 0.75, "A": 0.25},
    }


def pair_details():
    # Two changed pairs converging on V from ancestral A, one on a different clade.
    return [
        {"pair_id": 1, "node_id": 4, "focal_state": "A",
         "top_tip_mode": "V", "bottom_tip_mode": "A"},
        {"pair_id": 2, "node_id": 7, "focal_state": "A",
         "top_tip_mode": "T", "bottom_tip_mode": "A"},
        {"pair_id": 3, "node_id": 2, "focal_state": "A",
         "top_tip_mode": "V", "bottom_tip_mode": "A"},
    ]


def run():
    root, nodes = build_tree()
    node_index = ps.build_node_index(root)
    pnd = posteriors()
    pd = pair_details()

    ok = True
    for scheme in ("US", "GS1", "GS4"):
        base = ps.compute_asr_path_score(pd, pnd, node_index, scheme, False, "")
        cache = {}
        c1 = ps.compute_asr_path_score(pd, pnd, node_index, scheme, False, "",
                                       walk_cache=cache, site_key=42)
        n_after_first = len(cache)
        # score it again (same site) — must not grow the cache and must match
        c2 = ps.compute_asr_path_score(pd, pnd, node_index, scheme, False, "",
                                       walk_cache=cache, site_key=42)
        for k in ("asr_path_score", "core", "independence", "mrca_diversity",
                  "derived_agreement", "conservation_gate"):
            a, b = base[k], c1[k]
            same = (a == b) or (abs(a - b) < 1e-12)
            print(f"  {scheme:4s} {k:18s} uncached={a:.10f} cached={b:.10f} "
                  f"{'OK' if same else 'MISMATCH'}")
            ok &= same
        same_pairs = base["pair_scores"] == c1["pair_scores"]
        print(f"  {scheme:4s} pair_scores identical: {same_pairs}")
        ok &= same_pairs
        # POINT 3: raw derived/ancestral residue maps are returned and unaffected
        # by the walk memoisation.
        raw_ok = all(
            base.get(k) == c1.get(k) and isinstance(base.get(k), dict)
            for k in ("pair_ancestral", "pair_derived_top", "pair_derived_bot")
        )
        if scheme == "US":
            # pair 1 (node 4) changed top A->V; both pairs' bottom sides conserved
            raw_ok &= base["pair_derived_top"].get(1) == "V" and base["pair_ancestral"].get(1) == "A"
            raw_ok &= base["pair_derived_bot"] == {}
        print(f"  {scheme:4s} raw residue maps returned + cache-invariant: {raw_ok}")
        ok &= raw_ok
        reused = (len(cache) == n_after_first) and (c2["asr_path_score"] == c1["asr_path_score"])
        print(f"  {scheme:4s} cache reused on 2nd call (no growth): {reused} "
              f"({n_after_first} entries)")
        ok &= reused

    # conserved-pair path exercises the is_changed=False cache branch
    pd_cons = pd + []
    base_c = ps.compute_asr_path_score(pd, pnd, node_index, "US", True, "3")
    cache = {}
    cc = ps.compute_asr_path_score(pd, pnd, node_index, "US", True, "3",
                                   walk_cache=cache, site_key=7)
    same = abs(base_c["conservation_gate"] - cc["conservation_gate"]) < 1e-12
    print(f"  conserved-pair conservation_gate: uncached={base_c['conservation_gate']:.10f} "
          f"cached={cc['conservation_gate']:.10f} {'OK' if same else 'MISMATCH'}")
    ok &= same

    # POINT 2: conserved-pair identity maps are returned and non-empty when the
    # position has conserved pairs (pid "3" here), both early and normal returns.
    cps = base_c.get("conserved_pair_scores")
    cpn = base_c.get("conserved_pair_nodes")
    have_maps = isinstance(cps, dict) and isinstance(cpn, dict) and cps and \
        set(cps) == set(cpn)
    print(f"  conserved_pair_scores/nodes returned: {have_maps} ({cps})")
    ok &= bool(have_maps)
    same_maps = (cps == cc.get("conserved_pair_scores")
                 and cpn == cc.get("conserved_pair_nodes"))
    ok &= same_maps
    # A position with no conserved pairs -> empty maps present in the return.
    base_nc = ps.compute_asr_path_score(pd, pnd, node_index, "US", False, "")
    ok &= base_nc.get("conserved_pair_scores") == {} and \
        base_nc.get("conserved_pair_nodes") == {}

    print("\n" + ("ALL EQUIVALENCE CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(run())
