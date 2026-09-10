#!/usr/bin/env python3
"""Shared helpers for the T0 golden net (see docs roadmap, tier T0).

Both the fixture generator (``gen_golden.py``) and the golden tests
(``test_path_scores_golden.py``) import from here so the tree/posterior
(de)serialisation is defined exactly once. No production code depends on this
module.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

# ── module loading ───────────────────────────────────────────────────────────
_LOCAL = Path(__file__).resolve().parents[1]          # .../convergence
_SRC_ROOT = Path(__file__).resolve().parents[3]       # .../local  (has src/)
FIXTURE_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if str(_SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(_SRC_ROOT))
    spec = importlib.util.spec_from_file_location(name, _LOCAL / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── minimal tree node (same shape path_scores expects) ───────────────────────
class Node:
    __slots__ = ("node_id", "children", "parent")

    def __init__(self, nid: int):
        self.node_id = nid
        self.children: List["Node"] = []
        self.parent = None


def build_tree(edges: List[Tuple[int, int]], node_ids: List[int]):
    """edges: list of (parent_id, child_id). Returns (root, {id: Node})."""
    nodes: Dict[int, Node] = {i: Node(i) for i in node_ids}
    child_ids = set()
    for p, c in edges:
        nodes[c].parent = nodes[p]
        nodes[p].children.append(nodes[c])
        child_ids.add(c)
    root = next(nodes[i] for i in node_ids if i not in child_ids)
    return root, nodes


# ── numeric comparison ───────────────────────────────────────────────────────
def diff_report(expected: Any, got: Any, tol: float, path: str = "") -> List[str]:
    """Recursively compare nested dict/list/number structures."""
    out: List[str] = []
    if isinstance(expected, dict):
        if not isinstance(got, dict):
            return [f"{path}: type {type(got).__name__} != dict"]
        for k in sorted(set(expected) | set(got), key=str):
            if k not in expected:
                out.append(f"{path}.{k}: unexpected key")
            elif k not in got:
                out.append(f"{path}.{k}: missing key")
            else:
                out += diff_report(expected[k], got[k], tol, f"{path}.{k}")
    elif isinstance(expected, list):
        if not isinstance(got, list) or len(got) != len(expected):
            return [f"{path}: list mismatch ({got!r} != {expected!r})"]
        for i, (a, b) in enumerate(zip(expected, got)):
            out += diff_report(a, b, tol, f"{path}[{i}]")
    elif isinstance(expected, (int, float)) and isinstance(got, (int, float)):
        if abs(float(expected) - float(got)) > tol:
            out.append(f"{path}: {got!r} != {expected!r} (|d|>{tol})")
    else:
        if expected != got:
            out.append(f"{path}: {got!r} != {expected!r}")
    return out


def load_json(name: str) -> Any:
    return json.loads((FIXTURE_DIR / name).read_text())
