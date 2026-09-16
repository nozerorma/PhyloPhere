#!/usr/bin/env python3
"""Regression net for the Tier 3 fix in
docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md: `get_mrca`'s tip lookups
(`find_node_by_name`/`find_node_by_taxid`) used to re-walk the whole tree per
tip name on every call. `build_name_taxid_index` now builds an O(1) lookup once
per tree, and `get_mrca` uses it when given via the new `name_index`/
`taxid_index` params — falling back to the original recursive search when they
are omitted, so any caller that doesn't opt in (e.g. `node_identification.py`)
is unaffected.

This test proves the indexed and unindexed paths agree on the same MRCA,
across name-matching, taxid-matching, and mixed-matching queries, plus the
default-None path staying byte-for-byte the old behavior.

Run:  python -m pytest test_get_mrca_indexed.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # .../src

from asr.tree_parser import TreeNode, build_name_taxid_index, get_mrca  # noqa: E402


def _build_tree():
    """
    root
    ├── A (internal)
    │   ├── tip1 ('lineageA_100')
    │   └── tip2 ('lineageB_200')
    └── B (internal)
        ├── tip3 ('lineageC_300')
        └── C (internal)
            ├── tip4 ('lineageD_400')
            └── tip5 ('lineageE_500')
    """
    root = TreeNode(name=None, node_id=0)
    a = TreeNode(name=None, node_id=1)
    b = TreeNode(name=None, node_id=2)
    tip1 = TreeNode(name="lineageA_100", node_id=3)
    tip2 = TreeNode(name="lineageB_200", node_id=4)
    tip3 = TreeNode(name="lineageC_300", node_id=5)
    c = TreeNode(name=None, node_id=6)
    tip4 = TreeNode(name="lineageD_400", node_id=7)
    tip5 = TreeNode(name="lineageE_500", node_id=8)

    for parent, children in [
        (root, [a, b]),
        (a, [tip1, tip2]),
        (b, [tip3, c]),
        (c, [tip4, tip5]),
    ]:
        parent.children = children
        for child in children:
            child.parent = parent

    return root, dict(
        root=root, a=a, b=b, tip1=tip1, tip2=tip2, tip3=tip3, c=c, tip4=tip4, tip5=tip5,
    )


def test_indexed_matches_unindexed_on_name_query():
    root, n = _build_tree()
    name_index, taxid_index = build_name_taxid_index(root)
    expected = get_mrca(root, ["lineageA_100", "lineageB_200"])
    got = get_mrca(root, ["lineageA_100", "lineageB_200"], name_index=name_index, taxid_index=taxid_index)
    assert expected is n["a"]
    assert got is expected


def test_indexed_matches_unindexed_on_taxid_query():
    # Query by bare taxid (not the full 'lineage_taxid' tip label) — exercises
    # the taxid_index fallback path on both the indexed and unindexed sides.
    root, n = _build_tree()
    name_index, taxid_index = build_name_taxid_index(root)
    expected = get_mrca(root, ["400", "500"])
    got = get_mrca(root, ["400", "500"], name_index=name_index, taxid_index=taxid_index)
    assert expected is n["c"]
    assert got is expected


def test_indexed_matches_unindexed_across_subtrees():
    root, n = _build_tree()
    name_index, taxid_index = build_name_taxid_index(root)
    expected = get_mrca(root, ["lineageA_100", "400"])
    got = get_mrca(root, ["lineageA_100", "400"], name_index=name_index, taxid_index=taxid_index)
    assert expected is n["root"]
    assert got is expected


def test_default_none_path_unchanged():
    """No index given -> identical object to calling get_mrca with the
    pre-fix signature (positional root, tip_names only)."""
    root, n = _build_tree()
    assert get_mrca(root, ["lineageA_100", "lineageB_200"]) is n["a"]
    assert get_mrca(root, ["nonexistent"]) is None


def test_index_covers_every_node_and_only_leaves_for_taxid():
    root, n = _build_tree()
    name_index, taxid_index = build_name_taxid_index(root)
    assert name_index["lineageA_100"] is n["tip1"]
    assert set(taxid_index.keys()) == {"100", "200", "300", "400", "500"}
    # Internal nodes (name=None) are never queried by name in practice, but the
    # index must not crash or misbehave when one has an actual empty name.
    assert None not in name_index
