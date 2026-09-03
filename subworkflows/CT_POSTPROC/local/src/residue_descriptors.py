#!/usr/bin/env python3
"""Position-level raw-AA descriptors for the disambiguation table.

Computed in CT_POSTPROC's input-prep step, upstream of the
``filtered_discovery.tsv`` fork, so SCORING and VEP consume ONE canonical column
set instead of each re-deriving the ancestral / derived residue sets from the
``caas`` pattern string.

Data model
----------
The ``caas`` string is ``<top>/<bottom>`` and is **positional**: index *i* of the
top string and index *i* of the bottom string are the two contrasting-phenotype
members of the same disjoint species pair *i*. A pair whose two members carry the
same residue is *conserved* and does not contribute to the derived residues
(``conserved_<j>_cons`` is a numeric conservation score for pair *j*, not a
residue). The per-pair reconstructed states live in
``mrca_<i>_{anc,top,bot}_aa`` (empty on a side that did not substitute).

``change_side`` (``top`` / ``bottom`` / ``both`` / ``none``) is the
disambiguation's authoritative call for which clade carries the substantive
change. **Side assignment follows ``change_side``, never a heuristic.**

Columns produced (one value per ``(Gene, Position)``, broadcast to every row):

``derived_residues``
    ``"<top>/<bottom>"`` — same left/right convention as ``caas``. The
    ``change_side``-sanctioned side shows its **derived** residues
    (``mrca_*_top_aa`` / ``mrca_*_bot_aa`` at changed pairs); the other side shows
    the **ancestral** residue (``mrca_*_anc_aa``). ``change_side == "both"`` shows
    the derived residues on both sides. ``""`` when the position has no changed
    pair.
``top_residue_support`` / ``bottom_residue_support``
    ``"L:3,S:2"`` — per residue listed on that side, the count of DISTINCT MRCA
    nodes (pairs) supporting it, count-descending then alphabetical. ``""`` when
    that side has no residue.
``n_conserved_pairs``
    Count of DISTINCT ``conserved_<j>_node`` values across the position's rows
    (``""`` / ``0`` when the conserved-pair block is absent).

``fop_pool.R`` carries these four columns straight through; it still computes only
the scheme-dependent ``convergence_schemes`` itself.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Set, Tuple

import pandas as pd

DESCRIPTOR_COLUMNS = (
    "derived_residues",
    "top_residue_support",
    "bottom_residue_support",
    "n_conserved_pairs",
)

_EMPTY = {"derived_residues": "", "top_residue_support": "",
          "bottom_residue_support": "", "n_conserved_pairs": ""}

# change_side -> which side(s) show DERIVED residues (the other shows ancestral).
_DERIVED_SIDES: Dict[str, Set[str]] = {
    "top": {"top"},
    "bottom": {"bot"},
    "both": {"top", "bot"},
    "none": set(),
    "": set(),
}


def _pair_indices(columns) -> List[int]:
    idx = set()
    for col in columns:
        m = re.fullmatch(r"mrca_(\d+)_node", str(col))
        if m:
            idx.add(int(m.group(1)))
    return sorted(idx)


def _conserved_indices(columns) -> List[int]:
    idx = set()
    for col in columns:
        m = re.fullmatch(r"conserved_(\d+)_node", str(col))
        if m:
            idx.add(int(m.group(1)))
    return sorted(idx)


def _clean_aa(val) -> str:
    if val is None:
        return ""
    s = str(val).strip().upper()
    if s in ("", "NA", "NAN", "NONE"):
        return ""
    return s


def _derived_sides_for(change_side) -> Set[str]:
    return _DERIVED_SIDES.get(str(change_side or "").strip().lower(), set())


def _node_str(val) -> str:
    s = "" if val is None else str(val).strip()
    return "" if s in ("", "nan", "NA", "None") else s


def _collect(group: pd.DataFrame, pair_idx: List[int]):
    """Per side, ``{residue: {nodes}}`` for derived changes and for ancestral cells.

    A ``mrca_<i>_top_aa`` / ``mrca_<i>_bot_aa`` cell is a *derived* change on that
    side; ``mrca_<i>_anc_aa`` is the ancestral residue of pair *i*. Both are keyed
    by residue with the set of distinct ``mrca_<i>_node`` values supporting them.
    """
    # Dedup by (node, side): a physical pair has one derived residue on a side
    # (first non-empty wins, matching fop_pool.R::.collect_changed_pairs). Same
    # for the ancestral cell, keyed by node.
    d_seen: Dict[Tuple[str, str], str] = {}
    a_seen: Dict[str, str] = {}
    for i in pair_idx:
        ncol = f"mrca_{i}_node"
        if ncol not in group.columns:
            continue
        nodes = list(group[ncol])
        for side, acol in (("top", f"mrca_{i}_top_aa"), ("bot", f"mrca_{i}_bot_aa")):
            if acol not in group.columns:
                continue
            for node_val, aa_val in zip(nodes, group[acol]):
                node = _node_str(node_val)
                aa = _clean_aa(aa_val)
                if node and aa:
                    d_seen.setdefault((node, side), aa)
        acol = f"mrca_{i}_anc_aa"
        if acol in group.columns:
            for node_val, aa_val in zip(nodes, group[acol]):
                node = _node_str(node_val)
                aa = _clean_aa(aa_val)
                if node and aa:
                    a_seen.setdefault(node, aa)

    derived: Dict[str, Dict[str, Set[str]]] = {"top": {}, "bot": {}}
    for (node, side), aa in d_seen.items():
        derived[side].setdefault(aa, set()).add(node)
    ancestral: Dict[str, Set[str]] = {}
    for node, aa in a_seen.items():
        ancestral.setdefault(aa, set()).add(node)
    return derived, ancestral


def _fmt_support(res_nodes: Dict[str, Set[str]]) -> str:
    counts = {aa: len(nodes) for aa, nodes in res_nodes.items() if nodes}
    if not counts:
        return ""
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{aa}:{n}" for aa, n in ordered)


def _n_conserved(group: pd.DataFrame, cons_idx: List[int]) -> str:
    if not cons_idx:
        return ""
    nodes: Set[str] = set()
    for j in cons_idx:
        ncol = f"conserved_{j}_node"
        if ncol not in group.columns:
            continue
        for v in group[ncol]:
            n = _node_str(v)
            if n:
                nodes.add(n)
    return str(len(nodes))


def _descriptors_for_group(group: pd.DataFrame, pair_idx: List[int],
                           cons_idx: List[int]) -> Dict[str, str]:
    derived, ancestral = _collect(group, pair_idx)
    n_cons = _n_conserved(group, cons_idx)

    changed = bool(derived["top"]) or bool(derived["bot"])
    if not changed:
        out = dict(_EMPTY)
        out["n_conserved_pairs"] = n_cons
        return out

    # change_side is a per-position call; every row of the group agrees.
    cside = ""
    if "change_side" in group.columns:
        vals = [str(v).strip().lower() for v in group["change_side"] if _node_str(v)]
        if vals:
            cside = vals[0]
    if cside in _DERIVED_SIDES and cside not in ("", "none"):
        der_sides = _DERIVED_SIDES[cside]
    else:
        # No usable change_side -> infer from which sides actually substituted.
        der_sides = {s for s in ("top", "bot") if derived[s]}

    # Ancestral residues: only those attached to a node that actually carried a
    # change (keeps the descriptor tied to the CAAS pairs, not the whole column).
    changed_nodes: Set[str] = set()
    for side in ("top", "bot"):
        for nodes in derived[side].values():
            changed_nodes |= nodes
    anc_at_change = {
        aa: (nodes & changed_nodes) for aa, nodes in ancestral.items()
    }
    anc_at_change = {aa: nodes for aa, nodes in anc_at_change.items() if nodes}

    def side_map(side: str) -> Dict[str, Set[str]]:
        # Sanctioned side -> its derived residues; the other side -> ancestral.
        return derived[side] if side in der_sides else anc_at_change

    top_map = side_map("top")
    bot_map = side_map("bot")

    top_field = "".join(sorted(top_map)) or "?"
    bot_field = "".join(sorted(bot_map)) or "?"
    return {
        "derived_residues": f"{top_field}/{bot_field}",
        "top_residue_support": _fmt_support(top_map),
        "bottom_residue_support": _fmt_support(bot_map),
        "n_conserved_pairs": n_cons,
    }


def add_residue_descriptors(
    df: pd.DataFrame,
    gene_col: str = "Gene",
    position_col: str = "Position",
) -> pd.DataFrame:
    """Return ``df`` with ``DESCRIPTOR_COLUMNS`` added (or overwritten).

    No-op-safe: if the raw ``mrca_<i>_node`` / ``mrca_<i>_*_aa`` block is absent,
    every row gets empty strings so the output schema is stable.
    """
    out = df.copy()
    pair_idx = _pair_indices(out.columns)
    cons_idx = _conserved_indices(out.columns)
    have_block = pair_idx and any(
        f"mrca_{i}_top_aa" in out.columns or f"mrca_{i}_bot_aa" in out.columns
        for i in pair_idx
    )
    if not have_block or gene_col not in out.columns or position_col not in out.columns:
        for c in DESCRIPTOR_COLUMNS:
            out[c] = ""
        return out

    per_group = {
        key: _descriptors_for_group(g, pair_idx, cons_idx)
        for key, g in out.groupby([gene_col, position_col], sort=False)
    }
    keys = list(zip(out[gene_col], out[position_col]))
    for c in DESCRIPTOR_COLUMNS:
        out[c] = [per_group[k][c] for k in keys]
    return out
