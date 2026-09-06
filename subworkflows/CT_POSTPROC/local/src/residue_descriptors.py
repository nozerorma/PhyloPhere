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
    ``"L:3,S:2"`` — per residue listed on that side, the number of DISTINCT CAAS
    contrast pairs (``mrca_<i>`` blocks) that carry it, count-descending then
    alphabetical. This is the *actual* support the position has: it is bounded by
    the CAAS's pair count and is NOT inflated by the number of discovering
    hypotheses. ``""`` when that side has no residue.
``top_residue_support_detail`` / ``bottom_residue_support_detail``
    Same shape, but counting DISTINCT reconstructed ancestral nodes across every
    discovering hypothesis rather than physical pairs. A single pair can resolve
    to different ``mrca_<i>_node`` values under different hypotheses, so this
    number blends pair count with reconstruction/hypothesis multiplicity — kept
    as a secondary, finer-grained view, not an evidence count.
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
    "top_residue_support_detail",
    "bottom_residue_support_detail",
    "n_conserved_pairs",
)

_EMPTY = {"derived_residues": "", "top_residue_support": "",
          "bottom_residue_support": "", "top_residue_support_detail": "",
          "bottom_residue_support_detail": "", "n_conserved_pairs": ""}

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
    """Per side, ``{residue: {support units}}`` for derived changes and ancestral cells.

    A ``mrca_<i>_top_aa`` / ``mrca_<i>_bot_aa`` cell is a *derived* change on that
    side; ``mrca_<i>_anc_aa`` is the ancestral residue of pair *i*.

    Returns two (derived, ancestral) pairs of dicts:
      * ``*_pairs``: support unit = the CAAS contrast pair index *i* (physical
        support — bounded by the CAAS's pair count, hypothesis-invariant).
      * ``*_nodes``: support unit = the distinct reconstructed ``mrca_<i>_node``
        value (finer, but blends pair count with hypothesis multiplicity because
        one pair can reconstruct to different nodes under different hypotheses).
    """
    # (i, side) -> {derived residues seen for that pair across hypotheses};
    # (i) -> {ancestral residues seen}. A pair keeps its full residue set so
    # `derived_residues` stays complete under hypothesis disagreement, but each
    # residue counts that pair only ONCE (physical support, hypothesis-invariant).
    dp_seen: Dict[Tuple[int, str], Set[str]] = {}
    ap_seen: Dict[int, Set[str]] = {}
    # (node, side) -> derived residue; node -> ancestral residue.
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
                if not aa:
                    continue
                dp_seen.setdefault((i, side), set()).add(aa)
                if node:
                    d_seen.setdefault((node, side), aa)
        acol = f"mrca_{i}_anc_aa"
        if acol in group.columns:
            for node_val, aa_val in zip(nodes, group[acol]):
                node = _node_str(node_val)
                aa = _clean_aa(aa_val)
                if not aa:
                    continue
                ap_seen.setdefault(i, set()).add(aa)
                if node:
                    a_seen.setdefault(node, aa)

    derived_pairs: Dict[str, Dict[str, Set[int]]] = {"top": {}, "bot": {}}
    for (i, side), aas in dp_seen.items():
        for aa in aas:
            derived_pairs[side].setdefault(aa, set()).add(i)
    ancestral_pairs: Dict[str, Set[int]] = {}
    for i, aas in ap_seen.items():
        for aa in aas:
            ancestral_pairs.setdefault(aa, set()).add(i)

    derived_nodes: Dict[str, Dict[str, Set[str]]] = {"top": {}, "bot": {}}
    for (node, side), aa in d_seen.items():
        derived_nodes[side].setdefault(aa, set()).add(node)
    ancestral_nodes: Dict[str, Set[str]] = {}
    for node, aa in a_seen.items():
        ancestral_nodes.setdefault(aa, set()).add(node)

    return derived_pairs, ancestral_pairs, derived_nodes, ancestral_nodes


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
    derived_p, ancestral_p, derived_n, ancestral_n = _collect(group, pair_idx)
    n_cons = _n_conserved(group, cons_idx)

    changed = bool(derived_p["top"]) or bool(derived_p["bot"])
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
        der_sides = {s for s in ("top", "bot") if derived_p[s]}

    def _side_fields(derived, ancestral):
        # Ancestral residues restricted to the support units that actually
        # carried a change (keeps the descriptor tied to the CAAS pairs).
        changed_units: Set = set()
        for side in ("top", "bot"):
            for units in derived[side].values():
                changed_units |= units
        anc_at_change = {aa: (units & changed_units) for aa, units in ancestral.items()}
        anc_at_change = {aa: units for aa, units in anc_at_change.items() if units}

        def side_map(side):
            return derived[side] if side in der_sides else anc_at_change

        return side_map("top"), side_map("bot")

    top_p, bot_p = _side_fields(derived_p, ancestral_p)
    top_n, bot_n = _side_fields(derived_n, ancestral_n)

    top_field = "".join(sorted(top_p)) or "?"
    bot_field = "".join(sorted(bot_p)) or "?"
    return {
        "derived_residues": f"{top_field}/{bot_field}",
        "top_residue_support": _fmt_support(top_p),
        "bottom_residue_support": _fmt_support(bot_p),
        "top_residue_support_detail": _fmt_support(top_n),
        "bottom_residue_support_detail": _fmt_support(bot_n),
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
