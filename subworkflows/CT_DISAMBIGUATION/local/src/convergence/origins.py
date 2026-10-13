#!/usr/bin/env python3
# origins.py — Independent origins of a derived residue on the species tree.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/convergence/

"""
OriginCounter: how many times a residue arises independently among the species that carry it.

The residue is read as a character of the species tree: a tip is a carrier when it holds the residue at the
alignment column. An origin is a maximal clade whose tips (the ones with a residue at the column) all hold it,
so carriers that inherit the residue from a common ancestor are one origin, and an isolated carrier is one origin
by itself. A tip without a residue (gap or missing) does not decide: it neither breaks nor makes a clade.

Imported by: src/core/driver.py, src/convergence/disambiguate_single.py
Inputs:      the tree (src.asr.tree_parser.TreeNode), a function that gives the residue of a tip at a column
Outputs:     number of origins of a residue among a set of carrier tips
"""

from typing import Callable, Dict, Iterable, List, Optional, Tuple


class OriginCounter:
    """Counts the origins of residues at alignment columns on one fixed tree.

    `tip_residue(tip_label, column)` returns the residue of a tip (its `name` in the tree) at the column, or None.
    The clades of a (column, residue) are computed once and reused for every set of carriers asked about it.
    """

    def __init__(self, root, tip_residue: Callable[[str, int], Optional[str]]):
        self._tip_residue = tip_residue
        # post-order list of (node, children) so that no recursion is needed
        self._order: List = []
        self._parent: Dict[int, Optional[int]] = {}
        stack = [(root, None, False)]
        while stack:
            node, parent, seen = stack.pop()
            if seen:
                self._order.append(node)
                continue
            self._parent[id(node)] = id(parent) if parent is not None else None
            stack.append((node, parent, True))
            for child in node.children:
                stack.append((child, node, False))
        self._root = root
        self._clades: Dict[Tuple[int, str], Dict[str, int]] = {}

    def _tip_clades(self, column: int, residue: str) -> Dict[str, int]:
        """Clade id of every tip that sits in a maximal clade of carriers of `residue` at `column`."""
        key = (column, residue)
        if key in self._clades:
            return self._clades[key]
        # state of a node: True when the tips below it that have a residue all hold `residue` (and there is one),
        # False when one of them holds another residue, None when no tip below it has a residue
        state: Dict[int, Optional[bool]] = {}
        for node in self._order:
            if not node.children:
                res = self._tip_residue(node.name, column) if node.name else None
                state[id(node)] = None if res is None else (res == residue)
                continue
            kids = [state[id(c)] for c in node.children if state[id(c)] is not None]
            state[id(node)] = None if not kids else all(kids)
        # maximal clades: an all-carrier node whose parent is not an all-carrier node
        tip_clade: Dict[str, int] = {}
        next_id = 0
        stack = [self._root]
        while stack:
            node = stack.pop()
            parent_id = self._parent[id(node)]
            if state[id(node)] is True and (parent_id is None or state.get(parent_id) is not True):
                cid = next_id
                next_id += 1
                below = [node]
                while below:
                    n = below.pop()
                    if n.children:
                        below.extend(n.children)
                    elif n.name:
                        tip_clade[n.name] = cid
                continue            # everything below belongs to this clade
            stack.extend(node.children)
        self._clades[key] = tip_clade
        return tip_clade

    def count(self, column: int, residue: str, carriers: Iterable[str]) -> int:
        """Number of origins of `residue` at `column` among the carrier tips (tree tip labels)."""
        tip_clade = self._tip_clades(column, residue)
        return len({tip_clade[c] for c in carriers if c in tip_clade})


def build_origin_counter(root, taxid_index, seq_by_id, seq_by_species, species_to_taxid):
    """`origins(column, residue, species)`: origins of `residue` at the 0-based alignment `column` among the carrier
    species (names), on the tree `root`. `taxid_index` maps a taxid to its tip, `species_to_taxid` a species name to
    its taxid, and the sequences are the alignment lookups of the gene."""
    from src.convergence.convergence import _fetch_residue_from_lookup

    def tip_residue(label, column):
        taxid = label.strip().split("_")[-1]
        res = _fetch_residue_from_lookup(taxid, column, seq_by_id, seq_by_species)
        return res.upper() if res and res.isalpha() else None

    counter = OriginCounter(root, tip_residue)

    def origins(column, residue, species):
        labels = []
        for sp in species:
            node = taxid_index.get(str(species_to_taxid.get(sp, sp)))
            if node is not None and node.name:
                labels.append(node.name)
        return counter.count(column, residue, labels)

    return origins
