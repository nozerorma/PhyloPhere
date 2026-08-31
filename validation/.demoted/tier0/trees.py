"""Tree loading, depth-preserving pruning, pathological trees, foreground sampling.

Uses dendropy (a pipeline dependency). Everything downstream works on a light
``PhyloTree`` view: a preorder list of edges ``(parent_idx, child_idx, length)``
with node 0 the root, plus the tip labels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import dendropy
import numpy as np


@dataclass
class PhyloTree:
    labels: list[str]                       # node index -> label ("" for internal)
    edges: list[tuple[int, int, float]]     # (parent_idx, child_idx, branch_length), preorder
    root: int = 0

    @property
    def n_nodes(self) -> int:
        return len(self.labels)

    @property
    def tips(self) -> list[str]:
        return [self.labels[i] for i in self._leaves()]

    def _leaves(self) -> list[int]:
        has_child = {p for p, _, _ in self.edges}
        return [i for i in range(self.n_nodes) if i not in has_child]

    def total_length(self) -> float:
        return sum(bl for _, _, bl in self.edges)

    def children_of(self, node: int) -> list[int]:
        return [c for p, c, _ in self.edges if p == node]

    def descendant_tips(self, node: int) -> set[str]:
        leaves = set(self._leaves())
        out: set[str] = set()
        stack = [node]
        while stack:
            n = stack.pop()
            if n in leaves:
                out.add(self.labels[n])
            else:
                stack.extend(self.children_of(n))
        return out


def _from_dendropy(t: dendropy.Tree) -> PhyloTree:
    labels: list[str] = []
    edges: list[tuple[int, int, float]] = []
    idx: dict[int, int] = {}

    seed = t.seed_node
    idx[id(seed)] = 0
    labels.append(seed.taxon.label if seed.taxon else "")
    for node in t.preorder_node_iter():
        if node is seed:
            continue
        i = len(labels)
        idx[id(node)] = i
        labels.append(node.taxon.label if node.taxon else "")
        bl = node.edge.length if node.edge and node.edge.length is not None else 0.0
        edges.append((idx[id(node.parent_node)], i, float(bl)))
    return PhyloTree(labels=labels, edges=edges, root=0)


def load_tree(path: str | Path) -> PhyloTree:
    t = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    t.suppress_unifurcations()
    return _from_dendropy(t)


def prune_depth_preserving(path: str | Path, k: int, rng: np.random.Generator) -> PhyloTree:
    """Retain ``k`` tips chosen to keep the tree's depth range: always keep the
    two tips defining the maximum pairwise distance, then fill the rest by a
    farthest-point sweep so the subsample is spread across the tree rather than
    clustered in one dense clade.
    """
    t = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    t.suppress_unifurcations()
    pdm = t.phylogenetic_distance_matrix()
    taxa = list(t.taxon_namespace)
    if k >= len(taxa):
        return _from_dendropy(t)

    # seed with the most distant pair
    best = (None, None, -1.0)
    for i, a in enumerate(taxa):
        for b in taxa[i + 1:]:
            d = pdm.distance(a, b)
            if d > best[2]:
                best = (a, b, d)
    chosen = [best[0], best[1]]
    # farthest-point: repeatedly add the taxon maximising min-distance to the set
    while len(chosen) < k:
        cand, cand_d = None, -1.0
        for tx in taxa:
            if tx in chosen:
                continue
            dmin = min(pdm.distance(tx, c) for c in chosen)
            if dmin > cand_d:
                cand, cand_d = tx, dmin
        chosen.append(cand)

    keep = {tx.label for tx in chosen}
    t.retain_taxa_with_labels(sorted(keep))
    t.suppress_unifurcations()
    return _from_dendropy(t)


def star_tree(k: int, branch_length: float = 0.3, prefix: str = "sp") -> PhyloTree:
    """One internal node, ``k`` equidistant tips. The canonical empty-private-
    segment / degenerate-MRCA fixture from the ASR path-score history."""
    labels = [""] + [f"{prefix}{i}" for i in range(k)]
    edges = [(0, i + 1, branch_length) for i in range(k)]
    return PhyloTree(labels=labels, edges=edges, root=0)


def ladder_tree(k: int, branch_length: float = 0.15, prefix: str = "sp") -> PhyloTree:
    """Caterpillar / pectinate tree: every internal node has one tip child and
    one internal child. Maximises private-path-segment length asymmetry."""
    labels = [""]
    edges: list[tuple[int, int, float]] = []
    cur_internal = 0
    for i in range(k - 1):
        tip = len(labels)
        labels.append(f"{prefix}{i}")
        edges.append((cur_internal, tip, branch_length))
        if i < k - 2:
            nxt = len(labels)
            labels.append("")
            edges.append((cur_internal, nxt, branch_length))
            cur_internal = nxt
        else:
            last = len(labels)
            labels.append(f"{prefix}{k - 1}")
            edges.append((cur_internal, last, branch_length))
    return PhyloTree(labels=labels, edges=edges, root=0)


def scale_branches(tree: PhyloTree, factor: float) -> PhyloTree:
    return PhyloTree(
        labels=list(tree.labels),
        edges=[(p, c, bl * factor) for p, c, bl in tree.edges],
        root=tree.root,
    )


# --------------------------------------------------------------------------- #
# foreground sampling
# --------------------------------------------------------------------------- #
@dataclass
class Foreground:
    tips: list[str]                      # all foreground tips
    origin_nodes: list[int]              # nodes where the phenotype independently arises
    fg_edges: set[int]                   # indices into tree.edges that are "phenotype on"
    n_transitions: int = field(init=False)

    def __post_init__(self) -> None:
        self.n_transitions = len(self.origin_nodes)


def _node_depths(tree: PhyloTree) -> dict[int, float]:
    """Root-to-node distance for every node (preorder over tree.edges)."""
    d = {tree.root: 0.0}
    for p, c, bl in tree.edges:
        d[c] = d[p] + max(bl, 0.0)
    return d


def _clade_separation(tree: PhyloTree, a: int, b: int,
                      depths: dict[int, float], parent_of: dict[int, int]) -> float:
    """Patristic distance between two nodes = depth(a)+depth(b)-2*depth(lca)."""
    anc_a = set()
    x = a
    while True:
        anc_a.add(x)
        if x == tree.root:
            break
        x = parent_of[x]
    y = b
    while y not in anc_a:
        y = parent_of[y]
    return depths[a] + depths[b] - 2 * depths[y]


def _outgroup_tips(tree: PhyloTree, factor: float = 5.0) -> set[str]:
    """Tips whose terminal branch is a gross outlier (> factor x the 95th
    percentile of terminal branch lengths) — the deep outgroups. Foreground on
    these distorts contrast selection (an outgroup pairs trivially with anything)
    and is biologically odd. On the primate tree this is exactly Mus /
    Oryctolagus / Tupaia / Galeopterus."""
    bl = {c: b for _, c, b in tree.edges}
    leaves = tree._leaves()
    term = np.array([bl.get(l, 0.0) for l in leaves])
    if term.size == 0:
        return set()
    thr = factor * float(np.quantile(term, 0.95))
    return {tree.labels[l] for l in leaves if bl.get(l, 0.0) > thr}


def sample_foreground(tree: PhyloTree, n_transitions: int, rng: np.random.Generator,
                      *, min_clade: int = 1, max_clade: int = 3,
                      max_fg_fraction: float = 0.5,
                      exclude_tips: set[str] | None = None) -> Foreground:
    """Pick ``n_transitions`` disjoint, **phylogenetically well-separated** subtrees
    as independent origins of the foreground.

    The separation matters: production contrast selection (independent-contrasts +
    Dunn-gated pairing) merges origins that are too close, so poorly-spread planted
    clades yield < ``min_contrasts`` operative pairs and the run skips. Origins are
    chosen farthest-point: seed with a random candidate, then repeatedly add the
    disjoint candidate that maximises its minimum patristic distance to the ones
    already picked.
    """
    leaves = tree._leaves()
    n_tips = len(leaves)
    edge_index = {(p, c): i for i, (p, c, _) in enumerate(tree.edges)}
    parent_of = {c: p for p, c, _ in tree.edges}
    depths = _node_depths(tree)
    child_of = {c for _, c, _ in tree.edges}

    banned = set(exclude_tips or ()) | _outgroup_tips(tree)
    cand = []
    for node in range(tree.n_nodes):
        if node == tree.root and tree.root not in child_of:
            continue
        dt = tree.descendant_tips(node)
        if min_clade <= len(dt) <= max_clade and not (dt & banned):
            cand.append(node)
    rng.shuffle(cand)  # tie-break / seed randomness

    budget = int(n_tips * max_fg_fraction)
    chosen: list[int] = [cand[0]] if cand else []
    used_tips: set[str] = set(tree.descendant_tips(cand[0])) if cand else set()

    while len(chosen) < n_transitions:
        best, best_sep = None, -1.0
        for node in cand:
            if node in chosen:
                continue
            dt = tree.descendant_tips(node)
            if dt & used_tips or len(used_tips) + len(dt) > budget:
                continue
            sep = min(_clade_separation(tree, node, c, depths, parent_of) for c in chosen)
            if sep > best_sep:
                best, best_sep = node, sep
        if best is None:
            break
        chosen.append(best)
        used_tips |= tree.descendant_tips(best)

    if len(chosen) < n_transitions:
        raise ValueError(
            f"could only place {len(chosen)}/{n_transitions} disjoint foreground "
            f"clades on a {n_tips}-tip tree (clade size {min_clade}..{max_clade}, "
            f"<= {max_fg_fraction:.0%} of tips)"
        )

    # collect fg edges: for each origin, its subtending edge + all edges within
    fg_edges: set[int] = set()
    parent_of = {c: p for p, c, _ in tree.edges}
    for node in chosen:
        if node in parent_of:
            fg_edges.add(edge_index[(parent_of[node], node)])
        stack = [node]
        while stack:
            n = stack.pop()
            for c in tree.children_of(n):
                fg_edges.add(edge_index[(n, c)])
                stack.append(c)

    return Foreground(tips=sorted(used_tips), origin_nodes=chosen, fg_edges=fg_edges)
