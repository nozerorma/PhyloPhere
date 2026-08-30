"""Tier 0 phenotype construction — the two archetypes from DECISIONS.md D-T0-C.

    echo      : exact 0/1 presence-absence trait. Foreground = a set of non-nested
                clades sampled up front (trees.sample_foreground). RER runs binary.
    bodysize  : a Brownian-motion continuous trait simulated on the tree; the
                foreground is *emergent* = the top-quantile tips of that draw.
                RER runs continuous.

Both return a ``Phenotype``: the per-tip trait value(s), the operative foreground
tip set the molecular signal will be planted against, the independent "origin"
nodes, and the ``top_quantile`` / ``bottom_quantile`` to hand the pipeline so
that FADE's extreme-species pick and contrast selection's high/low cut both
land on the planted foreground.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .trees import Foreground, PhyloTree, sample_foreground


@dataclass
class Phenotype:
    archetype: str                       # "echo" | "bodysize"
    kind: str                            # "binary" | "continuous"
    values: dict[str, float]             # tip -> trait value
    foreground_tips: list[str]           # operative fg (signal planted on their lineages)
    origin_nodes: list[int]              # independent origins (for edge marking)
    fg_edges: set[int]                   # edges "phenotype on" (into simulate())
    top_quantile: float = 0.90
    bottom_quantile: float = 0.10
    bm_sigma: float = 1.0
    notes: dict = field(default_factory=dict)

    @property
    def n_transitions(self) -> int:
        return len(self.origin_nodes)

    def as_foreground(self) -> Foreground:
        return Foreground(tips=sorted(self.foreground_tips),
                          origin_nodes=list(self.origin_nodes),
                          fg_edges=set(self.fg_edges))


# --------------------------------------------------------------------------- #
# Brownian motion on a tree
# --------------------------------------------------------------------------- #
def simulate_bm(tree: PhyloTree, rng: np.random.Generator, sigma: float = 1.0,
                root: float = 0.0) -> dict[str, float]:
    """Exact BM: each child value = parent + N(0, sigma^2 * branch_length)."""
    val = {tree.root: root}
    for p, c, bl in tree.edges:
        val[c] = val[p] + rng.normal(0.0, sigma * np.sqrt(max(bl, 0.0)))
    leaves = tree._leaves()
    return {tree.labels[i]: float(val[i]) for i in leaves}


# --------------------------------------------------------------------------- #
# edge marking from a foreground tip set
# --------------------------------------------------------------------------- #
def _fg_edges_for_tips(tree: PhyloTree, fg_tips: set[str]) -> tuple[set[int], list[int]]:
    """A node is 'fg' if it is an fg tip or all its children are fg. An edge is fg
    if its child node is fg. Maximal fg subtrees are the independent origins."""
    leaves = set(tree._leaves())
    label = tree.labels
    children: dict[int, list[int]] = {}
    for p, c, _ in tree.edges:
        children.setdefault(p, []).append(c)

    is_fg: dict[int, bool] = {}

    def resolve(n: int) -> bool:
        if n in is_fg:
            return is_fg[n]
        if n in leaves:
            is_fg[n] = label[n] in fg_tips
        else:
            kids = children.get(n, [])
            is_fg[n] = bool(kids) and all(resolve(k) for k in kids)
        return is_fg[n]

    resolve(tree.root)

    edge_index = {(p, c): i for i, (p, c, _) in enumerate(tree.edges)}
    parent_of = {c: p for p, c, _ in tree.edges}
    fg_edges: set[int] = set()
    origins: list[int] = []
    for n, fg in is_fg.items():
        if not fg:
            continue
        # origin = fg node whose parent is not fg
        if n == tree.root or not is_fg.get(parent_of.get(n, tree.root), False):
            origins.append(n)
        if n in parent_of:
            fg_edges.add(edge_index[(parent_of[n], n)])
        for k in children.get(n, []):
            if is_fg.get(k):
                fg_edges.add(edge_index[(n, k)])
    return fg_edges, origins


# --------------------------------------------------------------------------- #
# archetypes
# --------------------------------------------------------------------------- #
def make_echo(tree: PhyloTree, n_transitions: int, rng: np.random.Generator,
              *, max_clade: int = 3) -> Phenotype:
    """Presence/absence archetype — **exact 0/1** ``--my_traits`` column.

    RER auto-detects this as binary (exactly two unique values) and runs
    ``foreground2Tree`` / the binary path. Contrast selection's production
    discrete path categorises it trivially (1s -> top, 0s -> bottom). The CAAS
    permulation works only because ``lean_contrast_selector.R`` was relaxed from
    ``trait > median`` to ``trait >= median`` (DECISIONS.md D-T0-E / the pipeline
    fix) — a minority-foreground 0/1 vector has ``median == 0`` and the strict
    form left ``low_sp`` empty.

    ``top_quantile`` / ``bottom_quantile`` are set from the foreground fraction so
    ``quantile(trait, top_quantile) == 1`` (picks exactly the 1s as "top") and
    ``quantile(trait, bottom_quantile) == 0`` (picks the 0s as "bottom"), which
    is also what FADE's ``EXTRACT_EXTREME_SPECIES`` uses.
    """
    fg = sample_foreground(tree, n_transitions, rng, max_clade=max_clade)
    fg_set = set(fg.tips)
    fg_frac = len(fg_set) / len(tree.tips)
    values = {t: (1.0 if t in fg_set else 0.0) for t in tree.tips}
    return Phenotype(
        archetype="echo", kind="binary", values=values,
        foreground_tips=sorted(fg_set), origin_nodes=list(fg.origin_nodes),
        fg_edges=set(fg.fg_edges),
        # top_quantile above (1 - fg_frac) => threshold falls on the 1s;
        # bottom_quantile below (1 - fg_frac) => threshold falls on the 0s.
        top_quantile=float(np.clip(1.0 - 0.5 * fg_frac, 0.55, 0.98)),
        bottom_quantile=float(np.clip((1.0 - fg_frac) * 0.5, 0.02, 0.45)),
        notes={"fg_fraction": fg_frac, "n_transitions_requested": n_transitions,
               "encoding": "exact 0/1"},
    )


def make_bodysize(tree: PhyloTree, rng: np.random.Generator, *, quantile: float = 0.15,
                  bm_sigma: float = 1.0) -> Phenotype:
    values = simulate_bm(tree, rng, sigma=bm_sigma)
    tips_sorted = sorted(values, key=lambda t: values[t])
    n = len(tips_sorted)
    k = max(2, round(n * quantile))
    fg_tips = set(tips_sorted[-k:])
    fg_edges, origins = _fg_edges_for_tips(tree, fg_tips)
    return Phenotype(
        archetype="bodysize", kind="continuous", values=values,
        foreground_tips=sorted(fg_tips), origin_nodes=origins, fg_edges=fg_edges,
        top_quantile=1.0 - quantile, bottom_quantile=quantile, bm_sigma=bm_sigma,
        notes={"quantile": quantile, "n_fg_tips": len(fg_tips),
               "n_emergent_origins": len(origins)},
    )


def make_null(tree: PhyloTree, rng: np.random.Generator, archetype: str,
              *, n_transitions: int = 3, quantile: float = 0.15) -> Phenotype:
    """A phenotype with NO planted molecular signal: the trait is drawn the same
    way as its archetype (so contrast selection / FADE / RER see a realistic
    input) but ``fg_edges`` is empty, so simulate() plants nothing."""
    if archetype == "echo":
        p = make_echo(tree, n_transitions, rng)
    else:
        p = make_bodysize(tree, rng, quantile=quantile)
    p.fg_edges = set()
    p.notes["null"] = True
    return p
