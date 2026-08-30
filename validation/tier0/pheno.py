"""Tier 0 phenotype construction.

The real pipeline's operative hypothesis is a handful (~3) of phylogenetically
independent foreground/background contrast pairs that CONTRAST SELECTION picks
from the trait (see the real cancer run: `traitfile.tab` = 3 pairs). Tier 0
mirrors that directly: pick `n_pairs` well-separated foreground anchor tips, pair
each with its nearest background neighbour, and build a trait whose contrast
selection recovers exactly those pairs.

Two archetypes, differing only in how the trait column is emitted:

    binary : a 0/1 presence/absence code (`--trait_type ordinal`). Anchors = 1,
             everyone else = 0.  Discrete top/bottom candidate path.
    rate   : a `c / n` proportion with `n_pop` / `n_cases` count columns
             (CLASS 1). Anchors get a high rate with a tight Jeffreys CI, their
             partners a low rate tight CI, everyone else a mid rate wide CI so
             only the anchor/partner pairs are CI-separated candidates.

The planted molecular signal is placed on the anchor lineages' terminal edges
(each anchor is its own single-tip origin — "one pair per non-nested foreground").
For a null replicate the trait and pairs are still built (so the pipeline runs
identically) but `fg_edges` is empty, so `simulate()` plants nothing.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .trees import (
    Foreground,
    PhyloTree,
    _clade_separation,
    _node_depths,
    _outgroup_tips,
)


@dataclass
class Phenotype:
    archetype: str                          # "binary" | "rate"
    kind: str                               # "binary" | "continuous"
    values: dict[str, float]                # tip -> trait value (species with data only)
    foreground_tips: list[str]              # anchors — signal planted on their lineages
    origin_nodes: list[int]                 # anchor leaf node ids
    fg_edges: set[int]                      # anchor terminal edges (empty for a null)
    pairs: list[tuple[str, str]]            # (foreground anchor, background partner)
    n_pop: dict[str, int] | None = None     # rate archetype: sample size per species
    n_cases: dict[str, int] | None = None   # rate archetype: event count per species
    top_quantile: float = 0.90
    bottom_quantile: float = 0.10
    notes: dict = field(default_factory=dict)

    @property
    def n_transitions(self) -> int:
        return len(self.origin_nodes)

    def as_foreground(self) -> Foreground:
        return Foreground(tips=sorted(self.foreground_tips),
                          origin_nodes=list(self.origin_nodes),
                          fg_edges=set(self.fg_edges))


def _tip_index(tree: PhyloTree) -> dict[str, int]:
    return {tree.labels[l]: l for l in tree._leaves()}


def make_paired_foreground(
    tree: PhyloTree,
    n_pairs: int,
    rng: np.random.Generator,
    *,
    kind: str = "binary",
    planted: bool = True,
    fg_rate: float = 0.25,
    bg_rate: float = 0.02,
    mid_rate: float = 0.08,
    tight_n: int = 150,
    wide_n: int = 18,
    n_mid_species: int = 30,
) -> Phenotype:
    if kind not in ("binary", "rate"):
        raise ValueError(f"kind must be 'binary' or 'rate', got {kind!r}")

    tips = tree.tips
    idx_of = _tip_index(tree)
    depths = _node_depths(tree)
    parent_of = {c: p for p, c, _ in tree.edges}
    edge_index = {(p, c): i for i, (p, c, _) in enumerate(tree.edges)}
    banned = _outgroup_tips(tree)

    def pd(a: str, b: str) -> float:
        return _clade_separation(tree, idx_of[a], idx_of[b], depths, parent_of)

    # ── foreground anchors: farthest-point over non-outgroup tips ────────────
    pool = [t for t in tips if t not in banned]
    if len(pool) < 2 * n_pairs:
        raise ValueError(f"{len(pool)} usable tips < 2 * n_pairs ({2 * n_pairs})")
    rng.shuffle(pool)
    anchors = [pool[0]]
    while len(anchors) < n_pairs:
        best, best_sep = None, -1.0
        for t in pool:
            if t in anchors:
                continue
            sep = min(pd(t, a) for a in anchors)
            if sep > best_sep:
                best, best_sep = t, sep
        anchors.append(best)

    # ── background partner = nearest unused tip to each anchor ──────────────
    used = set(anchors)
    partners: list[str] = []
    for a in anchors:
        cand = [t for t in tips if t not in used and t not in banned]
        p = min(cand, key=lambda t: pd(a, t))
        partners.append(p)
        used.add(p)
    pairs = list(zip(anchors, partners))

    origin_nodes = [idx_of[a] for a in anchors]
    fg_edges: set[int] = set()
    if planted:
        for a in anchors:
            n = idx_of[a]
            if n in parent_of:
                fg_edges.add(edge_index[(parent_of[n], n)])

    anchor_set, partner_set = set(anchors), set(partners)

    # ── trait values ───────────────────────────────────────────────────────
    if kind == "binary":
        values = {t: (1.0 if t in anchor_set else 0.0) for t in tips}
        n_pop = n_cases = None
        kind_out = "binary"
    else:
        others = [t for t in tips if t not in used]
        rng.shuffle(others)
        data_sp = anchor_set | partner_set | set(others[: max(0, n_mid_species)])
        values, n_pop, n_cases = {}, {}, {}
        for t in tips:
            if t not in data_sp:
                continue
            if t in anchor_set:
                r, nn = fg_rate * (1.0 + rng.uniform(-0.15, 0.15)), tight_n
            elif t in partner_set:
                r, nn = bg_rate * (1.0 + rng.uniform(-0.30, 0.30)), tight_n
            else:
                r, nn = mid_rate * (1.0 + rng.uniform(-0.4, 0.4)), int(wide_n * rng.uniform(0.7, 1.5))
            nn = max(int(nn), 5)
            c = int(round(np.clip(r, 0.0, 1.0) * nn))
            values[t] = c / nn
            n_pop[t] = nn
            n_cases[t] = c
        kind_out = "continuous"

    return Phenotype(
        archetype=("binary" if kind == "binary" else "rate"),
        kind=kind_out,
        values=values,
        foreground_tips=sorted(anchor_set),
        origin_nodes=origin_nodes,
        fg_edges=fg_edges,
        pairs=pairs,
        n_pop=n_pop,
        n_cases=n_cases,
        top_quantile=0.90,
        bottom_quantile=0.10,
        notes={
            "n_pairs": n_pairs,
            "planted": planted,
            "anchors": list(anchors),
            "partners": list(partners),
            "n_species_with_data": len(values),
        },
    )
