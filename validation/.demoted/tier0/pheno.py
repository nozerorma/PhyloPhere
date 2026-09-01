"""Tier 0 phenotype construction — evolutionary-model planting.

The circularity we removed: previous versions hand-picked artificially
independent foreground tips and engineered a trait around them, so contrast
selection "recovered" them by construction.

Instead: simulate a latent trait by **Brownian motion on the tree rescaled by
Pagel's lambda** — the same family the permulation draws its nulls from
(`permulations.R` `simpermvec`). lambda tunes phylogenetic clustering:

    lambda = 0   internal branches collapse -> star -> latent is iid ->
                 the trait extremes are scattered across the tree
    lambda = 1   full BM -> the trait extremes clump into clades
    lambda = 0.5 in between (the real cancer-prevalence trait fitted ~0.74)

The **true foreground** is whatever comes out at the tail of that latent draw,
reduced to phylogenetically independent origins (Dunn >= 1). The observed trait
handed to the pipeline is that latent + sampling noise. The molecular signal is
planted on the TRUE foreground lineages by a separate process (the Gillespie
substitution sim), so contrast-recovery (operative fg vs true pairs) is an
honest measurement that varies with lambda, not a tautology.

Each planted gene is planted in one direction, ``top`` or ``bottom`` — the
pipeline scores a convergent change in the high-trait species and in the
low-trait species symmetrically (`change_top` / `change_bottom`).

Three archetypes exercise the three contrast-selection candidate paths; each gets
data for EVERY (non-outgroup) tip of the (pruned) tree:
    binary     : threshold the latent -> 0/1 (`--trait_type ordinal`), a few bits
                 flipped -> discrete top/bottom candidate path
    rate       : c ~ Binomial(n, rate(latent percentile)), n ~ U(25, 70), with
                 n_pop / n_cases columns -> Jeffreys-CI non-overlap path (CLASS 1)
    continuous : the raw latent + observation noise, handed straight -> quantile
                 (discrete_method) discretisation path
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


class NotEnoughPairs(RuntimeError):
    """The latent draw produced < min_pairs independent tail origins."""


@dataclass
class Phenotype:
    archetype: str                          # "binary" | "rate" | "continuous"
    kind: str                               # "binary" | "continuous"
    lam: float                              # Pagel's lambda used for the latent
    values: dict[str, float]                # tip -> observed trait value
    pairs: list[tuple[str, str]]            # TRUE (top anchor, bottom partner) pairs
    anchor_nodes: list[int]
    partner_nodes: list[int]
    anchor_edges: set[int]
    partner_edges: set[int]
    n_pop: dict[str, int] | None = None
    n_cases: dict[str, int] | None = None
    notes: dict = field(default_factory=dict)

    @property
    def foreground_tips(self) -> list[str]:
        return sorted(a for a, _ in self.pairs)

    @property
    def n_transitions(self) -> int:
        return len(self.pairs)

    def foreground(self, direction: str) -> Foreground:
        if direction == "top":
            return Foreground(tips=[a for a, _ in self.pairs],
                              origin_nodes=list(self.anchor_nodes),
                              fg_edges=set(self.anchor_edges))
        if direction == "bottom":
            return Foreground(tips=[b for _, b in self.pairs],
                              origin_nodes=list(self.partner_nodes),
                              fg_edges=set(self.partner_edges))
        raise ValueError(f"direction must be 'top' or 'bottom', got {direction!r}")


def _tip_index(tree: PhyloTree) -> dict[str, int]:
    return {tree.labels[l]: l for l in tree._leaves()}


def _lambda_rescale(tree: PhyloTree, lam: float) -> PhyloTree:
    """Pagel's lambda: internal branches x lambda; each terminal branch extended
    so every root-to-tip distance is preserved. lambda in [0, 1]."""
    lam = float(np.clip(lam, 0.0, 1.0))
    leaves = set(tree._leaves())
    depth_parent = _node_depths(tree)  # root-to-node on the ORIGINAL tree
    new_edges = []
    for p, c, bl in tree.edges:
        if c in leaves:
            new_edges.append((p, c, bl + (1.0 - lam) * depth_parent[p]))
        else:
            new_edges.append((p, c, bl * lam))
    return PhyloTree(labels=list(tree.labels), edges=new_edges, root=tree.root)


def _bm_latent(tree: PhyloTree, rng: np.random.Generator) -> dict[str, float]:
    """Exact BM down the tree: child = parent + N(0, branch_length)."""
    v = {tree.root: 0.0}
    for p, c, bl in tree.edges:
        v[c] = v[p] + rng.normal(0.0, np.sqrt(max(bl, 0.0)))
    return {tree.labels[i]: float(v[i]) for i in tree._leaves()}


def _independent_pairs(
    tree: PhyloTree,
    top_ranked: list[str],
    bot_ranked: list[str],
    *,
    hard_cap: int = 20,
) -> list[tuple[str, str]]:
    """Greedily build ALL (top, bottom) pairs that are mutually Dunn>=1
    independent, up to hard_cap. The caller slices the fixed number it uses and
    keeps len(this) as the diagnostic "how many contrasts the trait supports".

    top_ranked / bot_ranked: tail tips, most-extreme first. Each top tip is
    paired with its nearest available bottom tip; a pair joins only if every
    already-placed cluster stays Dunn>=1 from it.
    """
    idx = _tip_index(tree)
    depths = _node_depths(tree)
    parent_of = {c: p for p, c, _ in tree.edges}

    def pd(a: str, b: str) -> float:
        return _clade_separation(tree, idx[a], idx[b], depths, parent_of)

    chosen: list[tuple[str, str]] = []
    used: set[str] = set()
    for a in top_ranked:
        if a in used or len(chosen) >= hard_cap:
            continue
        cand_b = [b for b in bot_ranked if b not in used and b != a]
        if not cand_b:
            break
        b = min(cand_b, key=lambda x: pd(a, x))
        new = {a, b}
        ok = True
        for (ta, tb) in chosen:
            intra = max(pd(ta, tb), pd(a, b))
            inter = min(pd(x, y) for x in new for y in {ta, tb})
            if intra > 0 and inter / intra < 1.0:
                ok = False
                break
        if ok:
            chosen.append((a, b))
            used |= new
    return chosen


def make_lambda_foreground(
    tree: PhyloTree,
    n_pairs: int,
    rng: np.random.Generator,
    *,
    kind: str = "binary",
    planted: bool = True,
    lam: float = 0.5,
    tail_frac: float = 0.25,
    max_tries: int = 60,
    fg_rate: float = 0.22,
    bg_rate: float = 0.03,
    n_min: int = 25,
    n_max: int = 70,
    n_mid_species: int = 30,
    binary_flip: int = 1,
) -> Phenotype:
    """`n_pairs` is FIXED (used = min(n_pairs, possible)); the number of
    independent pairs the latent could support is recorded as
    notes['n_possible_pairs'] — the interesting lambda-dependent diagnostic.
    Resamples the latent until it supports at least `n_pairs`."""
    if kind not in ("binary", "rate", "continuous"):
        raise ValueError(f"kind must be 'binary' | 'rate' | 'continuous', got {kind!r}")

    tips = tree.tips
    idx_of = _tip_index(tree)
    parent_of = {c: p for p, c, _ in tree.edges}
    edge_index = {(p, c): i for i, (p, c, _) in enumerate(tree.edges)}
    banned = _outgroup_tips(tree)
    rescaled = _lambda_rescale(tree, lam)

    # ── latent BM draw -> independent tail pairs (resample until it supports
    #    n_pairs, so every replicate contributes the SAME contrast count) ─────
    all_pairs: list[tuple[str, str]] = []
    latent: dict[str, float] = {}
    best: list[tuple[str, str]] = []
    for _try in range(max_tries):
        latent = _bm_latent(rescaled, rng)
        ranked = sorted((t for t in tips if t not in banned), key=lambda t: latent[t])
        k = max(n_pairs + 2, int(len(ranked) * tail_frac))
        bot_ranked = ranked[:k]
        top_ranked = list(reversed(ranked[-k:]))
        all_pairs = _independent_pairs(tree, top_ranked, bot_ranked)
        if len(all_pairs) > len(best):
            best, best_latent = all_pairs, dict(latent)
        if len(all_pairs) >= n_pairs:
            break
    if len(all_pairs) < n_pairs:
        all_pairs, latent = best, best_latent
    if len(all_pairs) < n_pairs:
        raise NotEnoughPairs(
            f"lambda={lam}: only {len(all_pairs)} independent tail pairs after "
            f"{max_tries} draws (need {n_pairs})"
        )
    n_possible = len(all_pairs)
    pairs = all_pairs[:n_pairs]                 # FIXED count -> comparable

    anchors = [a for a, _ in pairs]
    partners = [b for _, b in pairs]
    anchor_set, partner_set = set(anchors), set(partners)
    used = anchor_set | partner_set

    def _terms(names) -> set[int]:
        return {edge_index[(parent_of[idx_of[t]], idx_of[t])]
                for t in names if idx_of[t] in parent_of}

    anchor_edges = _terms(anchors) if planted else set()
    partner_edges = _terms(partners) if planted else set()

    # ── observed trait: latent + sampling noise ──────────────────────────
    lat_vals = np.array([latent[t] for t in tips])
    lat_rank = {t: r / (len(tips) - 1) for r, t in
                enumerate(sorted(tips, key=lambda t: latent[t]))}

    data_tips = [t for t in tips if t not in banned]   # every non-outgroup tip carries data
    n_pop = n_cases = None
    if kind == "binary":
        thr = np.quantile([latent[t] for t in data_tips], 1.0 - tail_frac)
        val = {t: (1.0 if latent[t] >= thr else 0.0) for t in data_tips}
        for a in list(rng.permutation(anchors))[:binary_flip]:      # missed foreground
            val[a] = 0.0
        pos_pool = [t for t in data_tips if t not in used and val[t] == 0.0]
        rng.shuffle(pos_pool)
        for t in pos_pool[:binary_flip]:                            # spurious foreground
            val[t] = 1.0
        values = val
        kind_out = "binary"
    elif kind == "rate":
        values, n_pop, n_cases = {}, {}, {}
        for t in data_tips:
            rate = bg_rate + (fg_rate - bg_rate) * (lat_rank[t] ** 2)
            nn = int(rng.integers(n_min, n_max + 1))
            c = int(rng.binomial(nn, min(rate, 1.0)))
            values[t] = c / nn
            n_pop[t] = nn
            n_cases[t] = c
        kind_out = "continuous"
    else:  # continuous — raw value + observation noise, quantile-discretised by the pipeline
        sd = float(np.std([latent[t] for t in data_tips])) or 1.0
        values = {t: round(latent[t] + rng.normal(0.0, 0.35 * sd), 6) for t in data_tips}
        kind_out = "continuous"

    return Phenotype(
        archetype=kind,
        kind=kind_out,
        lam=lam,
        values=values,
        pairs=pairs,
        anchor_nodes=[idx_of[a] for a in anchors],
        partner_nodes=[idx_of[b] for b in partners],
        anchor_edges=anchor_edges,
        partner_edges=partner_edges,
        n_pop=n_pop,
        n_cases=n_cases,
        notes={
            "lambda": lam, "planted": planted,
            "n_pairs": len(pairs),                 # FIXED, used in the analysis
            "n_possible_pairs": n_possible,        # how many the latent supported
            "tries": _try + 1,
            "anchors": anchors, "partners": partners,
            "n_species_with_data": len(values),
        },
    )
