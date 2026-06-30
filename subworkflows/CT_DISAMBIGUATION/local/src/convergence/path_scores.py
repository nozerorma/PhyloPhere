#!/usr/bin/env python3
"""ASR path scoring: MRCA->root background validation of CAAS changes.

This module collapses three previously-separate signals — the binary ASR
conservation gate, the convergence score, and the parallel score — into a
single continuous per-position score derived directly from ancestral state
reconstruction (ASR).

Conceptual model
----------------
For each CAAS pair we already know the tip states (per phenotype side) and the
modal ancestral state at the pair's MRCA (``focal_state``). The question this
module answers is: *was each change isolated from the deeper background, or was
the "derived" state already present above the MRCA?*

To answer it we walk the chain of internal nodes from the pair's MRCA **up to
the root** and, at each node ``j`` (hop ``k`` above the MRCA), weight its ASR
posterior by ``w_k = 1 / (k + 1)``. Nodes immediately above the MRCA dominate;
nodes near the root contribute little. This is topological (hop) weighting,
deliberately **not** branch-length weighting: gene-tree branch lengths are in
substitutions/site and are not comparable across genes, and long branches
conflate evolutionary time with rate heterogeneity and reconstruction
uncertainty.

Two regimes per side
--------------------
For each phenotype side of a pair we compare the tip state to the MRCA's modal
state (both encoded in the active grouping scheme):

* **Changed side** (tip != MRCA): *signed isolation*. At each node we score
  ``P(ancestral) - P(derived)``. If the derived state was already present above
  the MRCA the signal goes negative and the score collapses toward 0. A
  ``contaminated`` flag is raised when the **modal** state at hop+1 (the node
  directly above the MRCA) already equals the derived state — a qualitative
  disqualifier indicating the "change" predates the MRCA.
* **Conserved side** (tip == MRCA): *not scored*. A side that retains the
  ancestral state is the contrast, not the signal. (``side_path_score`` still
  supports an unsigned-conservation mode for possible future use, but the
  aggregation below scores only changed sides.)

Participating and conserved pairs
----------------------------------
Convergence support comes only from *participating* pairs — those whose tip
diverged from their MRCA. The conserved trait side of a participating pair is the
contrast, not the signal, and is not scored (scoring it would let a clade that
merely held ancestral inflate the position, since conservation out-scores
isolation). A **conserved pair** (``conserved_pair`` from metadata) is a pair
where both tips retained the ancestral residue and the expected derived amino acid
was not acquired. It is scored by conservation-to-root and folded into the
``conservation_gate`` factor (below), which confirms or weakens the contrast but
never inflates the convergence signal.

Aggregation
-----------
``pair_score`` = mean isolation of a participating pair's changed side(s) (the
per-pair MRCA→root walk). Across participating pairs the position score is the
mean of pair scores multiplied by **five factors**:

* ``independence`` (shared-origin axis): the changed pairs only converged
  *independently* if their shared ancestors did not already carry the derived
  state. We locate the **LAC** nodes — the merge points of the changed pairs'
  MRCAs (internal nodes of the minimal subtree connecting them, at most n-1) —
  and take ``∏_L (1 - P(derived pool at L))``. If any merge point already shows a
  derived state the pairs below it inherited it rather than converging, and the
  factor (and score) crater. Because it ranges over the few LAC nodes, not the
  full root paths, the product is **depth-independent**. Read straight from
  ``P(derived)`` — no ancestral term — so a third residue at a node simply lowers
  ``P(derived)`` instead of being misread as ancestral (the §i edge case).

* ``mrca_diversity`` (parallel axis, continuous): how distinct are the
  ancestral starting points of the changed pairs? Compares the **posteriors** at
  each changed pair's MRCA (in encoded group space), not just modal states.
  Heavy posterior overlap => shared ancestral state (parallel — weaker, a shared
  mutational tendency could explain it); disjoint posteriors => independent
  origins (non-parallel — stronger evidence of selection). Count-independent (a
  mean over pairwise overlaps), so replication — already in the CAAS p-value —
  is not double-counted::

      diversity      = 1 - mean_pairwise_overlap(MRCA posteriors)   # 0..1
      diversity_mult = floor + (1 - floor) * diversity   # floor (def 0.75) .. 1.0

* ``derived_agreement`` (convergence axis): within each phenotype side, do the
  changed pairs agree on **one** derived (encoded) residue? Multiple distinct
  derived states on the same side is divergence, not convergence, and the
  probability that this is a genuine convergent signal is divided across them::

      derived_agreement = mean over changed sides of (1 / n_distinct_derived)

  One derived residue per side -> 1.0; a side splitting across k residues -> 1/k.
  Computed **per side**, so a pair changing on *both* sides (e.g. trait1->S,
  trait0->T) is not penalised — each side is internally coherent (n=1) and a
  both-sides change is a strong, interesting signal.

* ``conservation_gate`` (conserved-pair axis): the mean conservation-to-root of
  conserved pairs (both tips retained ancestral), softened to
  ``0.5 + 0.5 * cons``. A deeply conserved ancestral state -> ~1.0 (confirms the
  contrast); a pair that drifted -> down to 0.5 (weakens but never annihilates
  the real convergence). Neutral (1.0) when there are no conserved pairs (novel
  position — all lineages changed). The gate can only confirm/undermine, never
  boost — convergence support is the primary signal.

* ``core`` (replication axis): ``P(≥2 independent changes)`` — the combinatorial
  probability that at least two pairs changed independently in at least one
  phenotype direction (top or bottom). Computed per side via inclusion-exclusion
  on the per-pair path scores (each treated as P(change is genuine)), then
  combined across directions: ``core = 1 − (1 − core_top)(1 − core_bottom)``.
  A pair that changed on *both* sides contributes to both halves; opposite-only
  changes (one pair top-only, another bottom-only) correctly yield ``core = 0``
  because neither direction reaches ≥2.

    asr_path_score = core * independence * diversity_mult
                     * derived_agreement * conservation_gate

Most-interesting case D->A & S->A maximizes the multipliers (diversity 1.0,
agreement 1.0); parallel S->A & S->A scores diversity 0.75; a false "S->A & S->V"
is penalized on agreement (0.5); a conserved-CAAS with a clean control keeps
gate~1, a drifting control drops it toward 0.5. GS encoding handles biochemistry,
so agreement is tested on the scheme-encoded residues.

The module is intentionally free of PAML/IO dependencies so it can be unit
tested with synthetic trees and posteriors.

Author
------
Miguel Ramon Alonso — Evolutionary Genomics Lab, IBE-UPF
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

from src.biochem.grouping import get_grouping_scheme


# ── Encoding helpers ─────────────────────────────────────────────────────────
def encode_aa(aa: Optional[str], scheme: Optional[str]) -> Optional[str]:
    """Encode an amino acid in the active grouping scheme.

    Under the unscored scheme ("US") the raw residue is returned unchanged so
    single-residue comparisons are exact. Under a GS scheme the group label is
    returned (falling back to the raw residue if the residue is unknown).
    """
    if not aa:
        return None
    raw = str(aa).strip().upper()
    if not raw:
        return None
    if not scheme or scheme.upper() == "US":
        return raw
    return get_grouping_scheme(raw, scheme) or raw


def group_probability(
    distribution: Optional[Dict[str, float]],
    target_enc: Optional[str],
    scheme: Optional[str],
) -> float:
    """Sum posterior mass over all residues that encode to ``target_enc``.

    Aggregating the 20-AA PAML posterior into the scheme's group space makes the
    signal more robust under GS schemes: a node that is solidly "polar" under
    GS1 gives a cleaner signal than one whose probability is split across
    S/T/N/Q individually.
    """
    if not distribution or not target_enc:
        return 0.0
    total = 0.0
    for aa, prob in distribution.items():
        if encode_aa(aa, scheme) == target_enc:
            try:
                total += float(prob)
            except (TypeError, ValueError):
                continue
    return total


def modal_encoded(
    distribution: Optional[Dict[str, float]], scheme: Optional[str]
) -> Optional[str]:
    """Return the scheme-encoded modal residue of a node's posterior."""
    if not distribution:
        return None
    try:
        modal_aa = max(distribution.items(), key=lambda kv: kv[1])[0]
    except ValueError:
        return None
    return encode_aa(modal_aa, scheme)


def encoded_distribution(
    distribution: Optional[Dict[str, float]], scheme: Optional[str]
) -> Dict[str, float]:
    """Collapse a node's 20-AA posterior into scheme-encoded group space.

    Returns ``{encoded_group: summed_posterior}``. Under US this is essentially
    the raw distribution; under a GS scheme residues in the same biochemical
    group are merged, so two MRCAs that are both "polar" (but split across
    S/T/N/Q individually) read as the *same* ancestral state.
    """
    enc: Dict[str, float] = {}
    if not distribution:
        return enc
    for aa, prob in distribution.items():
        key = encode_aa(aa, scheme)
        if key is None:
            continue
        try:
            enc[key] = enc.get(key, 0.0) + float(prob)
        except (TypeError, ValueError):
            continue
    return enc


def posterior_overlap(dist_a: Dict[str, float], dist_b: Dict[str, float]) -> float:
    """Shared posterior mass (dot product) between two encoded distributions.

    ``Σ_g P_a(g) * P_b(g)`` over all groups. 1.0 = identical distributions (two
    MRCAs in the same ancestral state → parallel origin); 0.0 = disjoint support
    (genuinely different ancestral starting points → non-parallel, independent).
    """
    if not dist_a or not dist_b:
        return 0.0
    return sum(
        dist_a.get(g, 0.0) * dist_b.get(g, 0.0)
        for g in set(dist_a) & set(dist_b)
    )


# ── Tree helpers ─────────────────────────────────────────────────────────────
def build_node_index(root) -> Dict[int, Any]:
    """Map node_id -> TreeNode for the whole tree (single traversal)."""
    index: Dict[int, Any] = {}
    if root is None:
        return index
    stack = [root]
    while stack:
        node = stack.pop()
        if node is None:
            continue
        if getattr(node, "node_id", None) is not None:
            index[int(node.node_id)] = node
        stack.extend(getattr(node, "children", []) or [])
    return index


def node_dist(
    per_node_dist: Dict[Any, Dict[str, float]], node_id: Optional[int]
) -> Dict[str, float]:
    """Look up a node's posterior, tolerating int- or str-keyed maps.

    PAML posteriors are loaded with integer node ids, but the same map can
    arrive JSON-decoded with string keys. Centralising the fallback here keeps
    the three call sites (per-node walk, LAC product, MRCA-diversity) consistent
    and the ``per_node_dist`` key type honestly ``Any``.
    """
    if node_id is None:
        return {}
    return per_node_dist.get(node_id) or per_node_dist.get(str(node_id)) or {}


def path_to_root_ids(node_index: Dict[int, Any], mrca_id: Optional[int]) -> List[int]:
    """Node ids from the parent of the MRCA up to the root (MRCA excluded).

    The MRCA itself is the ancestral reference, so the background being tested is
    strictly the nodes *above* it. Element 0 is hop+1 (parent of the MRCA).
    """
    if mrca_id is None:
        return []
    node = node_index.get(int(mrca_id))
    if node is None:
        return []
    path: List[int] = []
    current = getattr(node, "parent", None)
    while current is not None:
        if getattr(current, "node_id", None) is not None:
            path.append(int(current.node_id))
        current = getattr(current, "parent", None)
    return path


def find_lca(
    node_index: Dict[int, Any], id_a: Optional[int], id_b: Optional[int]
) -> Optional[int]:
    """Lowest common ancestor of two nodes (the node where their lineages merge).

    Used to locate the **LAC** (last ancestral common) nodes of a set of pair
    MRCAs — the shared merge points where, if the derived state is already
    present, the convergence is not independent across those pairs.
    """
    if id_a is None or id_b is None:
        return None
    chain_a = [int(id_a)] + path_to_root_ids(node_index, id_a)  # node itself .. root
    set_b = set([int(id_b)] + path_to_root_ids(node_index, id_b))
    for n in chain_a:
        if n in set_b:
            return n
    return None


# ── Core per-side score ──────────────────────────────────────────────────────
# Score returned when the MRCA sits at/near the root (no background above it):
# there is no evidence either way, so we return a neutral weak score.
EMPTY_PATH_SCORE = 0.5


def side_path_score(
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    mrca_id: Optional[int],
    ancestral_enc: str,
    derived_enc: Optional[str],
    scheme: Optional[str],
    is_changed: bool,
    stop_at_id: Optional[int] = None,
) -> Tuple[float, bool]:
    """Score private isolation (changed) or global conservation (conserved).

    For a changed side (is_changed=True), walks the private segment from the
    parent of the MRCA up to stop_at_id (exclusive) and computes the product
    of (1 - P(derived)). If the private segment is empty (count == 0, e.g. sibling
    merge), returns 1.0 (no private contamination).

    For a conserved side (is_changed=False), walks the entire path to the root
    and computes the unweighted mean of P(ancestral) (no stop_at_id).
    """
    path = path_to_root_ids(node_index, mrca_id)
    if not path:
        return EMPTY_PATH_SCORE, False

    score = 1.0
    total_anc = 0.0
    count = 0
    contaminated = False

    for k, node_id in enumerate(path, start=1):  # k=1 -> parent of MRCA (hop+1)
        if is_changed and stop_at_id is not None and node_id == stop_at_id:
            break  # reached the shared LAC; stop private segment walk

        dist = node_dist(per_node_dist, node_id)
        p_anc = group_probability(dist, ancestral_enc, scheme)

        if is_changed:
            p_der = group_probability(dist, derived_enc, scheme)
            score *= max(0.0, 1.0 - p_der)
            if k == 1 and modal_encoded(dist, scheme) == derived_enc and derived_enc:
                contaminated = True
        else:
            total_anc += p_anc

        count += 1

    if is_changed:
        if count == 0:
            return 1.0, contaminated # no private nodes to contradict
        return max(0.0, min(1.0, score)), contaminated
    else:
        if count == 0:
            return EMPTY_PATH_SCORE, False
        return max(0.0, min(1.0, total_anc / count)), False


# ── Conserved-pair id parsing ────────────────────────────────────────────────
def parse_conserved_ids(conserved_pair: Optional[str], n_pairs: int) -> List[int]:
    """Parse the comma-separated conserved-pair id string into validated ints."""
    if not conserved_pair:
        return []
    ids: List[int] = []
    for token in str(conserved_pair).split(","):
        token = token.strip()
        if not token:
            continue
        try:
            ids.append(int(float(token)))
        except (ValueError, TypeError):
            continue
    return [p for p in ids if 1 <= p <= n_pairs]


# ── Position-level aggregation ───────────────────────────────────────────────
def compute_asr_path_score(
    pair_details: Optional[List[Dict[str, Any]]],
    per_node_dist: Dict[Any, Dict[str, float]],
    node_index: Dict[int, Any],
    scheme: Optional[str],
    is_conserved_meta: bool,
    conserved_pair: Optional[str],
    diversity_floor: float = 0.75,
) -> Dict[str, Any]:
    """Compute the unified ASR path score for one CAAS position row.

    Args:
        pair_details: list of per-pair dicts (``pair_id``, ``node_id``,
            ``focal_state``, ``top_tip_mode``, ``bottom_tip_mode``).
        per_node_dist: ``node_id -> {aa: posterior}`` for the focal site.
        node_index: ``node_id -> TreeNode`` (from :func:`build_node_index`).
        scheme: active grouping scheme (``"US"`` or a GS label).
        is_conserved_meta: whether this position has conserved pairs listed in
            metadata (both tips show potential ancestral state conservation).
        conserved_pair: comma-separated id(s) of the conserved pair(s). These
            are scored by conservation-to-root and folded into the
            conservation_gate factor, not counted as convergence evidence.
        diversity_floor: lower bound of ``diversity_mult``. A pure-parallel
            position (identical MRCA posteriors) gets exactly this; a maximally
            independent one gets 1.0. Default 0.75.

    Returns:
        Dict with ``asr_path_score`` (position level, 0-1), ``independence``,
        ``mrca_diversity``, ``derived_agreement``, ``conservation_gate``,
        ``core``, ``pair_scores`` ({pair_id: score}) and ``pair_contaminated``.
    """
    pairs = pair_details or []
    n_pairs = len(pairs)
    conserved_ids = (
        set(parse_conserved_ids(conserved_pair, n_pairs)) if is_conserved_meta else set()
    )

    pair_contaminated: Dict[int, bool] = {}
    conserved_conservations: List[float] = []  # conservation-to-root of conserved pairs
    # Derived (encoded) residues of changed tips, kept per phenotype side so
    # within-side divergence (the non-convergent case) is assessed separately
    # from a pair changing on *both* sides (a strong convergent signal).
    derived_by_side: Dict[str, List[str]] = {"top_tip_mode": [], "bottom_tip_mode": []}
    derived_pool: set = set()  # all encoded derived states across changed pairs

    # ── Phase 1: classify pairs ──────────────────────────────────────────────
    # Conserved pairs (conserved_ids) did not acquire the expected derived amino
    # acid — scored by conservation-to-root and folded into conservation_gate.
    # Changed pairs are collected here (sides not scored yet: the per-pair core
    # walk needs the LAC merge points computed below to know where to stop).
    changed: List[Dict[str, Any]] = []  # {pid, mrca_id, anc_enc, sides:[(key,enc)]}
    for pair in pairs:
        pid = pair.get("pair_id")
        mrca_id = pair.get("node_id")
        anc_enc = encode_aa(pair.get("focal_state"), scheme)
        if pid is None or mrca_id is None or anc_enc is None:
            continue

        if pid in conserved_ids:
            cons, _ = side_path_score(
                node_index, per_node_dist, mrca_id, anc_enc, None, scheme,
                is_changed=False,
            )
            conserved_conservations.append(cons)
            continue

        sides: List[Tuple[str, str]] = []
        for side_key in ("top_tip_mode", "bottom_tip_mode"):
            tip_enc = encode_aa(pair.get(side_key), scheme)
            if tip_enc is None or tip_enc == anc_enc:
                continue  # missing or conserved side → not scored
            sides.append((side_key, tip_enc))
            derived_by_side[side_key].append(tip_enc)
            derived_pool.add(tip_enc)

        if sides:
            changed.append(
                {"pid": pid, "mrca_id": int(mrca_id), "anc_enc": anc_enc, "sides": sides}
            )

    # Conservation gate: conserved pairs that held ancestral deeply (high
    # conservation) confirm the contrast (gate→1); pairs that drifted weaken
    # but never annihilate the real convergence (softened floor 0.5).
    # Novel positions (all lineages changed, no conserved pairs) → gate=1.0.
    if conserved_conservations:
        conservation_gate = 0.5 + 0.5 * (
            sum(conserved_conservations) / len(conserved_conservations)
        )
    else:
        conservation_gate = 1.0  # novel position: all lineages changed

    if not changed:
        return {
            "asr_path_score": 0.0,
            "independence": 1.0,
            "mrca_diversity": 0.0,
            "derived_agreement": 0.0,
            "conservation_gate": conservation_gate,
            "core": 0.0,
            "pair_scores": {},
            "pair_contaminated": {},
        }

    # ── LAC merge points ─────────────────────────────────────────────────────
    # The LAC nodes are the merge points of the changed pairs' MRCAs (the
    # internal nodes of the minimal subtree connecting them). At most n-1 of
    # them, so a *product* over LAC nodes is depth-independent — unlike a product
    # over the full MRCA→root paths, which would collapse on deep trees.
    changed_mrcas = [c["mrca_id"] for c in changed]
    lac_nodes: set = set()
    for a, b in combinations(changed_mrcas, 2):
        lac = find_lca(node_index, a, b)
        if lac is not None:
            lac_nodes.add(lac)

    # ── Phase 2: per-pair core isolation (each MRCA → root) ───────────────────
    pair_scores: Dict[int, float] = {}
    # Per-side scores: needed for directional core (see below). A pair that
    # changed only on top contributes to top_pair_scores; one that changed on
    # both sides contributes to both. This prevents a pair changing top→S and a
    # different pair changing bottom→T from being counted as "2 convergent
    # changes" — they are changes in opposite directions.
    top_pair_scores: Dict[int, float] = {}
    bottom_pair_scores: Dict[int, float] = {}
    for c in changed:
        # Nearest LAC = deepest merge point on this pair's root-path. The private
        # segment (below it) is the pair's own, judged by the core walk; the
        # shared segment (at/above it) is judged by the independence product.
        stop_at = next(
            (n for n in path_to_root_ids(node_index, c["mrca_id"]) if n in lac_nodes),
            None,
        )
        side_scores: List[float] = []
        pair_contam = False
        for side_key, tip_enc in c["sides"]:
            score, contam = side_path_score(
                node_index, per_node_dist, c["mrca_id"], c["anc_enc"], tip_enc,
                scheme, is_changed=True, stop_at_id=stop_at,
            )
            side_scores.append(score)
            pair_contam = pair_contam or contam
            if side_key == "top_tip_mode":
                top_pair_scores[c["pid"]] = score
            else:
                bottom_pair_scores[c["pid"]] = score
        pair_scores[c["pid"]] = sum(side_scores) / len(side_scores)
        pair_contaminated[c["pid"]] = pair_contam

    # ── Independence: derived state must be ABSENT at the shared merge points ──
    # For each LAC node, (1 - P(derived_pool there)). Product across LAC nodes:
    # if any merge point already carries a derived state, the pairs below it did
    # not converge independently — they inherited it — and the score craters.
    # This is the direct "probability the pattern predates the shared ancestor"
    # check, read straight from the posteriors (no ancestral term needed, so a
    # third residue at a node simply lowers P(derived) rather than misreading).
    independence = 1.0
    for lac in lac_nodes:
        dist = node_dist(per_node_dist, lac)
        p_derived = sum(group_probability(dist, d, scheme) for d in derived_pool)
        independence *= max(0.0, 1.0 - p_derived)

    # Parallel axis (continuous): how distinct are the ancestral starting points
    # of the changed pairs? Compares the *posteriors* at each MRCA in encoded
    # group space. Heavy overlap → same ancestral state (parallel — weaker, a
    # shared mutational tendency could explain it); disjoint → independent origins
    # (non-parallel — stronger). Count-independent (mean over pairwise overlaps):
    # replication already lives in the CAAS p-value.
    if len(changed_mrcas) >= 2:
        mrca_dists = [
            encoded_distribution(node_dist(per_node_dist, m), scheme)
            for m in changed_mrcas
        ]
        overlaps = [posterior_overlap(a, b) for a, b in combinations(mrca_dists, 2)]
        mean_overlap = sum(overlaps) / len(overlaps) if overlaps else 1.0
        diversity = 1.0 - mean_overlap  # 0 = all parallel, 1 = all independent
    else:
        diversity = 0.5  # single changed pair: parallelism not assessable
    diversity_mult = diversity_floor + (1.0 - diversity_floor) * diversity

    # Convergence axis: within each phenotype side, do the changed pairs agree on
    # *one* derived residue? Multiple distinct derived states on the same side is
    # divergence, not convergence, and the "probability this is a real convergent
    # signal" is divided across them: derived_agreement = 1 / n_distinct_derived.
    # Computed per side and averaged, so a pair changing on BOTH sides (e.g.
    # trait1→S, trait0→T) is NOT penalised — each side is internally coherent
    # (n=1) and a both-sides change is a strong, interesting signal.
    side_agreements: List[float] = []
    for side_derived in derived_by_side.values():
        if side_derived:  # >=1 changed pair on this phenotype side
            n_distinct_derived = len(set(side_derived))
            side_agreements.append(1.0 / n_distinct_derived)
    if side_agreements:
        derived_agreement = sum(side_agreements) / len(side_agreements)
    else:
        derived_agreement = 1.0  # no changed sides (handled above by empty `changed`)

    # ── Combinatorial Core Score: P(>= 2 independent changes) ──────────────────
    # Computed *per phenotype side* and then combined via union. This is the
    # critical directional check: a pair changing only on the top side and a
    # different pair changing only on the bottom side are changes in opposite
    # directions — not convergence. Without the per-side split they would both
    # enter p_list and produce P(≥2) ≈ 1 even though neither direction has ≥2
    # independent events.

    def _p_at_least_2(p_list: List[float]) -> float:
        if len(p_list) < 2:
            return 0.0
        p0 = 1.0
        for p in p_list:
            p0 *= (1.0 - p)
        p1 = 0.0
        for i in range(len(p_list)):
            term = p_list[i]
            for j in range(len(p_list)):
                if j != i:
                    term *= (1.0 - p_list[j])
            p1 += term
        return max(0.0, min(1.0, 1.0 - p0 - p1))

    core_top = _p_at_least_2(list(top_pair_scores.values()))
    core_bottom = _p_at_least_2(list(bottom_pair_scores.values()))
    core = 1.0 - (1.0 - core_top) * (1.0 - core_bottom)

    asr_path_score = max(
        0.0,
        min(
            1.0,
            core
            * independence
            * diversity_mult
            * derived_agreement
            * conservation_gate,
        ),
    )

    return {
        "asr_path_score": asr_path_score,
        "independence": independence,
        "mrca_diversity": diversity,
        "derived_agreement": derived_agreement,
        "conservation_gate": conservation_gate,
        "core": core,
        "pair_scores": pair_scores,
        "pair_contaminated": pair_contaminated,
    }
