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

To answer it we walk internal nodes bounded by two regions, never the full
MRCA-to-root path: the pair's own **private segment** (from the node directly
above its MRCA up to, but excluding, the nearest point where its lineage
merges with another changed pair's — the LAC), and the **LAC nodes**
themselves. Nodes above the LAC are not examined. At each node visited, we
read its ASR posterior directly (no hop weighting — PAML's posterior
uncertainty already grows toward the root, so depth robustness instead comes
from bounding the walk to these two regions rather than weighting a longer one).

Two regimes per side
--------------------
For each phenotype side of a pair we compare the tip state to the MRCA's modal
state (both encoded in the active grouping scheme):

* **Changed side** (tip != MRCA): *signed isolation*. At each node we score
  ``1 - P(derived)``. If the derived state was already present above
  the MRCA the product collapses toward 0. A ``contaminated`` flag is raised
  when the **modal** state at hop+1 (the node directly above the MRCA) already
  equals the derived state — a qualitative disqualifier indicating the
  "change" predates the MRCA.
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
Five factors, computed per CAAS position and multiplied together:

* ``independence`` (shared-origin axis): the changed pairs only converged
  *independently* if their shared ancestors did not already carry the derived
  state. We locate the **LAC** nodes — the merge points of the changed pairs'
  MRCAs (internal nodes of the minimal subtree connecting them, at most n-1) —
  and take ``∏_L (1 - P(derived pool at L))``. If any merge point already shows a
  derived state the pairs below it inherited it rather than converging, and the
  factor (and score) crater. Because it ranges over the few LAC nodes, not the
  full root paths, the product is **depth-independent**. Read straight from
  ``P(derived)`` — no ancestral term — so a third residue at a node simply lowers
  ``P(derived)`` instead of being misread as ancestral.

* ``core`` (replication axis): ``P(≥2 independent changes)`` — the combinatorial
  probability that at least two pairs changed independently in at least one
  phenotype direction (top or bottom). Computed per side via inclusion-exclusion
  on the per-pair private-isolation scores (each treated as P(change is
  genuine)), then combined across directions:
  ``core = 1 − (1 − core_top)(1 − core_bottom)``. A pair that changed on *both*
  sides contributes to both halves; opposite-only changes (one pair top-only,
  another bottom-only) correctly yield ``core = 0`` because neither direction
  reaches ≥2.

  ``independence`` and ``core`` together form the **replication** tier —
  ``P(shared ancestor clean AND ≥2 private-clean events)``. Multiplying them is
  the correct joint probability, not two independent checks: private-segment
  nodes and LAC nodes are disjoint, conditionally-independent parts of the
  tree, so ``P(A) · P(B|A) = P(A ∩ B)``.

* ``mrca_diversity`` (parallel axis, continuous): did the changed pairs'
  private segments (the same MRCA→LAC walks ``core`` scores) pass through
  genuinely different ancestral backgrounds, or does one pair's own MRCA
  state show up somewhere along another pair's walk? For every pair of
  changed pairs, both directions are checked — does A's MRCA state appear
  anywhere in B's private segment, and does B's appear anywhere in A's —
  using the same worst-case bound as ``core``/``independence`` at each node
  visited, combined via the same "found somewhere" union logic. The two
  directions are unioned into one "these share a background" probability per
  pairwise comparison; the position's diversity is the mean over every
  pairwise comparison among changed pairs. Heavy shared-background evidence
  => parallel origin (weaker — a shared mutational tendency could explain
  it); genuinely absent from each other's walks => independent origins
  (non-parallel — stronger evidence of selection). Replication is already
  captured by ``core``, so it is not double-counted here::

      diversity      = 1 - mean_pairwise("shared background" probability)   # 0..1
      diversity_mult = floor + (1 - floor) * diversity   # floor (def 0.75) .. 1.0

* ``derived_agreement`` (convergence axis): within each phenotype side that has
  at least two changed pairs, what fraction land on the position's single most
  common (plurality) derived residue? A side with exactly one changed pair has
  nothing to agree or disagree with and is excluded from the average entirely
  — "agreement" is a ≥2-observer concept. Computed from plain pair counts, not
  weighted by any ancestral-reconstruction quantity: this axis is a fact about
  observed tip residues (directly sequenced, no posterior), and weighting it by
  a posterior-derived score would make it redundant with ``core``/
  ``independence``, which already own the ancestral-uncertainty question::

      concentration_s = max_r( count of pairs landing on residue r )
                         / count of changed pairs on side s          # in (0, 1]
      derived_agreement = mean over qualifying sides of concentration_s

  ``diversity_mult`` and ``derived_agreement`` together form the **strength**
  tier — two independent discounts on the *quality* of a confirmed replication
  event: did the lineages start from different places, and did they land on
  the same place.

* ``conservation_gate`` (conserved-pair axis): the mean conservation-to-root of
  conserved pairs (both tips retained ancestral), softened to
  ``0.5 + 0.5 * cons``. A deeply conserved ancestral state -> ~1.0 (confirms the
  contrast); a pair that drifted -> down to 0.5 (weakens but never annihilates
  the real convergence). Neutral (1.0) when there are no conserved pairs (novel
  position — all lineages changed). The gate can only confirm/undermine, never
  boost, and is computed from a disjoint set of pairs — it stays outside the
  replication/strength grouping entirely::

      replication = independence * core
      strength    = diversity_mult * derived_agreement
      asr_path_score = replication * strength * conservation_gate

Most-interesting case D->A & S->A maximizes the multipliers (diversity 1.0,
agreement 1.0); parallel S->A & S->A scores diversity 0.75; a false "S->A & S->V"
is penalized on agreement; a conserved-CAAS with a clean control keeps gate~1, a
drifting control drops it toward 0.5. GS encoding handles biochemistry, so
agreement is tested on the scheme-encoded residues.

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


def worst_case_group_probability(
    distribution: Optional[Dict[str, float]],
    target_enc: Optional[str],
    scheme: Optional[str],
) -> float:
    """Upper bound on posterior mass for ``target_enc`` given a possibly
    incomplete distribution.

    Returns the exact recorded mass if ``target_enc`` is present. Otherwise
    returns the total *unrecorded* mass (``1 - sum(recorded)``) as a ceiling:
    everything not accounted for could, in the worst case, belong to
    ``target_enc``. This requires no assumption about how the unrecorded mass
    is actually distributed among the residues that aren't recorded, unlike a
    uniform-fill estimate -- it is a guaranteed bound, not a guess.
    """
    if not target_enc:
        return 0.0
    enc = encoded_distribution(distribution, scheme)
    if target_enc in enc:
        return enc[target_enc]
    return max(0.0, 1.0 - sum(enc.values()))


def worst_case_any_group_probability(
    distribution: Optional[Dict[str, float]],
    target_encs,
    scheme: Optional[str],
) -> float:
    """Upper bound on posterior mass for *any* of several target groups.

    Sums the exact recorded mass for targets that are present, plus the
    shared unrecorded remainder *once* if at least one target is absent from
    the distribution -- not once per absent target, which would double-count
    the same unknown pool of probability mass.
    """
    if not target_encs:
        return 0.0
    enc = encoded_distribution(distribution, scheme)
    known = 0.0
    any_unrecorded = False
    for t in target_encs:
        if t in enc:
            known += enc[t]
        else:
            any_unrecorded = True
    remainder = max(0.0, 1.0 - sum(enc.values()))
    return known + (remainder if any_unrecorded else 0.0)


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


def _changed_side_walk(
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    mrca_id: Optional[int],
    derived_enc: Optional[str],
    scheme: Optional[str],
) -> Tuple[List[int], List[float], bool]:
    """The labeling-invariant part of a changed-side score.

    Returns ``(path_ids, cumprod, contam_hop1)`` where ``path_ids`` is the node
    walk (parent-of-MRCA .. root), ``cumprod[i] = product_{0..i} (1 - P(derived))``
    with P(derived) the worst-case bound, and ``contam_hop1`` is whether the
    node directly above the MRCA already has the derived state as its mode.

    None of these depend on the phenotype labeling or on which hypothesis is
    being scored — only on ``(mrca_id, site, scheme, derived_enc)`` — so the
    FOP null and multi-hypothesis observed runs memoise this once per
    ``(mrca_id, derived_enc)`` within a ``(site, scheme)`` and index into
    ``cumprod`` by each hypothesis's LAC cut (see :func:`_apply_changed_stop`).
    """
    path = path_to_root_ids(node_index, mrca_id)
    cumprod: List[float] = []
    running = 1.0
    contam_hop1 = False
    for k, node_id in enumerate(path, start=1):
        dist = node_dist(per_node_dist, node_id)
        p_der = worst_case_group_probability(dist, derived_enc, scheme)
        running *= max(0.0, 1.0 - p_der)
        cumprod.append(running)
        if k == 1 and derived_enc and modal_encoded(dist, scheme) == derived_enc:
            contam_hop1 = True
    return path, cumprod, contam_hop1


def _apply_changed_stop(
    path: List[int], cumprod: List[float], contam_hop1: bool,
    stop_at_id: Optional[int],
) -> Tuple[float, bool]:
    """Index a precomputed changed-side walk at one hypothesis's LAC cut.

    Bit-identical to the inline loop in :func:`side_path_score`: the walk stops
    *before* ``stop_at_id`` (the shared LAC), so ``n_walk`` nodes are scored and
    the score is ``cumprod[n_walk - 1]`` (or 1.0 when nothing is walked).
    """
    if not path:
        return EMPTY_PATH_SCORE, False
    n_walk = len(path)
    if stop_at_id is not None:
        for j, nid in enumerate(path):
            if nid == stop_at_id:
                n_walk = j
                break
    if n_walk == 0:
        return 1.0, False  # count == 0: no private nodes; hop+1 never reached
    return max(0.0, min(1.0, cumprod[n_walk - 1])), contam_hop1


def side_path_score(
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    mrca_id: Optional[int],
    ancestral_enc: str,
    derived_enc: Optional[str],
    scheme: Optional[str],
    is_changed: bool,
    stop_at_id: Optional[int] = None,
    walk_cache: Optional[Dict[Any, Any]] = None,
    cache_scope: Optional[Tuple[Any, Optional[str]]] = None,
) -> Tuple[float, bool]:
    """Score private isolation (changed) or global conservation (conserved).

    For a changed side (is_changed=True), walks the private segment from the
    parent of the MRCA up to stop_at_id (exclusive) and computes the product
    of (1 - P(derived)), where P(derived) is the worst-case bound (exact if
    the derived residue is the node's recorded state, otherwise the total
    unrecorded mass) rather than a naive lookup that silently reads 0 for any
    residue that is not the recorded state. If the private segment is empty
    (count == 0, e.g. sibling merge), returns 1.0 (no private contamination).

    For a conserved side (is_changed=False), walks the entire path to the root
    and computes the unweighted mean of P(ancestral) (no stop_at_id).

    ``walk_cache`` + ``cache_scope`` (``(site_key, scheme)``): when both are
    supplied the labeling-invariant walk is memoised, so scoring N FOP
    hypotheses for one (site, scheme) walks each distinct pair once instead of
    N times. Results are bit-identical to the uncached path.
    """
    if walk_cache is not None and cache_scope is not None:
        if is_changed:
            key = (cache_scope[0], cache_scope[1], "chg", mrca_id, derived_enc)
            entry = walk_cache.get(key)
            if entry is None:
                entry = _changed_side_walk(
                    node_index, per_node_dist, mrca_id, derived_enc, scheme
                )
                walk_cache[key] = entry
            return _apply_changed_stop(entry[0], entry[1], entry[2], stop_at_id)
        else:
            key = (cache_scope[0], cache_scope[1], "cons", mrca_id, ancestral_enc)
            if key not in walk_cache:
                walk_cache[key] = _conserved_side_score(
                    node_index, per_node_dist, mrca_id, ancestral_enc, scheme
                )
            return walk_cache[key], False

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
            p_der = worst_case_group_probability(dist, derived_enc, scheme)
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


def _conserved_side_score(
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    mrca_id: Optional[int],
    ancestral_enc: str,
    scheme: Optional[str],
) -> float:
    """Conserved-side score (mean P(ancestral) over the full MRCA→root path).

    Split out so it can be memoised; identical arithmetic to the
    ``is_changed=False`` branch of :func:`side_path_score`.
    """
    path = path_to_root_ids(node_index, mrca_id)
    if not path:
        return EMPTY_PATH_SCORE
    total_anc = 0.0
    for node_id in path:
        total_anc += group_probability(node_dist(per_node_dist, node_id), ancestral_enc, scheme)
    return max(0.0, min(1.0, total_anc / len(path)))


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


# ── Shared combinatorics ──────────────────────────────────────────────────────
def _p_at_least_2(p_list: List[float]) -> float:
    """Exact P(>= 2 successes) via inclusion-exclusion over independent Bernoullis."""
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


# ── Position-level aggregation ───────────────────────────────────────────────
def compute_asr_path_score(
    pair_details: Optional[List[Dict[str, Any]]],
    per_node_dist: Dict[Any, Dict[str, float]],
    node_index: Dict[int, Any],
    scheme: Optional[str],
    is_conserved_meta: bool,
    conserved_pair: Optional[str],
    diversity_floor: float = 0.75,
    walk_cache: Optional[Dict[Any, Any]] = None,
    site_key: Optional[Any] = None,
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
        Dict with ``asr_path_score`` (position level, 0-1), ``replication``
        (independence * core), ``strength`` (diversity_mult *
        derived_agreement), ``independence``, ``mrca_diversity``,
        ``derived_agreement``, ``conservation_gate``, ``core``,
        ``pair_scores`` ({pair_id: score}), ``pair_contaminated``,
        ``conserved_pair_scores`` ({pair_id: conservation-to-root}) and
        ``conserved_pair_nodes`` ({pair_id: mrca_node_id}) for the conserved
        pairs folded into ``conservation_gate`` (both empty when none), and
        ``pair_ancestral`` / ``pair_derived_top`` / ``pair_derived_bot``
        ({pair_id: raw AA}) — the un-encoded ancestral and per-side derived
        residues of each changed pair, for the FOP harvest-wide, per-scheme
        derived_agreement rebuild (POINT 3).
    """
    pairs = pair_details or []
    n_pairs = len(pairs)
    conserved_ids = (
        set(parse_conserved_ids(conserved_pair, n_pairs)) if is_conserved_meta else set()
    )

    pair_contaminated: Dict[int, bool] = {}
    conserved_conservations: List[float] = []  # conservation-to-root of conserved pairs
    # Per-conserved-pair record, analogous to pair_scores / pair_contaminated for
    # changed pairs. Lets the FOP domain-pooler (fop_pool.R / fop_pool.py)
    # reconstruct conservation_gate from the DISTINCT conserved pairs shared
    # across hypotheses instead of averaging already-transformed per-hypothesis
    # gates. Keyed by pair_id.
    conserved_pair_scores: Dict[int, float] = {}
    conserved_pair_nodes: Dict[int, Any] = {}
    # Derived (encoded) residues of changed tips, kept per phenotype side so
    # within-side divergence (the non-convergent case) is assessed separately
    # from a pair changing on *both* sides (a strong convergent signal). Each
    # entry also carries the pair id so derived_agreement can count pairs per
    # residue without a second pass over `pairs`.
    derived_by_side: Dict[str, List[Tuple[int, str]]] = {
        "top_tip_mode": [], "bottom_tip_mode": [],
    }
    derived_pool: set = set()  # all encoded derived states across changed pairs

    # Raw (un-encoded) residues per changed pair, kept so the FOP domain-pooler
    # (fop_pool.R / fop_pool.py POINT 3) can recompute derived_agreement
    # HARVEST-WIDE and PER SCHEME over the pooled, node-deduplicated changed-pair
    # set — a position that is unanimous within each hypothesis but split BETWEEN
    # hypotheses (US: V/I/L) then gets a low harvest-wide da under US and 1.0
    # under a scheme that co-encodes those residues, automatically. Keyed by
    # pair_id. ``pair_derived_top`` / ``pair_derived_bot`` only carry a residue
    # for the side that actually changed (empty otherwise).
    pair_ancestral: Dict[int, Optional[str]] = {}
    pair_derived_top: Dict[int, str] = {}
    pair_derived_bot: Dict[int, str] = {}

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
                walk_cache=walk_cache, cache_scope=(site_key, scheme),
            )
            conserved_conservations.append(cons)
            conserved_pair_scores[pid] = cons
            conserved_pair_nodes[pid] = mrca_id
            continue

        sides: List[Tuple[str, str]] = []
        for side_key in ("top_tip_mode", "bottom_tip_mode"):
            raw_tip = pair.get(side_key)
            tip_enc = encode_aa(raw_tip, scheme)
            if tip_enc is None or tip_enc == anc_enc:
                continue  # missing or conserved side → not scored
            sides.append((side_key, tip_enc))
            derived_by_side[side_key].append((pid, tip_enc))
            derived_pool.add(tip_enc)
            raw_tip_u = str(raw_tip).strip().upper()
            if side_key == "top_tip_mode":
                pair_derived_top[pid] = raw_tip_u
            else:
                pair_derived_bot[pid] = raw_tip_u

        if sides:
            fs = pair.get("focal_state")
            pair_ancestral[pid] = str(fs).strip().upper() if fs else None

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
            "replication": 0.0,
            "strength": 0.0,
            "independence": 1.0,
            "mrca_diversity": 0.0,
            "derived_agreement": 0.0,
            "conservation_gate": conservation_gate,
            "core": 0.0,
            "pair_scores": {},
            "top_pair_scores": {},
            "bottom_pair_scores": {},
            "pair_contaminated": {},
            "conserved_pair_scores": conserved_pair_scores,
            "conserved_pair_nodes": conserved_pair_nodes,
            "pair_ancestral": pair_ancestral,
            "pair_derived_top": pair_derived_top,
            "pair_derived_bot": pair_derived_bot,
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
    # Each pair's own MRCA plus its private-segment node ids (parent-of-MRCA
    # up to, excluding, its nearest LAC) -- used below by mrca_diversity as
    # the set of places another pair's background state might be found. The
    # MRCA itself must be included: for a "sibling merge" (LAC is directly
    # the MRCA's parent) the private segment is empty, and without the MRCA
    # as a fallback search point, two pairs with IDENTICAL MRCA states would
    # wrongly read as maximally diverse (nothing to find them in) rather than
    # maximally similar.
    diversity_search_nodes: Dict[int, List[int]] = {}
    for c in changed:
        # Nearest LAC = deepest merge point on this pair's root-path. The private
        # segment (below it) is the pair's own, judged by the core walk; the
        # shared segment (at/above it) is judged by the independence product.
        full_path = path_to_root_ids(node_index, c["mrca_id"])
        stop_at = next((n for n in full_path if n in lac_nodes), None)
        private_segment = (
            full_path if stop_at is None else full_path[: full_path.index(stop_at)]
        )
        diversity_search_nodes[c["pid"]] = [c["mrca_id"]] + private_segment
        side_scores: List[float] = []
        pair_contam = False
        for side_key, tip_enc in c["sides"]:
            score, contam = side_path_score(
                node_index, per_node_dist, c["mrca_id"], c["anc_enc"], tip_enc,
                scheme, is_changed=True, stop_at_id=stop_at,
                walk_cache=walk_cache, cache_scope=(site_key, scheme),
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
    # For each LAC node, (1 - P(any residue in derived_pool present there)).
    # Product across LAC nodes: if any merge point already carries a derived
    # state, the pairs below it did not converge independently — they
    # inherited it — and the score craters. P(any of derived_pool) uses the
    # same worst-case bound as core (exact for a recorded residue, the
    # unrecorded remainder added once -- not once per candidate residue,
    # which would double-count the same unknown mass -- for however many of
    # derived_pool aren't the node's recorded state).
    independence = 1.0
    for lac in lac_nodes:
        dist = node_dist(per_node_dist, lac)
        p_derived = worst_case_any_group_probability(dist, derived_pool, scheme)
        independence *= max(0.0, 1.0 - p_derived)

    # Parallel axis (continuous): did the changed pairs' private segments (the
    # same MRCA->LAC walks core just scored) pass through genuinely different
    # ancestral backgrounds, or does one pair's own MRCA state show up
    # somewhere along another pair's walk (suggesting a shared background)?
    # For every pair of changed pairs, check both directions -- does pair A's
    # MRCA state appear anywhere in pair B's private segment, and does pair
    # B's appear anywhere in A's -- using the same worst-case bound at every
    # node visited, combined into "found somewhere" via the same P(>=1) union
    # logic as core's P(>=2): the chance a state is absent at every node in a
    # segment, flipped. The two directions are then unioned into one
    # "these share a background" probability for that pair-comparison, and
    # the position's diversity is the mean over every pairwise comparison
    # among changed pairs -- the same aggregation the old MRCA-only version
    # used, just built from richer, multi-node comparisons.
    def _found_probability(
        target_enc: Optional[str], segment_nodes: List[int],
        owner_mrca_id: Optional[int] = None,
    ) -> float:
        if target_enc is None or not segment_nodes:
            return 0.0
        # diversity_search_nodes[pid] is always [mrca_id] + a prefix of that
        # pair's MRCA->root path, so the "absent at every node" product is a
        # prefix of the cumulative product over [mrca_id] + full_path. Memoise
        # that once per (site, scheme, owner_mrca, target_enc) and slice by
        # segment length -- same labeling-invariance argument as core.
        if walk_cache is not None and site_key is not None and owner_mrca_id is not None:
            key = (site_key, scheme, "found", owner_mrca_id, target_enc)
            cp = walk_cache.get(key)
            if cp is None:
                full = [owner_mrca_id] + path_to_root_ids(node_index, owner_mrca_id)
                cp = []
                running = 1.0
                for nid in full:
                    running *= max(0.0, 1.0 - worst_case_group_probability(
                        node_dist(per_node_dist, nid), target_enc, scheme))
                    cp.append(running)
                walk_cache[key] = cp
            n = len(segment_nodes)
            if n == 0 or not cp:
                return 0.0
            return max(0.0, min(1.0, 1.0 - cp[min(n, len(cp)) - 1]))
        p_absent = 1.0
        for node_id in segment_nodes:
            dist = node_dist(per_node_dist, node_id)
            p_absent *= max(0.0, 1.0 - worst_case_group_probability(dist, target_enc, scheme))
        return max(0.0, min(1.0, 1.0 - p_absent))

    if len(changed) >= 2:
        pairwise_diversities: List[float] = []
        for ca, cb in combinations(changed, 2):
            a_state = modal_encoded(node_dist(per_node_dist, ca["mrca_id"]), scheme)
            b_state = modal_encoded(node_dist(per_node_dist, cb["mrca_id"]), scheme)
            p_a_in_b = _found_probability(
                a_state, diversity_search_nodes.get(cb["pid"], []), cb["mrca_id"])
            p_b_in_a = _found_probability(
                b_state, diversity_search_nodes.get(ca["pid"], []), ca["mrca_id"])
            shared = 1.0 - (1.0 - p_a_in_b) * (1.0 - p_b_in_a)
            pairwise_diversities.append(1.0 - shared)
        diversity = sum(pairwise_diversities) / len(pairwise_diversities)
    else:
        diversity = 0.5  # single changed pair: parallelism not assessable
    diversity_mult = diversity_floor + (1.0 - diversity_floor) * diversity

    # Convergence axis: within each phenotype side that has at least two
    # changed pairs, what fraction land on the position's single most common
    # (plurality) derived residue? A side with exactly one changed pair has
    # nothing to agree or disagree with and is excluded from the average
    # entirely — "agreement" is a >=2-observer concept, and crediting a lone
    # change as trivial (1/1) agreement would dilute a real disagreement
    # elsewhere on the position. Computed from plain pair counts: this axis is
    # a fact about observed tip residues (no posterior involved), kept
    # independent of the ancestral-uncertainty question `core`/`independence`
    # already own.
    concentration_by_side: Dict[str, float] = {}
    for side_key, pid_residues in derived_by_side.items():
        if len(pid_residues) < 2:
            continue
        count_by_residue: Dict[str, int] = {}
        for _pid, r in pid_residues:
            count_by_residue[r] = count_by_residue.get(r, 0) + 1
        total = len(pid_residues)
        plurality_r = max(count_by_residue, key=count_by_residue.get)
        concentration_by_side[side_key] = count_by_residue[plurality_r] / total

    if concentration_by_side:
        derived_agreement = sum(concentration_by_side.values()) / len(concentration_by_side)
    else:
        derived_agreement = 1.0  # no side has >=2 changed pairs to assess agreement on

    # ── Combinatorial Core Score: P(>= 2 independent changes) ──────────────────
    # Computed *per phenotype side* and then combined via union. This is the
    # critical directional check: a pair changing only on the top side and a
    # different pair changing only on the bottom side are changes in opposite
    # directions — not convergence. Without the per-side split they would both
    # enter p_list and produce P(≥2) ≈ 1 even though neither direction has ≥2
    # independent events.
    core_top = _p_at_least_2(list(top_pair_scores.values()))
    core_bottom = _p_at_least_2(list(bottom_pair_scores.values()))
    core = 1.0 - (1.0 - core_top) * (1.0 - core_bottom)

    # replication: P(shared ancestor clean AND >=2 private-clean events). A
    # genuine joint probability, not two independent checks — private-segment
    # nodes and LAC nodes are disjoint, conditionally-independent parts of the
    # tree, so P(A) * P(B|A) = P(A ∩ B).
    replication = independence * core
    # strength: two independent discounts on the quality of a confirmed
    # replication event — did the lineages start from different places, and
    # did they land on the same place.
    strength = diversity_mult * derived_agreement

    asr_path_score = max(0.0, min(1.0, replication * strength * conservation_gate))

    return {
        "asr_path_score": asr_path_score,
        "replication": replication,
        "strength": strength,
        "independence": independence,
        "mrca_diversity": diversity,
        "derived_agreement": derived_agreement,
        "conservation_gate": conservation_gate,
        "core": core,
        "pair_scores": pair_scores,
        "top_pair_scores": top_pair_scores,
        "bottom_pair_scores": bottom_pair_scores,
        "pair_contaminated": pair_contaminated,
        "conserved_pair_scores": conserved_pair_scores,
        "conserved_pair_nodes": conserved_pair_nodes,
        "pair_ancestral": pair_ancestral,
        "pair_derived_top": pair_derived_top,
        "pair_derived_bot": pair_derived_bot,
    }
