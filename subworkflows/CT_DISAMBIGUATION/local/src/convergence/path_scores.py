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
merges with another changed pair's — the LCA), and the **LCA nodes**
themselves. Nodes above the LCA are not examined. At each node visited, we
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
was not acquired. It is scored by conservation-to-root and reported in
``conserved_pair_scores`` / ``conserved_pair_nodes`` (latent — no longer a
score multiplier in T1; see Aggregation).

Aggregation (scoring_v2 T3a — per-side pairwise core)
----------------------------------------------------
One ``core_s`` per ``(Gene, Position, scheme, side)``, ``s ∈ {top, bottom}``,
built as a pairwise aggregation that folds in the former ``independence`` and
``derived_agreement`` axes per pair instead of as global multipliers (see
``docs/scoring_v2_T3_core_pareado.md``):

* ``D_s`` — the design of the side: pairs that changed on ``s`` (``P_s``) plus
  the conserved-metadata pairs (``C``, both sides). A pair that changed *only on
  the other side* is **not** in ``D_s``.
* ``s_c^s`` — private-segment isolation of pair ``c`` (unchanged primitive:
  ``∏ (1 − P_wc(derived))`` from the parent of the MRCA up to the nearest
  same-side LCA merge point).
* ``contrib(c, d) = s_c^s · s_d^s · [enc(der_c) == enc(der_d)]
  · (1 − P_wc(enc(der_c) @ LCA(mrca_c, mrca_d)))`` — probability that ``c`` and
  ``d`` form a genuine convergent pair: both changed for real, at the same
  encoded residue, with a clean shared ancestor. ``independence`` now enters
  here, once per pair, keyed to that pair's own merge node and residue.
* ``score_c = noisy_or_{d ∈ partners(c)} contrib(c, d)`` — best evidence that
  ``c`` converged with at least one same-residue partner (``0`` with no partner).
* ``core_s = ( Σ_{c ∈ D_s} score_c ) / |D_s|`` — mean over the whole design.
  Conserved pairs and orphan-residue pairs contribute ``score_c = 0`` but count
  in the denominator, so a design of ``n`` pairs where only ``k`` converge scores
  ``≈ k/n``.
* ``asr_path_score[side = s] = clamp01(core_s)``. Sides never recombine — a
  "both" position is two independent rows.

``agree(c, d)`` is hard 0/1 on the encoded residue; the biochemical gradient is
supplied later by averaging the five schemes in ``scoring_compute.R`` §2g, not by
a soft rule here. ``convergence_type`` is derived from
``(agree_num_s, agree_den_s)`` (see :func:`_convergence_type`), replacing the
separate categorical classifier.

Conserved pairs are still walked (``conserved_pair_scores`` /
``conserved_pair_nodes`` flow out for diagnostics and the FOP pooler); they no
longer multiply anything — they are ``score_c = 0`` members of ``D_s``.

Return shape: ``native_side_split=True`` → ``{"top": {...}, "bottom": {...}}``.
Default (feature flag ``--native_side_split`` off) → one flat T1-shaped dict
whose ``asr_path_score`` is ``max(core_top, core_bottom)`` so ``disambiguate_
single`` keeps emitting a single row until T3b.

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
    the three call sites (per-node walk, LCA product, MRCA-diversity) consistent
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

    Used to locate the **LCA** (lowest common ancestor) nodes of a set of pair
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
    ``cumprod`` by each hypothesis's LCA cut (see :func:`_apply_changed_stop`).
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
    """Index a precomputed changed-side walk at one hypothesis's LCA cut.

    Bit-identical to the inline loop in :func:`side_path_score`: the walk stops
    *before* ``stop_at_id`` (the shared LCA), so ``n_walk`` nodes are scored and
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
            break  # reached the shared LCA; stop private segment walk

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


def noisy_or(probs) -> float:
    """``1 - ∏(1 - p_i)`` — probability at least one of independent events fires.

    Empty input returns ``0.0`` (a pair with no same-residue partner has no
    convergence evidence). Each ``p_i`` is clamped to ``[0, 1]`` first.
    """
    acc = 1.0
    for p in probs:
        acc *= (1.0 - max(0.0, min(1.0, float(p))))
    return max(0.0, min(1.0, 1.0 - acc))


def aggregate_core_side(
    participants: List[Dict[str, Any]],
    n_conserved: int,
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    scheme: Optional[str],
    iso_override: Optional[Dict[Any, float]] = None,
    walk_cache: Optional[Dict[Any, Any]] = None,
    cache_scope: Optional[Tuple[Any, Optional[str]]] = None,
) -> Dict[str, Any]:
    """The per-side pairwise ``core_s`` aggregation (T3-doc §5), factored out.

    One ``core_s`` over ``D_s`` = ``participants`` ∪ ``n_conserved`` conserved
    slots. Each participant is ``{"pid", "mrca_id", "der_enc", "anc_enc"}`` — a
    pair that changed on this side to encoded residue ``der_enc``.

    * ``iso[pid] = s_c^s`` — private-segment isolation, walk stopped at the
      nearest same-side LCA merge point. Computed here unless ``iso_override``
      supplies it (the FOP pooler pools ``s_c^s`` per domain across hypotheses
      first, then passes it in keyed by ``pid``).
    * ``contrib(c, d) = iso_c · iso_d · [der_c == der_d]
      · (1 − P_wc(der_c @ LCA(mrca_c, mrca_d)))`` over same-residue participant
      pairs; ``score_c = noisy_or`` of a participant's contribs; ``core_s`` =
      mean of ``score_c`` over ``|participants| + n_conserved``.

    Returns the numeric per-side fields (no position-level conserved/derived
    maps — the caller attaches those). Bit-identical to the pre-refactor
    ``_side_result`` inner block.
    """
    P = participants
    n_part = len(P)
    n = n_part + max(0, int(n_conserved))

    # L_s — merge points of THIS side's participants only (T3-doc §5.2).
    lca_s: set = set()
    for a, b in combinations(P, 2):
        lca = find_lca(node_index, a["mrca_id"], b["mrca_id"])
        if lca is not None:
            lca_s.add(lca)

    # s_c^s — private-segment isolation; walk stops at the nearest L_s node.
    iso: Dict[Any, float] = {}
    contam: Dict[Any, bool] = {}
    for p in P:
        if iso_override is not None and p["pid"] in iso_override:
            iso[p["pid"]] = float(iso_override[p["pid"]])
            contam[p["pid"]] = False
            continue
        full_path = path_to_root_ids(node_index, p["mrca_id"])
        stop_at = next((x for x in full_path if x in lca_s), None)
        sc, ct = side_path_score(
            node_index, per_node_dist, p["mrca_id"],
            p["anc_enc"], p["der_enc"], scheme,
            is_changed=True, stop_at_id=stop_at,
            walk_cache=walk_cache, cache_scope=cache_scope,
        )
        iso[p["pid"]] = sc
        contam[p["pid"]] = ct

    # contrib(c, d) over unordered same-residue participant pairs.
    contribs: Dict[Any, List[float]] = {p["pid"]: [] for p in P}
    for a, b in combinations(P, 2):
        if a["der_enc"] != b["der_enc"]:
            continue  # agree == 0 → contrib == 0
        lca_ab = find_lca(node_index, a["mrca_id"], b["mrca_id"])
        p_shared = worst_case_group_probability(
            node_dist(per_node_dist, lca_ab), a["der_enc"], scheme
        )
        k = iso[a["pid"]] * iso[b["pid"]] * max(0.0, 1.0 - p_shared)
        contribs[a["pid"]].append(k)
        contribs[b["pid"]].append(k)

    partner_scores = {pid: noisy_or(v) for pid, v in contribs.items()}
    core_s = (sum(partner_scores.values()) / n) if n else 0.0
    core_s = max(0.0, min(1.0, core_s))

    by_res: Dict[str, int] = {}
    for p in P:
        by_res[p["der_enc"]] = by_res.get(p["der_enc"], 0) + 1
    agree_den = n_part
    agree_num = max(by_res.values()) if by_res else 0
    concentration = (agree_num / agree_den) if agree_den else 0.0

    return {
        "asr_path_score": core_s,
        "core": core_s,
        "n_pairs_side": n,
        "n_participating": n_part,
        "n_conserved": max(0, int(n_conserved)),
        "derived_agreement": concentration,
        "agree_num": agree_num,
        "agree_den": agree_den,
        "convergence_type": _convergence_type(agree_num, agree_den),
        "pair_scores": dict(iso),
        "pair_partner_scores": partner_scores,
        "pair_contaminated": contam,
        # Per-participant structure so a downstream pooler (the FOP null,
        # SC2b) can rebuild ``participants`` for a fresh ``aggregate_core_side``
        # call over the node-deduped union of pairs across hypotheses.
        "pair_mrca": {p["pid"]: p["mrca_id"] for p in P},
        "pair_der_enc": {p["pid"]: p["der_enc"] for p in P},
        "pair_anc_enc": {p["pid"]: p["anc_enc"] for p in P},
    }


def _convergence_type(agree_num: int, agree_den: int) -> str:
    """Categorical label from the agreement numerator/denominator (T3-doc §7).

    ``agree_den`` = number of pairs that changed on the side; ``agree_num`` =
    size of the largest single-residue group among them. Derived here (not from
    the separate ``convergence.py`` classifier) so T4a can retire that path.
    """
    if agree_den >= 2 and agree_num >= 2:
        return "convergent"
    if agree_den >= 2:
        return "divergent"
    if agree_den == 1:
        return "single"
    return "no_change"


# ── Position-level aggregation ───────────────────────────────────────────────
def compute_asr_path_score(
    pair_details: Optional[List[Dict[str, Any]]],
    per_node_dist: Dict[Any, Dict[str, float]],
    node_index: Dict[int, Any],
    scheme: Optional[str],
    is_conserved_meta: bool,
    conserved_pair: Optional[str],
    walk_cache: Optional[Dict[Any, Any]] = None,
    site_key: Optional[Any] = None,
    native_side_split: bool = False,
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
            enter ``D_s`` on both sides as ``score_c = 0`` members of the mean's
            denominator (T3-doc §5.1) and are still reported in
            ``conserved_pair_scores`` / ``conserved_pair_nodes``.
        native_side_split: when True return ``{"top": {...}, "bottom": {...}}``
            — one row per ``(Gene, Position, scheme, side)`` (T3-doc §4). When
            False (default, feature flag ``--native_side_split`` off) collapse to
            one flat T1-shaped dict whose ``asr_path_score`` is the stronger
            side, so ``disambiguate_single`` keeps emitting a single row until
            T3b wires the two-row plumbing.

    Returns:
        ``native_side_split=False``: flat dict — ``asr_path_score`` (=
        ``max(core_top, core_bottom)``), ``core`` (alias), ``core_top`` /
        ``core_bottom``, ``derived_agreement`` (diagnostic concentration),
        ``pair_scores`` / ``top_pair_scores`` / ``bottom_pair_scores`` (per-pair
        private isolation ``s_c^s``), ``pair_partner_scores`` ({pid: score_c}),
        ``pair_contaminated``, ``conserved_pair_scores`` / ``conserved_pair_nodes``,
        ``pair_ancestral`` / ``pair_derived_top`` / ``pair_derived_bot``.

        ``native_side_split=True``: ``{"top": <row>, "bottom": <row>}`` where
        each ``<row>`` carries ``asr_path_score`` = ``core_s`` = ``clamp01`` of
        the mean of ``score_c`` (noisy-OR of a pair's same-residue partners) over
        ``D_s`` = participants ∪ conserved-metadata pairs, plus ``n_pairs_side``,
        ``n_participating``, ``n_conserved``, ``agree_num`` / ``agree_den``,
        ``convergence_type``, and the same pair-level / conserved-pair maps.
    """
    pairs = pair_details or []
    n_pairs = len(pairs)
    conserved_ids = (
        set(parse_conserved_ids(conserved_pair, n_pairs)) if is_conserved_meta else set()
    )

    pair_contaminated: Dict[int, bool] = {}
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
    # walk needs the LCA merge points computed below to know where to stop).
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

    # ── T3 per-side pairwise core (T3-doc §5) ────────────────────────────────
    # One ``core_s`` per (Gene, Position, scheme, side). For side s:
    #
    #   D_s   = P_s ∪ C            P_s = pairs that changed on s
    #                              C   = conserved-metadata pairs (both sides)
    #   L_s   = { LCA(mrca_a, mrca_b) : a, b ∈ P_s }        merge points of s
    #   s_c^s = ∏ (1 − P_wc(der_c @ k)) over the private segment below stop_c
    #   contrib(c,d) = s_c^s · s_d^s · [enc(der_c)==enc(der_d)]
    #                  · (1 − P_wc(enc(der_c) @ LCA(mrca_c, mrca_d)))
    #   score_c = noisy_or_{d ∈ partners(c)} contrib(c,d)   0 if no partner
    #   core_s  = ( Σ_{c ∈ D_s} score_c ) / |D_s|           0 if |D_s| = 0
    #
    # Conserved pairs and pairs that landed on an orphan residue contribute
    # score_c = 0 but still count in the denominator — a design of n pairs where
    # only k converge scores ≈ k/n (T3-doc §6.1). Lados never mix.
    anc_by_pid: Dict[int, str] = {c["pid"]: c["anc_enc"] for c in changed}
    conserved_present: List[int] = list(conserved_pair_scores.keys())

    # P_s per side: {pid, mrca_id, der_enc} for each changed side of each pair.
    participants: Dict[str, List[Dict[str, Any]]] = {
        "top_tip_mode": [], "bottom_tip_mode": [],
    }
    for c in changed:
        for side_key, tip_enc in c["sides"]:
            participants[side_key].append(
                {"pid": c["pid"], "mrca_id": c["mrca_id"], "der_enc": tip_enc,
                 "anc_enc": anc_by_pid[c["pid"]]}
            )

    def _side_result(side_key: str) -> Dict[str, Any]:
        agg = aggregate_core_side(
            participants[side_key], len(conserved_present),
            node_index, per_node_dist, scheme,
            walk_cache=walk_cache, cache_scope=(site_key, scheme),
        )
        return {
            **agg,
            "conserved_pair_scores": conserved_pair_scores,
            "conserved_pair_nodes": conserved_pair_nodes,
            "pair_ancestral": pair_ancestral,
            "pair_derived_top": pair_derived_top,
            "pair_derived_bot": pair_derived_bot,
        }

    top = _side_result("top_tip_mode")
    bottom = _side_result("bottom_tip_mode")

    if native_side_split:
        # T3-doc §4: one row per (Gene, Position, scheme, side). Downstream
        # (T3b) decides whether each row is emitted (per change_side).
        return {"top": top, "bottom": bottom}

    # ── Transitional flat return — feature flag --native_side_split OFF ───────
    # main and the single-row null still consume one flat T1-shaped dict. The
    # score collapses to the stronger side (a "both" position is not split into
    # two rows until T3b). ``independence`` / ``replication`` are dropped
    # (T3-doc §13); consumers read them with ``.get`` defaults.
    for pid, ct in top["pair_contaminated"].items():
        pair_contaminated[pid] = pair_contaminated.get(pid, False) or ct
    for pid, ct in bottom["pair_contaminated"].items():
        pair_contaminated[pid] = pair_contaminated.get(pid, False) or ct

    flat_pair_scores: Dict[int, float] = {}
    for c in changed:
        vals = [s[c["pid"]] for s in (top["pair_scores"], bottom["pair_scores"])
                if c["pid"] in s]
        if vals:
            flat_pair_scores[c["pid"]] = sum(vals) / len(vals)

    flat_partner: Dict[int, float] = {}
    for src in (top["pair_partner_scores"], bottom["pair_partner_scores"]):
        for pid, v in src.items():
            flat_partner[pid] = max(flat_partner.get(pid, 0.0), v)

    qual = [s["derived_agreement"] for s in (top, bottom) if s["agree_den"] >= 2]
    derived_agreement = (sum(qual) / len(qual)) if qual else 1.0

    core = max(top["asr_path_score"], bottom["asr_path_score"])

    return {
        "asr_path_score": core,
        "core": core,
        "core_top": top["asr_path_score"],
        "core_bottom": bottom["asr_path_score"],
        "derived_agreement": derived_agreement,
        "pair_scores": flat_pair_scores,
        "top_pair_scores": dict(top["pair_scores"]),
        "bottom_pair_scores": dict(bottom["pair_scores"]),
        "pair_partner_scores": flat_partner,
        "pair_contaminated": pair_contaminated,
        "conserved_pair_scores": conserved_pair_scores,
        "conserved_pair_nodes": conserved_pair_nodes,
        "pair_ancestral": pair_ancestral,
        "pair_derived_top": pair_derived_top,
        "pair_derived_bot": pair_derived_bot,
    }
