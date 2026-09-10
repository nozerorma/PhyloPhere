#!/usr/bin/env python3
"""CAAS convergence scoring on the Voronoi domain (scoring_v2 core v3).

This module scores one CAAS position as a **mean over its K fixed Voronoi
domains**. Each domain is a candidate fg/bg pair extracted by the canonical Dunn
selection of the reference hypothesis; the domain index *is* the unit of the
score (there is no per-pair private-segment walk any more — that machinery, and
the ``walk_cache`` / ``L_s`` / ``iso_override`` scaffolding around it, was
removed in core v3).

Conceptual model
----------------
For a domain ``d`` we know its MRCA node, the modal ancestral state there
(``focal_state``, encoded in the active scheme as ``anc_enc``), and the observed
fg tip residue per phenotype side. A domain **changed on side ``s``** iff its tip
residue is defined and differs from ``anc_enc``.

Two changed domains ``d_i``, ``d_j`` on the same side that landed on the *same*
encoded residue form a contrast:

    p_shared      = worst_case_any_group_probability(dist @ LCA(mrca_i, mrca_j),
                                                     {der_i, der_j}, scheme)
    contrib(i, j) = 1 - p_shared

The distribution is read **only at the LCA node** — no walk. ``score_d`` is the
``noisy_or`` of a domain's contribs with its same-residue partners (``0`` with no
partner). A domain that did not converge is simply ``score_d = 0``; a domain
without a valid reconstruction (no modal ASR at its MRCA, unmappable tip, or
missing ``focal_state``) is not "changed" and scores ``0``, but still occupies a
slot in ``domain_meta`` so the pooler keeps it in the denominator.

``agree`` is hard 0/1 on the encoded residue; the biochemical gradient is
supplied later by averaging the five schemes in ``scoring_compute.R`` section 2g,
not by a soft rule here. ``convergence_type`` is derived from
``(agree_num, agree_den)`` (see :func:`_convergence_type`) for the single
hypothesis; the harvest-wide label is recomputed in ``fop_pool.pool_domains``.

The PSS weights do **not** enter this layer — they weight domains only in the
pooling step (``fop_pool.pool_domains``). ``compute_domain_scores`` returns one
per-hypothesis record; ``pool_domains`` collapses ``M >= 1`` of them (``M = 1``
degenerates to the plain PSS-weighted mean over the K domains).

Return shape of :func:`compute_domain_scores`::

    {"top": <side dict>, "bottom": <side dict>, "domain_meta": {d: {...}}}

See ``docs/scoring_v3_core.md`` section 2 and Appendix A for the exact contract
and the hand-worked golden arithmetic.

The module is intentionally free of PAML/IO dependencies so it can be unit
tested with synthetic trees and posteriors.

Author
------
Miguel Ramon Alonso — Evolutionary Genomics Lab, IBE-UPF
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Optional

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
    the call sites consistent and the ``per_node_dist`` key type honestly
    ``Any``.
    """
    if node_id is None:
        return {}
    return per_node_dist.get(node_id) or per_node_dist.get(str(node_id)) or {}


def path_to_root_ids(node_index: Dict[int, Any], mrca_id: Optional[int]) -> List[int]:
    """Node ids from the parent of the MRCA up to the root (MRCA excluded).

    Kept for :func:`find_lca`, which walks each node's ancestor chain to locate
    the merge point of two domain MRCAs.
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

    In core v3 this locates the shared ancestor of two changed domains' MRCAs;
    the ASR posterior read there is the only tree lookup the score makes.
    """
    if id_a is None or id_b is None:
        return None
    chain_a = [int(id_a)] + path_to_root_ids(node_index, id_a)  # node itself .. root
    set_b = set([int(id_b)] + path_to_root_ids(node_index, id_b))
    for n in chain_a:
        if n in set_b:
            return n
    return None


# ── Shared combinatorics ─────────────────────────────────────────────────────
def noisy_or(probs) -> float:
    """``1 - ∏(1 - p_i)`` — probability at least one of independent events fires.

    Empty input returns ``0.0`` (a domain with no same-residue partner has no
    convergence evidence). Each ``p_i`` is clamped to ``[0, 1]`` first.
    """
    acc = 1.0
    for p in probs:
        acc *= (1.0 - max(0.0, min(1.0, float(p))))
    return max(0.0, min(1.0, 1.0 - acc))


def _convergence_type(agree_num: int, agree_den: int) -> str:
    """Categorical label from the agreement numerator/denominator.

    ``agree_den`` = number of domains that changed on the side; ``agree_num`` =
    size of the largest single-encoded-residue group among them.
    """
    if agree_den >= 2 and agree_num >= 2:
        return "convergent"
    if agree_den >= 2:
        return "divergent"
    if agree_den == 1:
        return "single"
    return "no_change"


# ── Core v3: domain scoring ──────────────────────────────────────────────────
def score_domains_side(
    domains: List[Dict[str, Any]],
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    scheme: Optional[str],
) -> Dict[str, Any]:
    """Score one phenotype side from its **changed domains only**.

    ``domains`` = ``list[{d, mrca_id, anc_enc, der_enc}]``. Domains pair up when
    their ``der_enc`` match; for a pair the contribution is
    ``1 - worst_case_any_group_probability(dist @ LCA(mrca_a, mrca_b),
    {der_a, der_b}, scheme)`` (``dist`` read only at the LCA), and
    ``score_d = noisy_or`` over a domain's partner contributions (``0`` with no
    partner).

    Returns ``{domain_scores, agree_num, agree_den, n_changed}`` — ``agree_den``
    and ``n_changed`` are both ``len(domains)``.
    """
    contribs: Dict[Any, List[float]] = {dm["d"]: [] for dm in domains}
    for a, b in combinations(domains, 2):
        if a["der_enc"] != b["der_enc"]:
            continue  # agree == 0 → no contribution
        lca = find_lca(node_index, a["mrca_id"], b["mrca_id"])
        p_shared = worst_case_any_group_probability(
            node_dist(per_node_dist, lca), {a["der_enc"], b["der_enc"]}, scheme
        )
        contrib = max(0.0, 1.0 - p_shared)
        contribs[a["d"]].append(contrib)
        contribs[b["d"]].append(contrib)

    domain_scores = {d: noisy_or(v) for d, v in contribs.items()}

    by_res: Dict[str, int] = {}
    for dm in domains:
        by_res[dm["der_enc"]] = by_res.get(dm["der_enc"], 0) + 1
    agree_den = len(domains)
    agree_num = max(by_res.values()) if by_res else 0
    return {
        "domain_scores": domain_scores,
        "agree_num": agree_num,
        "agree_den": agree_den,
        "n_changed": len(domains),
    }


def compute_domain_scores(
    pair_details: Optional[List[Dict[str, Any]]],
    per_node_dist: Dict[Any, Dict[str, float]],
    node_index: Dict[int, Any],
    scheme: Optional[str],
) -> Dict[str, Any]:
    """Per-hypothesis domain scores for one ``(Gene, Position, scheme)``.

    Args:
        pair_details: list of per-domain dicts (``pair_id``, ``node_id``,
            ``focal_state``, ``top_tip_mode``, ``bottom_tip_mode``). ``pair_id``
            is the domain id ``d``; ``node_id`` its MRCA.
        per_node_dist: ``node_id -> {aa: posterior}`` for the focal site.
        node_index: ``node_id -> TreeNode`` (from :func:`build_node_index`).
        scheme: active grouping scheme (``"US"`` or a GS label).

    Returns:
        ``{"top": <side>, "bottom": <side>, "domain_meta": {d: meta}}``.

        ``<side>`` carries ``domain_scores`` (only changed domains),
        ``domain_der`` / ``domain_der_enc`` / ``domain_anc`` (raw tip, encoded
        tip, raw ancestral — only changed domains), ``agree_num`` / ``agree_den``
        / ``n_changed``, and ``convergence_type``.

        ``domain_meta`` carries **all K domains** (every ``pair_detail``),
        including those without a reconstruction::

            {d: {"mrca_id": int,
                 "state": raw focal_state | None,
                 "posterior": P(anc_enc @ mrca_id) under the scheme, 0.0 if
                              the domain has no reconstruction}}
    """
    pairs = pair_details or []

    domain_meta: Dict[Any, Dict[str, Any]] = {}
    changed: Dict[str, List[Dict[str, Any]]] = {"top": [], "bottom": []}
    raw_der: Dict[str, Dict[Any, str]] = {"top": {}, "bottom": {}}
    raw_anc: Dict[str, Dict[Any, str]] = {"top": {}, "bottom": {}}
    der_enc: Dict[str, Dict[Any, str]] = {"top": {}, "bottom": {}}

    for pair in pairs:
        d = pair.get("pair_id")
        mrca_id = pair.get("node_id")
        if d is None or mrca_id is None:
            continue
        focal_raw = pair.get("focal_state")
        anc_enc = encode_aa(focal_raw, scheme)
        dist_mrca = node_dist(per_node_dist, mrca_id)

        if anc_enc is None:
            domain_meta[d] = {
                "mrca_id": int(mrca_id), "state": None, "posterior": 0.0,
            }
            continue

        domain_meta[d] = {
            "mrca_id": int(mrca_id),
            "state": str(focal_raw).strip().upper(),
            "posterior": encoded_distribution(dist_mrca, scheme).get(anc_enc, 0.0),
        }

        for side_key, side in (("top_tip_mode", "top"), ("bottom_tip_mode", "bottom")):
            raw_tip = pair.get(side_key)
            tip_enc = encode_aa(raw_tip, scheme)
            if tip_enc is None or tip_enc == anc_enc:
                continue  # missing or conserved side → not scored
            changed[side].append(
                {"d": d, "mrca_id": int(mrca_id),
                 "anc_enc": anc_enc, "der_enc": tip_enc}
            )
            raw_der[side][d] = str(raw_tip).strip().upper()
            raw_anc[side][d] = str(focal_raw).strip().upper()
            der_enc[side][d] = tip_enc

    out: Dict[str, Any] = {}
    for side in ("top", "bottom"):
        sc = score_domains_side(changed[side], node_index, per_node_dist, scheme)
        out[side] = {
            "domain_scores": sc["domain_scores"],
            "domain_der": raw_der[side],
            "domain_der_enc": der_enc[side],
            "domain_anc": raw_anc[side],
            "agree_num": sc["agree_num"],
            "agree_den": sc["agree_den"],
            "n_changed": sc["n_changed"],
            "convergence_type": _convergence_type(sc["agree_num"], sc["agree_den"]),
        }
    out["domain_meta"] = domain_meta
    return out
