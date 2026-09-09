#!/usr/bin/env python3
"""FOP multi-hypothesis -> domain-pooled position score (permulation-null side).

Python twin of ``subworkflows/SCORING/local/src/fop_pool.R``. The observed side
pools the FOP hypothesis harvest inside scoring_compute.R section 2b; the
permulation null must hold the SAME statistic or the FCS p.perm comparison is
miscalibrated, so this module reproduces that algebra for the null's per-cycle
replay.

Given, for one ``(base_cycle, Gene, Position, scheme)``, the per-hypothesis
records H1..Hn (each carrying ``asr_path_score`` plus ``pair_scores`` =
``{domain_id -> s(p, site)}``) and the per-(hypothesis, domain) PSS weight from
``resample_fop_pairs.tsv``:

* within each Voronoi domain d, pool ``s(p, site)`` over the distinct pairs that
  domain contributed across the harvest -> ``c_d`` (PSS-weighted mean; equal
  weights when PSS missing);

TWO-JOB PSS WEIGHTING (refined 2026-09-03; mirrors fop_pool.R). A single PSS per
hypothesis collapses to equal-weight in practice (every FOP hypothesis carries
the harvest's globally weakest domain), so the split does two jobs:

* Job A -- the per-domain ``c_d`` pool: weight each pair by its OWN PSS in domain
  d, ``pss_by_hyp_domain[(hyp, d)]``, not by the containing hypothesis's summary.
* Job B -- the axis pools (independence / mrca_diversity / derived_agreement /
  conservation_gate): per-hypothesis weight = MEAN pair PSS across its domains.

Residual R/py difference (``c_d`` only): the null record carries no MRCA-node
identity, so ``c_d`` here has no node dedup — it wmean's the per-hypothesis
domain scores directly. Same weight *scheme*, not bit-identical. (The POINT 3
derived_agreement no longer has this gap — R and py both build a PSS-weighted
residue distribution per (domain, side) with no dedup; see below.)

DIRECTIONAL core (fixed 2026-09-04; mirrors path_scores.py's own
``core_top``/``core_bottom`` split and fop_pool.R::pool_group). ``pair_scores``
is only the per-pair TOP/BOTTOM *average* — pooling it directly let a domain
whose corroborating pair changed on the top phenotype side and a different
domain whose pair changed on the bottom side count as ">=2 domains carrying a
clean change", when they are two unrelated, opposite-direction events (exactly
what ``_p_at_least_2`` requiring >=2 *same-side* observations exists to rule
out). When hypothesis records carry ``pair_top_scores``/``pair_bottom_scores``
(``{domain_id -> s}``, the un-averaged per-side scores), each side is pooled
SEPARATELY over its own distinct domains and combined the way
``compute_asr_path_score`` does:
* ``core_top    = P(>=2 domains carry a clean TOP-side change)``
* ``core_bottom = P(>=2 domains carry a clean BOTTOM-side change)``
* ``core = 1 - (1 - core_top) * (1 - core_bottom)``
Older records without the directional fields fall back to pooling the
side-averaged ``pair_scores`` in one undifferentiated ``core`` (the prior,
known-approximate behaviour) so a stale cache still scores instead of erroring.
* ``independence`` / ``mrca_diversity`` / ``conservation_gate`` are pooled
  across hypotheses by the same PSS weights (carried on the null's PositionAxes
  record); ``derived_agreement`` is rebuilt harvest-wide from the PSS-weighted
  per-(domain, side) residue distribution (POINT 3, ``_domain_side_dists`` +
  ``_da_frac_from_dists``);
* recombine with the path_scores.py algebra, verbatim from fop_pool.R (T1: the
  mrca_diversity and conservation_gate multipliers are gone; ``div`` / ``cg``
  are still pooled and reported, latent)::

      replication = independence_pooled * core_pooled
      asr_pooled  = replication * derived_agreement_pooled

Single hypothesis (or non-FOP) -> the lone ``asr_path_score`` passes through.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

try:  # src on path at runtime (gene_wrapper import); absent under isolated unit test
    from src.biochem.grouping import get_grouping_scheme as _get_grouping_scheme
except Exception:  # pragma: no cover
    _get_grouping_scheme = None

try:  # scoring_v2 T3c SC2b — the real per-side pairwise core for the FOP null
    from src.convergence.path_scores import aggregate_core_side as _aggregate_core_side
except Exception:  # pragma: no cover
    _aggregate_core_side = None


def _encode_aa(aa: Optional[str], scheme: Optional[str]) -> Optional[str]:
    """Encode a raw residue in ``scheme`` (mirror of path_scores.encode_aa).

    US / unknown scheme -> raw upper-cased residue; GS scheme -> group label,
    falling back to the raw residue for an unknown residue.
    """
    if not aa:
        return None
    raw = str(aa).strip().upper()
    if not raw:
        return None
    if not scheme or scheme.upper() == "US" or _get_grouping_scheme is None:
        return raw
    return _get_grouping_scheme(raw, scheme) or raw


def _rebuild_derived_agreement(
    changed: List[Tuple[str, str]], scheme: Optional[str]
) -> float:
    """Legacy MODAL harvest-wide derived_agreement (superseded 2026-09-07 by the
    PSS-weighted distribution form, ``_domain_side_dists`` + ``_da_frac_from_dists``).

    Kept for reference / older callers. ``changed`` is a list of ``(side, raw_aa)``
    over the null's changed pairs; plain plurality-count concentration per side.
    """
    concentrations: List[float] = []
    for want_side in ("top", "bot"):
        rs = [aa for sd, aa in changed if sd == want_side and aa]
        if len(rs) < 2:
            continue
        enc = [_encode_aa(r, scheme) for r in rs]
        counts: Dict[Optional[str], int] = {}
        for e in enc:
            counts[e] = counts.get(e, 0) + 1
        concentrations.append(max(counts.values()) / len(enc))
    if not concentrations:
        return 1.0
    return sum(concentrations) / len(concentrations)


def _domain_side_dists(
    hyps: List[Dict],
    pss_by_hyp_domain: Optional[Dict[Tuple[str, int], float]],
) -> Dict[Tuple[int, str], Dict[str, float]]:
    """PSS-weighted residue DISTRIBUTION per (Voronoi domain, side) — Python twin
    of fop_pool.R::.collect_changed_pair_dists, for the permulation null.

    For (domain d, side s): over the hypothesis records that recorded a change at
    that cell (``pair_derived_top`` / ``pair_derived_bot`` = ``{domain -> raw}``,
    only the changed side populated), accumulate weight
    ``pss_by_hyp_domain[(hyp, d)]`` (1.0 when missing/NaN), then normalise to
    sum 1. No ``(domain, side, raw)`` dedup — every record contributes its weight.
    """
    acc: Dict[Tuple[int, str], Dict[str, float]] = {}
    for r in hyps:
        h = r.get("hyp")
        for side, fld in (("top", "pair_derived_top"), ("bot", "pair_derived_bot")):
            for d, raw in (r.get(fld) or {}).items():
                if not raw:
                    continue
                raw_u = str(raw).strip().upper()
                if not raw_u or raw_u in ("NA", "NAN", "NONE"):
                    continue
                w = pss_by_hyp_domain.get((h, d)) if pss_by_hyp_domain else None
                if w is None or w != w:  # None or NaN -> equal weight
                    w = 1.0
                bucket = acc.setdefault((d, side), {})
                bucket[raw_u] = bucket.get(raw_u, 0.0) + float(w)
    out: Dict[Tuple[int, str], Dict[str, float]] = {}
    for k, resid in acc.items():
        tot = sum(resid.values())
        if tot > 0:
            out[k] = {res: v / tot for res, v in resid.items()}
    return out


def _da_frac_from_dists(
    dists: Dict[Tuple[int, str], Dict[str, float]], scheme: Optional[str]
) -> float:
    """Harvest-wide, per-scheme derived_agreement from PSS-weighted residue
    distributions — Python twin of fop_pool.R::rebuild_derived_agreement_frac.

        concentration_s = max_g( sum over domains d on side s of p_d^enc(g) ) / |D_s|
        da              = mean over sides with |D_s| >= 2 (1.0 when none qualify)

    A domain unanimous across the harvest is a point mass, so this equals the
    old modal plurality exactly on runs without split domains.
    """
    if not dists:
        return 1.0
    concentrations: List[float] = []
    for want_side in ("top", "bot"):
        keys = [k for k in dists if k[1] == want_side]
        if len(keys) < 2:
            continue
        gtot: Dict[Optional[str], float] = {}
        for k in keys:
            for res, p in dists[k].items():
                g = _encode_aa(res, scheme)
                gtot[g] = gtot.get(g, 0.0) + p
        concentrations.append(max(gtot.values()) / len(keys))
    if not concentrations:
        return 1.0
    return sum(concentrations) / len(concentrations)


def p_at_least_2(ps: Sequence[float]) -> float:
    """Exact P(>=2 successes) over independent Bernoullis (path_scores.py algebra)."""
    vals = [p for p in ps if p is not None and p == p]  # drop None / NaN
    if len(vals) < 2:
        return 0.0
    p0 = 1.0
    for p in vals:
        p0 *= (1.0 - p)
    p1 = 0.0
    for i, pi in enumerate(vals):
        term = pi
        for j, pj in enumerate(vals):
            if j != i:
                term *= (1.0 - pj)
        p1 += term
    return max(0.0, min(1.0, 1.0 - p0 - p1))


def _wmean(xs: Sequence[float], ws: Sequence[Optional[float]]) -> Optional[float]:
    pairs = [(x, w) for x, w in zip(xs, ws) if x is not None and x == x]
    if not pairs:
        return None
    xs2 = [x for x, _ in pairs]
    ws2 = [w if (w is not None and w == w and w > 0) else 0.0 for _, w in pairs]
    if sum(ws2) <= 0:
        return sum(xs2) / len(xs2)
    return sum(x * w for x, w in zip(xs2, ws2)) / sum(ws2)


def pool_hypotheses(
    hyp_records: List[Dict],
    pss_by_hyp_domain: Optional[Dict[Tuple[str, int], float]] = None,
    scheme: Optional[str] = None,
) -> Dict[str, float]:
    """Domain-pool one position's FOP hypothesis records.

    Args:
        hyp_records: list of dicts, one per hypothesis, each with keys
            ``hyp`` (e.g. "H3"), ``asr_path_score`` (float), ``pair_scores``
            ({domain_id:int -> s(p, site):float}, the per-pair side-average),
            and — when available — ``pair_top_scores`` / ``pair_bottom_scores``
            (same shape, the un-averaged per-side scores `core` should
            actually be built from; see the DIRECTIONAL core note above).
        pss_by_hyp_domain: {(hyp, domain) -> pss_score}; None -> equal weights.

    Returns:
        {"asr_path_score", "core", "independence", "mrca_diversity",
         "derived_agreement", "conservation_gate", "n_hypotheses"}.
    """
    hyps = [r for r in hyp_records if r.get("hyp")]
    if len(hyps) <= 1:
        r = hyp_records[0] if hyp_records else {}
        return {
            "asr_path_score": float(r.get("asr_path_score", 0.0) or 0.0),
            "core": float(r.get("core", 0.0) or 0.0),
            "independence": r.get("independence"),
            "mrca_diversity": r.get("mrca_diversity"),
            "derived_agreement": r.get("derived_agreement"),
            "conservation_gate": r.get("conservation_gate"),
            "n_hypotheses": len(hyps),
        }

    # Job B weight: per-hypothesis MEAN pair PSS across its domains (overall
    # credibility of the hypothesis; used for the axis pools). Mirrors fop_pool.R.
    def hyp_weight(h: str) -> Optional[float]:
        if not pss_by_hyp_domain:
            return None
        ws = [v for (hh, _d), v in pss_by_hyp_domain.items()
              if hh == h and v is not None and v == v]
        return sum(ws) / len(ws) if ws else None

    hw = {r["hyp"]: hyp_weight(r["hyp"]) for r in hyps}
    row_w = [hw.get(r["hyp"]) for r in hyps]

    # Job A: c_d = PSS-weighted mean of s(p, site) over the pairs domain d
    # contributed, each weighted by ITS OWN PSS in domain d (not the hypothesis
    # summary weight). No node-dedup here — the null record has no pair identity.
    domains = sorted({
        d for r in hyps
        for d in list((r.get("pair_scores") or {}).keys())
        + list((r.get("pair_top_scores") or {}).keys())
        + list((r.get("pair_bottom_scores") or {}).keys())
    })

    def _domain_pool(field: str) -> List[float]:
        out: List[float] = []
        for d in domains:
            xs, ws = [], []
            for r in hyps:
                s = (r.get(field) or {}).get(d)
                if s is None:
                    continue
                xs.append(float(s))
                ws.append(pss_by_hyp_domain.get((r["hyp"], d)) if pss_by_hyp_domain else None)
            m = _wmean(xs, ws)
            if m is not None:
                out.append(m)
        return out

    # DIRECTIONAL core (mirrors path_scores.py::compute_asr_path_score and
    # fop_pool.R::pool_group). `pair_scores` is only the per-pair TOP/BOTTOM
    # average — pooling it directly lets a domain that changed on top and a
    # different domain that changed on bottom count as "2 corroborating
    # domains" when they are two unrelated, opposite-direction events. Pool
    # each side's per-domain scores SEPARATELY and combine the way
    # core_top/core_bottom do, using `pair_top_scores`/`pair_bottom_scores`
    # when the hop record carries them.
    has_side_scores = any(
        (r.get("pair_top_scores") or r.get("pair_bottom_scores")) for r in hyps
    )
    if has_side_scores:
        core_top = p_at_least_2(_domain_pool("pair_top_scores"))
        core_bottom = p_at_least_2(_domain_pool("pair_bottom_scores"))
        core_pooled = 1.0 - (1.0 - core_top) * (1.0 - core_bottom)
    else:
        # Fallback for older null records without the directional fields:
        # collapses both directions into one pool (the known approximation —
        # can over-credit opposite-direction domains as if they corroborated).
        # Kept only so stale/older cached inputs still score instead of erroring.
        core_pooled = p_at_least_2(_domain_pool("pair_scores"))

    def _pool(field: str, default: float) -> float:
        v = _wmean([_num(r.get(field)) for r in hyps], row_w)
        return v if v is not None else default

    indep = _pool("independence", 1.0)
    div = _pool("mrca_diversity", 0.0)

    # derived_agreement (POINT 3): recompute HARVEST-WIDE under the null's active
    # scheme when the records carry raw per-side derived residues, mirroring
    # fop_pool.R::pool_group. Per (domain, side) build a PSS-weighted residue
    # distribution over the hypotheses and take a weighted plurality
    # concentration (_domain_side_dists + _da_frac_from_dists). Two null
    # hypotheses can each be internally unanimous yet land on different residues
    # in the same domain; a distribution keeps that split visible. A domain
    # unanimous across the harvest is a point mass -> equals the old modal value
    # exactly. Fall back to the per-hypothesis wmean when the raw residues are
    # absent (older null records / non-FOP).
    has_raw = any(
        ("pair_derived_top" in r or "pair_derived_bot" in r) for r in hyps
    )
    if has_raw:
        da = _da_frac_from_dists(
            _domain_side_dists(hyps, pss_by_hyp_domain), scheme)
    else:
        da = _pool("derived_agreement", 1.0)

    # conservation_gate (POINT 2): if the per-hypothesis records carry conserved
    # pair identity (conserved_pair_scores / conserved_pair_nodes, attached to the
    # null PositionAxes from compute_asr_path_score), rebuild the gate from the
    # DISTINCT conserved pairs shared across hypotheses (dedup by MRCA node) so a
    # pair conserved under several hypotheses is counted ONCE, symmetric with the
    # node-deduped changed-pair c_d pool. Otherwise fall back to pooling the
    # already-transformed per-hypothesis gates (residual: mean(0.5+0.5 cons_h)
    # != 0.5+0.5 mean(cons_p) when hypotheses differ in conserved-pair count).
    cons_by_node: Dict[Any, float] = {}
    for r in hyps:
        cs = r.get("conserved_pair_scores") or {}
        cn = r.get("conserved_pair_nodes") or {}
        for pid, sc in cs.items():
            node = cn.get(pid, ("_nokey_", pid))
            v = _num(sc)
            if v is None:
                continue
            # PSS weight for a conserved pair is not resolvable from the null
            # record (no domain mapping) -> equal weight. Keep the max score per
            # node if two hypotheses disagree (same pair, same posteriors -> same).
            if node not in cons_by_node or v > cons_by_node[node]:
                cons_by_node[node] = v
    if cons_by_node:
        cg = 0.5 + 0.5 * (sum(cons_by_node.values()) / len(cons_by_node))
    else:
        # "columns present" test mirrors fop_pool.R's have_cons_cols: the KEY is
        # on the record (even as an empty dict) -> the conserved block exists for
        # this run. Empty -> novel position -> gate 1.0. Key absent entirely
        # (older null records) -> fall back to pooling per-hypothesis gates.
        has_cons_cols = any(
            ("conserved_pair_scores" in r or "conserved_pair_nodes" in r)
            for r in hyps
        )
        cg = 1.0 if has_cons_cols else _pool("conservation_gate", 1.0)

    # T1 collapse (mirrors path_scores.compute_asr_path_score): drop the
    # mrca_diversity (diversity_mult) and conservation_gate multipliers.
    # ``div`` and ``cg`` / ``cons_by_node`` stay computed above — latent, kept
    # for a later tier and still reported below — but no longer enter the score.
    asr_pooled = max(0.0, min(1.0, indep * core_pooled * da))

    return {
        "asr_path_score": asr_pooled,
        "core": core_pooled,
        "independence": indep,
        "mrca_diversity": div,
        "derived_agreement": da,
        "conservation_gate": cg,
        "n_hypotheses": len(hyps),
    }


def _modal_str(vals: Sequence[Optional[str]]) -> Optional[str]:
    """Most frequent non-empty string (first-seen breaks ties). Used to settle a
    node's derived/ancestral residue when several hypotheses reconstruct the same
    MRCA node (they normally agree — tip residues are labeling-invariant)."""
    counts: Dict[str, int] = {}
    order: List[str] = []
    for v in vals:
        if not v:
            continue
        s = str(v)
        if s not in counts:
            order.append(s)
        counts[s] = counts.get(s, 0) + 1
    if not order:
        return None
    return max(order, key=lambda s: (counts[s], -order.index(s)))


def pool_hypotheses_pairwise(
    hyp_records: List[Dict],
    node_index: Dict[int, Any],
    per_node_dist: Dict[Any, Dict[str, float]],
    pss_by_hyp_domain: Optional[Dict[Tuple[str, int], float]] = None,
    scheme: Optional[str] = None,
) -> Dict[str, Any]:
    """Per-side pairwise ``core_s`` over the FOP hypothesis harvest (T3c SC2b).

    The real pairwise aggregation of ``path_scores.aggregate_core_side``, run over
    the **node-deduped union** of the harvest's participating pairs — not the
    ``p_at_least_2`` per-domain approximation of :func:`pool_hypotheses`. Needs
    the gene's tree (``node_index``) and the focal site's posteriors
    (``per_node_dist``), both available in ``_perms_worker``.

    Args:
        hyp_records: one dict per hypothesis for a single
            ``(base_cycle, Position, scheme)``, each with ``hyp`` and
            ``sides`` = ``{"top": <side_dict>, "bottom": <side_dict>}`` where a
            ``<side_dict>`` is a ``compute_asr_path_score(native_side_split=True)``
            side row (carries ``pair_scores`` = {pid: s_c^s}, ``pair_mrca``,
            ``pair_der_enc``, ``pair_anc_enc``, ``conserved_pair_nodes``).
        pss_by_hyp_domain: ``{(hyp, domain) -> pss}``; domain == pid by FOP
            construction. ``None`` / missing → equal weight.

    Returns:
        ``{"top": <agg>, "bottom": <agg>, "n_hypotheses": n}`` — each ``<agg>``
        the full :func:`aggregate_core_side` output for that side (its
        ``asr_path_score`` is ``core_s``; no ``1-(1-t)(1-b)`` recombination,
        no ``independence``/``derived_agreement`` multiplier).
    """
    hyps = [r for r in hyp_records if r.get("hyp")]
    if _aggregate_core_side is None:  # pragma: no cover — import guard
        raise RuntimeError("path_scores.aggregate_core_side unavailable")

    def _side(side: str) -> Dict[str, Any]:
        # node -> [{der_enc, anc_enc, iso, w}]  over every hypothesis that put a
        # participating pair on that MRCA node on this side.
        by_node: Dict[Any, List[Dict[str, Any]]] = {}
        cons_nodes: set = set()
        for r in hyps:
            sd = (r.get("sides") or {}).get(side) or {}
            iso_map = sd.get("pair_scores") or {}
            mrca_map = sd.get("pair_mrca") or {}
            der_map = sd.get("pair_der_enc") or {}
            anc_map = sd.get("pair_anc_enc") or {}
            for pid, iso in iso_map.items():
                node = mrca_map.get(pid)
                if node is None:
                    continue
                w = None
                if pss_by_hyp_domain:
                    w = pss_by_hyp_domain.get((r["hyp"], pid))
                    if w is None:
                        try:
                            w = pss_by_hyp_domain.get((r["hyp"], int(pid)))
                        except (TypeError, ValueError):
                            w = None
                by_node.setdefault(node, []).append({
                    "der_enc": der_map.get(pid), "anc_enc": anc_map.get(pid),
                    "iso": _num(iso), "w": w,
                })
            for cn in (sd.get("conserved_pair_nodes") or {}).values():
                if cn is not None:
                    cons_nodes.add(cn)

        participants: List[Dict[str, Any]] = []
        iso_override: Dict[Any, float] = {}
        for node in sorted(by_node, key=str):
            rows = by_node[node]
            iso = _wmean([x["iso"] for x in rows], [x["w"] for x in rows])
            der = _modal_str([x["der_enc"] for x in rows])
            anc = _modal_str([x["anc_enc"] for x in rows])
            if der is None or anc is None:
                continue  # cannot score a pair without its residues
            participants.append(
                {"pid": node, "mrca_id": node, "der_enc": der, "anc_enc": anc}
            )
            iso_override[node] = 0.0 if iso is None else float(iso)

        # A conserved node that also hosts a participating pair on this side is
        # already counted as a participant — don't double it in the denominator.
        n_conserved = len(cons_nodes - set(iso_override))
        return _aggregate_core_side(
            participants, n_conserved, node_index, per_node_dist, scheme,
            iso_override=iso_override,
        )

    return {
        "top": _side("top"),
        "bottom": _side("bottom"),
        "n_hypotheses": len(hyps),
    }


def _num(x) -> Optional[float]:
    try:
        f = float(x)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def base_cycle(tag: str) -> str:
    """'b_5~H3' -> 'b_5'; a plain 'b_5' passes through."""
    return tag.split("~", 1)[0]
