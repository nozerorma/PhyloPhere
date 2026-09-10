#!/usr/bin/env python3
"""FOP multi-hypothesis -> domain-pooled position score (scoring_v2 core v3).

``pool_domains`` collapses ``M >= 1`` per-hypothesis ``compute_domain_scores``
records for one ``(Gene, Position, scheme)`` into one score per phenotype side.
It is **treeless** — every tree lookup already happened in
``path_scores.compute_domain_scores``; here we only average scalars.

Per side ``s`` over the universe of K fixed Voronoi domains
(``⋃_h keys(sides["domain_meta"]) ∪ keys(pss)``):

    M       = number of hypothesis records
    s̄_d    = ( Σ_h  domain_scores_h.get(d, 0.0) ) / M
    w̄_d    = ( Σ_h  pss(h, d, default 1.0) )      / M
    core_s  = clamp01( Σ_d w̄_d·s̄_d / Σ_d w̄_d )      (0 if Σ w̄_d == 0)

``M = 1`` degenerates exactly to the plain PSS-weighted mean over the K domains
(arithmetic mean when no PSS file). The harvest-wide ``agree_num`` / ``agree_den``
/ ``convergence_type`` are recomputed here from the modal encoded derived residue
per domain across the harvest (a domain that splits V/I/L between hypotheses
scores low under US and high under a scheme that co-encodes them — the encoding
is already baked into each ``domain_der_enc`` upstream, ``pool_domains`` never
sees ``scheme``).

See ``docs/scoring_v3_core.md`` section 3 and Appendix B for the contract and the
hand-worked golden arithmetic.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from src.convergence.path_scores import _convergence_type


def _modal_str(vals: Sequence[Optional[str]]) -> Optional[str]:
    """Most frequent non-empty string (first-seen breaks ties). Settles a
    domain's derived / ancestral residue when several hypotheses reconstruct the
    same domain (they normally agree — tip residues are labeling-invariant)."""
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


def _num(x) -> Optional[float]:
    try:
        f = float(x)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def _wmean(xs: Sequence[float], ws: Sequence[Optional[float]]) -> Optional[float]:
    pairs = [(x, w) for x, w in zip(xs, ws) if x is not None and x == x]
    if not pairs:
        return None
    xs2 = [x for x, _ in pairs]
    ws2 = [w if (w is not None and w == w and w > 0) else 0.0 for _, w in pairs]
    if sum(ws2) <= 0:
        return sum(xs2) / len(xs2)
    return sum(x * w for x, w in zip(xs2, ws2)) / sum(ws2)


def base_cycle(tag: str) -> str:
    """'b_5~H3' -> 'b_5'; a plain 'b_5' passes through."""
    return tag.split("~", 1)[0]


# ── domain id matching (int vs str keys) ─────────────────────────────────────
# ``compute_domain_scores`` keys everything by the int ``pair_id``; a fixture
# round-tripped through JSON arrives str-keyed; the PSS map may be either. Match
# tolerantly, but keep the representative id (int when available) in the output.
def _cands(d):
    out = [d]
    if not isinstance(d, str):
        out.append(str(d))
    try:
        out.append(int(d))
    except (TypeError, ValueError):
        pass
    return out


def _dget(mapping: Optional[Dict], d, default=None):
    if not mapping:
        return default
    for k in _cands(d):
        if k in mapping:
            return mapping[k]
    return default


def pool_domains(
    hyp_records: List[Dict],
    pss_by_hyp_domain: Optional[Dict[Tuple[str, Any], float]] = None,
) -> Dict[str, Any]:
    """Treeless mean-of-means pooler over the K fixed Voronoi domains.

    Args:
        hyp_records: ``list[{"hyp": str, "sides": <compute_domain_scores return>}]``,
            ``M >= 1``.
        pss_by_hyp_domain: ``{(hyp, domain): pss}``; ``None`` / missing key -> 1.0.

    Returns:
        ``{"top": <agg>, "bottom": <agg>, "n_hypotheses": M}`` where ``<agg>`` is
        ``{asr_path_score, core, domain_scores, domain_weights, domain_der,
        domain_anc, agree_num, agree_den, n_participating, convergence_type}``;
        ``core == asr_path_score``. ``domain_der`` / ``domain_anc`` carry only
        domains changed in >= 1 hypothesis (modal residue); ``agree_den ==
        n_participating``.
    """
    hyps = [r for r in hyp_records if r.get("hyp")]
    M = len(hyps)

    # Canonical domain universe: domain_meta ids (all K) + any pss-only ids.
    rep: Dict[str, Any] = {}
    for r in hyps:
        for d in ((r.get("sides") or {}).get("domain_meta") or {}):
            rep.setdefault(str(d), d)
    if pss_by_hyp_domain:
        for (_h, d) in pss_by_hyp_domain:
            rep.setdefault(str(d), d)

    def _pss(h, d) -> Optional[float]:
        if not pss_by_hyp_domain:
            return None
        for k in _cands(d):
            if (h, k) in pss_by_hyp_domain:
                return pss_by_hyp_domain[(h, k)]
        return None

    def _agg(side: str) -> Dict[str, Any]:
        s_bar: Dict[Any, float] = {}
        w_bar: Dict[Any, float] = {}
        for sd, d in rep.items():
            num = 0.0
            for r in hyps:
                scores = ((r.get("sides") or {}).get(side) or {}).get("domain_scores")
                num += float(_dget(scores, d, 0.0) or 0.0)
            s_bar[d] = (num / M) if M else 0.0
            wsum = 0.0
            for r in hyps:
                w = _pss(r["hyp"], d)
                wsum += 1.0 if w is None else float(w)
            w_bar[d] = (wsum / M) if M else 0.0

        denom = sum(w_bar.values())
        core = (
            sum(w_bar[d] * s_bar[d] for d in w_bar) / denom if denom > 0 else 0.0
        )
        core = max(0.0, min(1.0, core))

        # Harvest-wide agreement over domains changed in >= 1 hypothesis.
        der_enc: Dict[str, List[str]] = {}
        der_raw: Dict[str, List[str]] = {}
        anc_raw: Dict[str, List[str]] = {}
        for r in hyps:
            srow = (r.get("sides") or {}).get(side) or {}
            for dd, e in (srow.get("domain_der_enc") or {}).items():
                der_enc.setdefault(str(dd), []).append(e)
            for dd, e in (srow.get("domain_der") or {}).items():
                der_raw.setdefault(str(dd), []).append(e)
            for dd, e in (srow.get("domain_anc") or {}).items():
                anc_raw.setdefault(str(dd), []).append(e)

        changed = list(der_enc.keys())
        agree_den = len(changed)
        modal_enc = {sd: _modal_str(v) for sd, v in der_enc.items()}
        grp: Dict[Optional[str], int] = {}
        for sd in changed:
            grp[modal_enc[sd]] = grp.get(modal_enc[sd], 0) + 1
        agree_num = max(grp.values()) if grp else 0

        return {
            "asr_path_score": core,
            "core": core,
            "domain_scores": dict(s_bar),
            "domain_weights": dict(w_bar),
            "domain_der": {rep.get(sd, sd): _modal_str(der_raw.get(sd, [])) for sd in changed},
            "domain_anc": {rep.get(sd, sd): _modal_str(anc_raw.get(sd, [])) for sd in changed},
            "agree_num": agree_num,
            "agree_den": agree_den,
            "n_participating": agree_den,
            "convergence_type": _convergence_type(agree_num, agree_den),
        }

    return {"top": _agg("top"), "bottom": _agg("bottom"), "n_hypotheses": M}
