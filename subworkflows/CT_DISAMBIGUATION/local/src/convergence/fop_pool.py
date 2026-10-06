#!/usr/bin/env python3
# fop_pool.py — Pool the per-hypothesis domain records of a position and scheme into one score per side.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/convergence/

"""
Multi-hypothesis (FOP) pooling of domain scores into one position score per phenotype side.

``pool_domains`` collapses ``M >= 1`` per-hypothesis ``compute_domain_scores``
records for one ``(Gene, Position, scheme)`` into one score per phenotype side.
It is **treeless**: every tree lookup already happened in
``path_scores.compute_domain_scores``; here we only average scalars.

Imported by: src/convergence/disambiguate_single.py, src/core/driver.py, disambiguation_perms_main.py and
reaggregate_perm_scores.py (`base_cycle`)
Inputs: lists of {"hyp", "sides"} records and an optional {(hyp, domain): pss} map (in memory)
Outputs: the pooled per-side aggregates described in `pool_domains`

Per side ``s`` over the universe of K fixed Voronoi domains
(``⋃_h keys(sides["domain_meta"]) ∪ keys(pss)``):

    M       = number of hypothesis records
    s̄_d    = ( Σ_h  domain_scores_h.get(d, 0.0) ) / M
    w̄_d    = ( Σ_h  pss(h, d, default 1.0) )      / M
    core_s  = clamp01( Σ_d w̄_d·s̄_d / Σ_d w̄_d )      (0 if Σ w̄_d == 0)

``M = 1`` degenerates exactly to the plain PSS-weighted mean over the K domains
(arithmetic mean when no PSS file). ``agree_num`` / ``agree_den`` and
``convergence_type`` are computed here, across all the pooled hypotheses, from the
modal encoded derived residue per domain (a domain that splits V/I/L between
hypotheses scores low under US and high under a scheme that co-encodes them; the
encoding is already applied in each ``domain_der_enc`` upstream, so ``pool_domains``
never sees the scheme).

See ``docs/scoring_v3_core.md`` for the scoring contract.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

from src.convergence.path_scores import _convergence_type
from src.convergence.support_fmt import fmt_support


def _tally_str(vals: Sequence[Optional[str]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for v in vals:
        if v:
            counts[str(v)] = counts.get(str(v), 0) + 1
    return counts


def _modal_str(vals: Sequence[Optional[str]]) -> Optional[str]:
    """Most frequent non-empty string. Settles a domain's derived / ancestral residue when several
    hypotheses reconstruct the same domain (they normally agree: tip residues are labeling-invariant).
    A tie goes to the smallest string, so the result does not depend on the order of the hypotheses;
    the choice is a convention, not a reading of the data (see :func:`_has_tie`)."""
    counts = _tally_str(vals)
    return min(counts, key=lambda s: (-counts[s], s)) if counts else None


def _has_tie(vals: Sequence[Optional[str]]) -> bool:
    """True when two or more strings share the highest count."""
    counts = _tally_str(vals)
    return bool(counts) and sum(1 for n in counts.values() if n == max(counts.values())) > 1


def _support_str(vals: Sequence[Optional[str]]) -> str:
    """'L:3,S:2'-style tally over a domain's per-hypothesis raw residues."""
    counts: Dict[str, int] = {}
    for v in vals:
        if not v:
            continue
        s = str(v)
        counts[s] = counts.get(s, 0) + 1
    return fmt_support(counts)


def base_cycle(tag: str) -> str:
    """'b_5~H3' -> 'b_5'; a plain 'b_5' passes through."""
    return tag.split("~", 1)[0]


# ── Domain id matching (int vs str keys) ──────────────────────────────────────
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
    """Treeless mean-of-means pooler over the fixed Voronoi domains.

    Args:
        hyp_records: ``list[{"hyp": str, "sides": <compute_domain_scores return>}]``,
            ``M >= 1``.
        pss_by_hyp_domain: ``{(hyp, domain): pss}``; ``None`` / missing key -> 1.0.

    Returns:
        ``{"top": <agg>, "bottom": <agg>, "n_hypotheses": M}`` where ``<agg>`` is
        ``{asr_path_score, domain_scores, domain_weights, domain_der,
        domain_anc, agree_num, agree_den, n_participating, convergence_type}``.
        ``domain_der`` / ``domain_anc`` carry only domains changed in >= 1
        hypothesis (modal residue); ``agree_den == n_participating``; ``agreement_tie`` flags a changed
        domain whose derived residue was tied; ``participating_hyps`` lists the
        hypotheses that changed >= 1 domain on the side.
    """
    hyps = [r for r in hyp_records if r.get("hyp")]
    M = len(hyps)

    # Domain universe: the domain_meta ids (all K domains) plus any ids only the PSS map has.
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
        # Sums are correctly rounded (math.fsum): the pooled score does not depend on the
        # order the hypotheses or the domains arrive in.
        for sd, d in rep.items():
            num = math.fsum(
                float(_dget(((r.get("sides") or {}).get(side) or {}).get("domain_scores"), d, 0.0) or 0.0)
                for r in hyps
            )
            s_bar[d] = (num / M) if M else 0.0
            wsum = math.fsum(
                1.0 if (w := _pss(r["hyp"], d)) is None else float(w) for r in hyps
            )
            w_bar[d] = (wsum / M) if M else 0.0

        denom = math.fsum(w_bar.values())
        core = (
            math.fsum(w_bar[d] * s_bar[d] for d in w_bar) / denom if denom > 0 else 0.0
        )
        core = max(0.0, min(1.0, core))

        # Agreement across the pooled hypotheses, over the domains changed in >= 1 of them.
        der_enc: Dict[str, List[str]] = {}
        der_raw: Dict[str, List[str]] = {}
        anc_raw: Dict[str, List[str]] = {}
        hyp_ids_by_domain: Dict[str, List[str]] = {}
        for r in hyps:
            srow = (r.get("sides") or {}).get(side) or {}
            for dd, e in (srow.get("domain_der_enc") or {}).items():
                der_enc.setdefault(str(dd), []).append(e)
                hyp_ids_by_domain.setdefault(str(dd), []).append(r["hyp"])
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
        # A tied domain was settled by the convention of _modal_str: agreement and convergence type may rest on it.
        agreement_tie = any(_has_tie(v) for v in der_enc.values())

        return {
            "asr_path_score": core,
            "domain_scores": dict(s_bar),
            "domain_weights": dict(w_bar),
            "domain_der": {rep.get(sd, sd): _modal_str(der_raw.get(sd, [])) for sd in changed},
            "domain_anc": {rep.get(sd, sd): _modal_str(anc_raw.get(sd, [])) for sd in changed},
            "domain_der_support": {rep.get(sd, sd): _support_str(der_raw.get(sd, [])) for sd in changed},
            "domain_anc_support": {rep.get(sd, sd): _support_str(anc_raw.get(sd, [])) for sd in changed},
            "agree_num": agree_num,
            "agree_den": agree_den,
            "agreement_tie": agreement_tie,
            "n_participating": agree_den,
            "convergence_type": _convergence_type(agree_num, agree_den),
            "participating_hyps": sorted({h for hs in hyp_ids_by_domain.values() for h in hs}),
        }

    return {"top": _agg("top"), "bottom": _agg("bottom"), "n_hypotheses": M}
