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
    w_d,h   = (1 - δ) + 2δ·(log(1+pss) - min) / (max - min) with δ = 1 / (2K)
    w̄_d    = ( Σ_h  w_d,h )                       / M
    score_d = ( Σ_d w̄_d·s̄_d / Σ_d w̄_d )            (0 if Σ w̄_d == 0)
    core_s  = clamp01( score_d * (M / (M + 2.0)) ) (shrinkage when the PSS map holds more than one hypothesis)

Without a PSS map, w_d,h = 1 and the shrinkage factor is 1, so ``M = 1`` is the arithmetic mean over the
K domains. With a PSS map the weights are the bounded ones above; the score is also shrunk by M / (M + 2)
when the map holds more than one hypothesis, and not shrunk when it holds a single one. ``agree_num`` / ``agree_den`` and
``convergence_type`` are computed here, across all the pooled hypotheses, from the
modal encoded derived residue per domain (a domain that splits V/I/L between
hypotheses scores low under US and high under a scheme that co-encodes them; the
encoding is already applied in each ``domain_der_enc`` upstream, so ``pool_domains``
never sees the scheme).

See ``docs/scoring_v3_core.md`` for the scoring contract.
"""

from __future__ import annotations

import collections
import os

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


_VARIANT = os.environ.get("PHYLOPHERE_POOL_VARIANT", "base")   # experiment switch (not committed)


def effective_number(supports) -> float:
    """Effective number of independent hypotheses among the supports of the hypotheses that detect a position on one side.

    A support is the set of (domain, MRCA node of the pair) that carry the change in one hypothesis, so two
    hypotheses that confirm the position with the same pairs are one piece of evidence, not two. The number is
    M^2 / sum_ij Z_ij (the Hill number of order 2 with a similarity matrix, Z_ij the Jaccard similarity of the
    supports): M when no two supports share a pair, 1 when all are identical, 0 without a support.
    """
    sups = [s for s in supports if s]
    m = len(sups)
    if m == 0:
        return 0.0
    groups = list(collections.Counter(sups).items())   # identical supports collapse into one group
    total = 0.0
    for a, ca in groups:
        for b, cb in groups:
            total += ca * cb * (len(a & b) / len(a | b))
    return (m * m) / total


def pool_domains(
    hyp_records: List[Dict],
    pss_by_hyp_domain: Optional[Dict[Tuple[str, Any], float]] = None,
    origin_fn=None,
) -> Dict[str, Any]:
    """Treeless mean-of-means pooler over the fixed Voronoi domains.

    Args:
        hyp_records: ``list[{"hyp": str, "sides": <compute_domain_scores return>}]``,
            ``M >= 1``.
        origin_fn: optional ``origin_fn(residue, species) -> int``, the number of independent origins of a residue among
            carrier species on the species tree (the pooler itself stays treeless); it gives ``n_orig`` of a side.
        pss_by_hyp_domain: ``{(hyp, domain): pss}`` of the design (every hypothesis it holds, not only the ones
            that detect the position); ``None`` / empty / missing key -> weight 1.0. The number of distinct
            hypotheses of this map is the design size H: the support shrinkage applies only when H > 1.

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

    # ── Component A: K-bounded PSS domain weights (delta = 1 / (2*K)) ──────────
    # Within each hypothesis h, compress raw PSS with log(1 + PSS) and linearly
    # rescale to [1 - delta, 1 + delta]. Guarantees that m+1 domains strictly
    # outweigh m domains for any K >= 2.
    K = len(rep)
    delta = (1.0 / (2.0 * K)) if K >= 1 else 0.0
    w_by_hyp_domain: Dict[Tuple[str, Any], float] = {}
    has_pss = bool(pss_by_hyp_domain)
    if has_pss and M >= 1:
        for r in hyps:
            h = r["hyp"]
            log_vals: Dict[Any, float] = {}
            for sd, d in rep.items():
                raw_w = _pss(h, d)
                p = 1.0 if raw_w is None else float(raw_w)
                log_vals[d] = math.log(1.0 + max(0.0, p))
            min_l = min(log_vals.values()) if log_vals else 0.0
            max_l = max(log_vals.values()) if log_vals else 0.0
            span = max_l - min_l
            for d, lv in log_vals.items():
                if span < 1e-12:
                    w_by_hyp_domain[(h, d)] = 1.0
                else:
                    w_by_hyp_domain[(h, d)] = (1.0 - delta) + 2.0 * delta * ((lv - min_l) / span)
    else:
        for r in hyps:
            h = r["hyp"]
            for sd, d in rep.items():
                w_by_hyp_domain[(h, d)] = 1.0

    if "nopss" in _VARIANT:
        for key in list(w_by_hyp_domain):
            w_by_hyp_domain[key] = 1.0
    w_bar: Dict[Any, float] = {}
    for sd, d in rep.items():
        wsum = math.fsum(w_by_hyp_domain.get((r["hyp"], d), 1.0) for r in hyps)
        w_bar[d] = (wsum / M) if M else 0.0
    denom = math.fsum(w_bar.values())

    # ── Component B: FOP hypothesis shrinkage f_FOP(M) = M / (M + 2) ───────────
    # Penalizes positions detected in only 1 or 2 isolated FOP hypotheses while
    # preserving robust positions detected across many hypotheses. It is a property of the design,
    # read from the number of hypotheses H of the PSS map (the same map the observed labeling and every
    # null cycle carry), so a single-hypothesis design is never shrunk and the observed score and the
    # null scores sit on one scale.
    k_0 = 2.0
    n_hyp_design = len({h for (h, _d) in pss_by_hyp_domain}) if has_pss else 0
    f_fop = (M / (M + k_0)) if (n_hyp_design > 1 and M >= 1) else 1.0

    def _redundancy_weights(side: str):
        """Weight of each hypothesis, 1 / sum_h' Z_hh' (Z the Jaccard similarity of the supports, an empty support
        being identical to another empty one only), normalized to 1: hypotheses with the same support share one unit."""
        sups = []
        for r in hyps:
            rs = r.get("sides") or {}
            ch = (rs.get(side) or {}).get("domain_der_enc") or {}
            meta_h = rs.get("domain_meta") or {}
            sups.append(frozenset((str(dd), (meta_h.get(dd) or {}).get("mrca_id")) for dd in ch))
        groups = collections.Counter(sups)
        def z(a, b):
            if not a and not b:
                return 1.0
            if not a or not b:
                return 0.0
            return len(a & b) / len(a | b)
        wgt = {g: 1.0 / sum(cb * z(g, b) for b, cb in groups.items()) for g in groups}
        v = [wgt[g] for g in sups]
        tot = math.fsum(v)
        return [x / tot for x in v]

    def _agg(side: str) -> Dict[str, Any]:
        s_bar: Dict[Any, float] = {}
        vw = _redundancy_weights(side) if ("redund" in _VARIANT and M) else None
        w_bar_use, denom_use = w_bar, denom
        if vw:
            w_bar_use = {d: math.fsum(v * w_by_hyp_domain.get((r["hyp"], d), 1.0) for v, r in zip(vw, hyps)) for sd, d in rep.items()}
            denom_use = math.fsum(w_bar_use.values())
        # Sums are correctly rounded (math.fsum): the pooled score does not depend on the
        # order the hypotheses or the domains arrive in.
        for sd, d in rep.items():
            vals = [float(_dget(((r.get("sides") or {}).get(side) or {}).get("domain_scores"), d, 0.0) or 0.0) for r in hyps]
            if vw:
                s_bar[d] = math.fsum(v * x for v, x in zip(vw, vals))
            else:
                s_bar[d] = (math.fsum(vals) / M) if M else 0.0

        score_dom = (
            math.fsum(w_bar_use[d] * s_bar[d] for d in w_bar_use) / denom_use if denom_use > 0 else 0.0
        )
        if _VARIANT.startswith("hier") and M:
            # two levels: weighted mean over the domains inside each hypothesis (weights and scores of the same
            # pairs), then a mean over the hypotheses with redundancy-discounted weights
            vh = _redundancy_weights(side)
            per_h = []
            for r in hyps:
                dsc = ((r.get("sides") or {}).get(side) or {}).get("domain_scores")
                num = math.fsum(w_by_hyp_domain.get((r["hyp"], d), 1.0) * float(_dget(dsc, d, 0.0) or 0.0) for sd, d in rep.items())
                den = math.fsum(w_by_hyp_domain.get((r["hyp"], d), 1.0) for sd, d in rep.items())
                per_h.append(num / den if den > 0 else 0.0)
            score_dom = math.fsum(v * x for v, x in zip(vh, per_h))
        if _VARIANT.startswith("union") and M:
            # no hypotheses: the distinct pairs of species of two domains that carry the same derived residue (each
            # counted once whatever the number of hypotheses that use it) give the independence of the two domains;
            # a domain scores the noisy-OR of its partners and the position the plain mean over the domains
            comb = {}
            pss_max = max(pss_by_hyp_domain.values()) if has_pss else None
            for r in hyps:
                srow = (r.get("sides") or {}).get(side) or {}
                spd = srow.get("domain_der_species") or {}
                for a, b, _lca, ct in (srow.get("pair_lca") or []):
                    sa, sb = spd.get(a), spd.get(b)
                    if sa is not None and sb is not None:
                        if "pssin" in _VARIANT and pss_max:
                            # the PSS of the two pairs scales the contribution of the pair of domains, before the noisy-OR
                            wa, wb = [math.log1p(_pss(r["hyp"], d) or pss_max) / math.log1p(pss_max) for d in (a, b)]
                            ct = ct * math.sqrt(wa * wb)
                        comb[(a, b, sa, sb)] = ct
            by_pair = collections.defaultdict(list)
            for (a, b, _sa, _sb), ct in comb.items():
                by_pair[(a, b)].append(ct)
            cdom = collections.defaultdict(list)
            for (a, b), cts in by_pair.items():
                c = max(cts) if _VARIANT.endswith("max") else math.fsum(cts) / len(cts)
                cdom[a].append(c)
                cdom[b].append(c)
            s_dom = [(1.0 - math.prod(1.0 - c for c in cdom[d])) if cdom.get(d) else 0.0 for d in rep.values()]
            score_dom = math.fsum(s_dom) / len(s_dom) if s_dom else 0.0
            if _VARIANT.endswith("_pss") and denom > 0:
                # the mean of the domains weighted by the bounded PSS of each domain
                score_dom = math.fsum(w_bar[d] * x for d, x in zip(rep.values(), s_dom)) / denom
        core = max(0.0, min(1.0, score_dom * f_fop))

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

        supports = []
        for r in hyps:
            rs = r.get("sides") or {}
            changed_h = (rs.get(side) or {}).get("domain_der_enc") or {}
            if changed_h:
                meta_h = rs.get("domain_meta") or {}
                supports.append(frozenset((str(dd), (meta_h.get(dd) or {}).get("mrca_id")) for dd in changed_h))

        # Species that carry the derived residue in a domain that converges (a positive domain score) in at least
        # one hypothesis: the independent observations behind the position, whatever the hypotheses that used them.
        species = set()
        for r in hyps:
            srow = (r.get("sides") or {}).get(side) or {}
            ds = srow.get("domain_scores") or {}
            for dd, sp in (srow.get("domain_der_species") or {}).items():
                if float(ds.get(dd, 0.0) or 0.0) > 0.0:
                    species.add(sp)

        n_orig = None
        if origin_fn is not None:
            carriers_by_residue = collections.defaultdict(set)
            for r in hyps:
                srow = (r.get("sides") or {}).get(side) or {}
                ds = srow.get("domain_scores") or {}
                enc = srow.get("domain_der_enc") or {}
                for dd, sp in (srow.get("domain_der_species") or {}).items():
                    if float(ds.get(dd, 0.0) or 0.0) > 0.0 and enc.get(dd):
                        carriers_by_residue[enc[dd]].add(sp)
            n_orig = max((origin_fn(res, sps) for res, sps in carriers_by_residue.items()), default=0)

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
            "n_eff": effective_number(supports),
            "n_sp": len(species),
            "n_orig": n_orig,
            "convergence_type": _convergence_type(agree_num, agree_den),
            "participating_hyps": sorted({h for hs in hyp_ids_by_domain.values() for h in hs}),
        }

    return {"top": _agg("top"), "bottom": _agg("bottom"), "n_hypotheses": M}
