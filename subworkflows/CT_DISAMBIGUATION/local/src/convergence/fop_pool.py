#!/usr/bin/env python3
# fop_pool.py — Pool the per-hypothesis domain records of a position and scheme into one score per side.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/convergence/

"""
Multi-hypothesis (FOP) pooling of domain scores into one position score per phenotype side.

``pool_domains`` collapses ``M >= 1`` per-hypothesis ``compute_domain_scores``
records for one ``(Gene, Position, scheme)`` into one score per phenotype side.
It is **treeless**: every tree lookup already happened in
``path_scores.compute_domain_scores``; here we only combine scalars.

Imported by: src/convergence/disambiguate_single.py, src/core/driver.py, disambiguation_perms_main.py and
reaggregate_perm_scores.py (`base_cycle`)
Inputs: lists of {"hyp", "sides"} records, the (hypothesis, domain) keys of the design and its species count S per side
Outputs: the pooled per-side aggregates described in `pool_domains`

Per side, over the K fixed Voronoi domains of the design and the pooled hypotheses taken together (a hypothesis
is only a way of drawing the pairs, it carries no weight of its own):

    c(a, b)  = 1 - P(the derived residue was already at the LCA of a and b), the contribution of a pair of
               domains that carry the same derived residue; the distinct species combinations of the same two
               domains give one value each and the pair of domains takes their mean
    s(d)     = 1 - prod over the partners b of d of (1 - c(d, b))          (0 without a partner)
    strength = sum_d s(d) / K
    breadth  = n_sp / S_side, n_sp the distinct species that carry the derived residue in a domain with a positive
               score in at least one hypothesis, S_side the distinct species of that side's group (species1 for
               top, species2 for bottom) in the design of the cycle
    score    = clamp01(strength * breadth)                                 (breadth = 1 when S_side is not given)

A species counts once however many hypotheses use it. ``agree_num`` / ``agree_den`` and ``convergence_type`` come
from the modal encoded derived residue per domain across the pooled hypotheses (a domain that splits V/I/L
between hypotheses scores low under US and high under a scheme that co-encodes them; the encoding is already
applied in each ``domain_der_enc`` upstream, so ``pool_domains`` never sees the scheme).

See ``docs/scoring_v3_core.md`` for the scoring contract.
"""

from __future__ import annotations

import collections

import math
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

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


def _pair_contributions(hyps: List[Dict], side: str) -> Dict[Tuple[str, str], List[float]]:
    """{(domain_a, domain_b): [c, ...]} over the distinct (domain, species) combinations of the pooled hypotheses.

    A combination that several hypotheses share counts once; a pair of domains without the species of both
    ends (no row of `domain_der_species`) cannot be told apart from a repeat and is left out."""
    comb: Dict[Tuple, float] = {}
    for r in hyps:
        srow = (r.get("sides") or {}).get(side) or {}
        spd = {str(d): sp for d, sp in (srow.get("domain_der_species") or {}).items()}
        for a, b, _lca, c in srow.get("pair_lca") or []:
            a, b = str(a), str(b)
            if a in spd and b in spd:
                ends = tuple(sorted(((a, spd[a]), (b, spd[b]))))
                comb[ends] = float(c)
    by_pair: Dict[Tuple[str, str], List[float]] = collections.defaultdict(list)
    for ((a, _sa), (b, _sb)), c in comb.items():
        by_pair[(a, b)].append(c)
    return by_pair


def pool_domains(
    hyp_records: List[Dict],
    design_domains: Optional[Iterable[Tuple[str, Any]]] = None,
    n_design_species: Optional[Dict[str, int]] = None,
) -> Dict[str, Any]:
    """Treeless pooler over the fixed Voronoi domains.

    Args:
        hyp_records: ``list[{"hyp": str, "sides": <compute_domain_scores return>}]``,
            ``M >= 1``.
        design_domains: the ``(hypothesis, domain)`` keys of the design (every hypothesis it holds, not only the
            ones that detect the position); its domains complete the K of the design when a record lacks them.
        n_design_species: S per side ``{"top", "bottom"}``, the distinct species of that phenotype group in the
            design of the cycle; a side without S (``None``, missing or <= 0) takes breadth 1.

    Returns:
        ``{"top": <agg>, "bottom": <agg>, "n_hypotheses": M}`` where ``<agg>`` is
        ``{asr_path_score, domain_scores, domain_der, domain_anc, agree_num, agree_den, n_participating,
        convergence_type, n_sp}``. ``domain_scores`` is s(d) of every domain of the design;
        ``domain_der`` / ``domain_anc`` carry only domains changed in >= 1
        hypothesis (modal residue); ``agree_den == n_participating``; ``agreement_tie`` flags a changed
        domain whose derived residue was tied; ``participating_hyps`` lists the
        hypotheses that changed >= 1 domain on the side.
    """
    hyps = [r for r in hyp_records if r.get("hyp")]
    M = len(hyps)

    # Domain universe: the domain_meta ids (all K domains) plus any ids only the design keys have.
    rep: Dict[str, Any] = {}
    for r in hyps:
        for d in ((r.get("sides") or {}).get("domain_meta") or {}):
            rep.setdefault(str(d), d)
    for (_h, d) in (design_domains or ()):
        rep.setdefault(str(d), d)
    K = len(rep)

    def _agg(side: str) -> Dict[str, Any]:
        # Noisy-OR of the contributions of each domain's partners; the pair of two domains counts the mean of its
        # distinct species combinations.
        partners: Dict[str, List[float]] = collections.defaultdict(list)
        for (a, b), cs in _pair_contributions(hyps, side).items():
            c = math.fsum(cs) / len(cs)
            partners[a].append(c)
            partners[b].append(c)
        s_dom = {sd: (1.0 - math.prod(1.0 - c for c in partners[sd])) if partners.get(sd) else 0.0 for sd in rep}
        strength = (math.fsum(s_dom.values()) / K) if K else 0.0

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

        # Species that carry the derived residue in a domain that converges (a positive domain score) in at least
        # one hypothesis: the independent observations behind the position, whatever the hypotheses that used them.
        species = set()
        for r in hyps:
            srow = (r.get("sides") or {}).get(side) or {}
            ds = srow.get("domain_scores") or {}
            for dd, sp in (srow.get("domain_der_species") or {}).items():
                if float(ds.get(dd, 0.0) or 0.0) > 0.0:
                    species.add(sp)
        S = int((n_design_species or {}).get(side) or 0)
        breadth = min(1.0, len(species) / S) if S > 0 else 1.0
        core = max(0.0, min(1.0, strength * breadth))

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
            "domain_scores": {rep[sd]: v for sd, v in s_dom.items()},
            "domain_der": {rep.get(sd, sd): _modal_str(der_raw.get(sd, [])) for sd in changed},
            "domain_anc": {rep.get(sd, sd): _modal_str(anc_raw.get(sd, [])) for sd in changed},
            "domain_der_support": {rep.get(sd, sd): _support_str(der_raw.get(sd, [])) for sd in changed},
            "domain_anc_support": {rep.get(sd, sd): _support_str(anc_raw.get(sd, [])) for sd in changed},
            "agree_num": agree_num,
            "agree_den": agree_den,
            "agreement_tie": agreement_tie,
            "n_participating": agree_den,
            "n_sp": len(species),
            "convergence_type": _convergence_type(agree_num, agree_den),
            "participating_hyps": sorted({h for hs in hyp_ids_by_domain.values() for h in hs}),
        }

    return {"top": _agg("top"), "bottom": _agg("bottom"), "n_hypotheses": M}
