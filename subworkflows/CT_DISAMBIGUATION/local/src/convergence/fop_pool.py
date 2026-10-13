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
Inputs: lists of {"hyp", "sides"} records and the design pairs {(hypothesis, domain): (top species, bottom species, pss)}
Outputs: the pooled per-side aggregates described in `pool_domains`

Per side, over the K fixed Voronoi domains of the design. A hypothesis only draws which pair represents each
domain, so the unit of a domain is the species that represent it, not the hypotheses:

    s_x    = mean, over the hypotheses of the design in which the species x represents the domain, of the domain
             score of the hypothesis (0 where the hypothesis did not detect the position)
    w_x    = mean PSS of the pairs in which x represents the domain (1 without PSS)
    S_d    = sum_x w_x s_x / sum_x w_x          over the distinct species x of the side that represent the domain
    score  = clamp01( sum_d S_d / K )

A species counts once however many hypotheses use it, and a domain weighs the same as any other: the weights only
distribute the weight of a domain among its species, so a constant factor of the PSS does not change the score.
``agree_num`` / ``agree_den`` and ``convergence_type`` come from the modal encoded derived residue per domain across
the pooled hypotheses (a domain that splits V/I/L between hypotheses scores low under US and high under a scheme that
co-encodes them; the encoding is already applied in each ``domain_der_enc`` upstream, so ``pool_domains`` never sees
the scheme).

See ``docs/scoring_v3_core.md`` for the scoring contract.
"""

from __future__ import annotations

import collections

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


def _domain_scores(hyps: List[Dict], side: str, rep: Dict[str, Any], design: Dict[Tuple[str, str], Tuple]) -> Dict[str, float]:
    """S_d of every domain of the side (see the module docstring); 0 for a domain without a species in the design."""
    col = 0 if side == "top" else 1
    by_hyp: Dict[str, Dict[str, float]] = {}
    for r in hyps:
        ds = ((r.get("sides") or {}).get(side) or {}).get("domain_scores") or {}
        by_hyp[r["hyp"]] = {str(d): float(v or 0.0) for d, v in ds.items()}
    # species -> ([score of each hypothesis where it represents the domain], [pss of each of those pairs])
    per: Dict[str, Dict[str, Tuple[List[float], List[float]]]] = collections.defaultdict(lambda: collections.defaultdict(lambda: ([], [])))
    for (h, d), pair in design.items():
        pss = pair[2]
        scores, weights = per[d][pair[col]]
        scores.append(by_hyp.get(h, {}).get(d, 0.0))
        weights.append(1.0 if pss is None else max(float(pss), 0.0))
    out: Dict[str, float] = {}
    for sd in rep:
        species = per.get(sd) or {}
        if not species:
            out[sd] = 0.0
            continue
        s_x = {x: math.fsum(sc) / len(sc) for x, (sc, _w) in species.items()}
        w_x = {x: math.fsum(w) / len(w) for x, (_sc, w) in species.items()}
        total = math.fsum(w_x.values())
        if total <= 0.0:   # no pair with a positive PSS: the species weigh the same
            w_x, total = {x: 1.0 for x in species}, float(len(species))
        out[sd] = math.fsum(w_x[x] * s_x[x] for x in species) / total
    return out


def pool_domains(
    hyp_records: List[Dict],
    design_pairs: Dict[Tuple[str, Any], Tuple[str, str, Optional[float]]],
) -> Dict[str, Any]:
    """Treeless pooler over the fixed Voronoi domains.

    Args:
        hyp_records: ``list[{"hyp": str, "sides": <compute_domain_scores return>}]``,
            ``M >= 1``: the hypotheses that detected the position.
        design_pairs: ``{(hypothesis, domain): (top species, bottom species, pss)}`` of the design (every hypothesis it
            holds, not only the ones that detect the position); ``pss`` may be ``None`` (every pair weighs the same).
            When the design has a single hypothesis, every record belongs to it whatever its label.

    Returns:
        ``{"top": <agg>, "bottom": <agg>, "n_hypotheses": M}`` where ``<agg>`` is
        ``{asr_path_score, domain_scores, domain_der, domain_anc, agree_num, agree_den, n_participating,
        convergence_type, n_sp}``. ``domain_scores`` is S_d of every domain of the design;
        ``domain_der`` / ``domain_anc`` carry only domains changed in >= 1
        hypothesis (modal residue); ``agree_den == n_participating``; ``agreement_tie`` flags a changed
        domain whose derived residue was tied; ``participating_hyps`` lists the
        hypotheses that changed >= 1 domain on the side; ``n_sp`` counts the distinct species that carry the derived
        residue in a domain with a positive score in at least one hypothesis.
    """
    if not design_pairs:
        raise ValueError("pool_domains needs the design pairs {(hypothesis, domain): (top species, bottom species, pss)}")
    design = {(h, str(d)): pair for (h, d), pair in design_pairs.items()}
    design_hyps = {h for h, _d in design}
    sole = next(iter(design_hyps)) if len(design_hyps) == 1 else None
    # a record of a single-hypothesis design may carry any label (the observed rows of a one-contrast run have none)
    hyps = [dict(r, hyp=sole) if (sole is not None and r.get("hyp") not in design_hyps) else r for r in hyp_records if r.get("hyp")]
    M = len(hyps)

    # Domain universe: the domain_meta ids (all K domains) plus any ids only the design has.
    rep: Dict[str, Any] = {}
    for r in hyps:
        for d in ((r.get("sides") or {}).get("domain_meta") or {}):
            rep.setdefault(str(d), d)
    for (_h, d) in design:
        rep.setdefault(d, d)
    K = len(rep)

    def _agg(side: str) -> Dict[str, Any]:
        s_dom = _domain_scores(hyps, side, rep, design)
        core = max(0.0, min(1.0, math.fsum(s_dom.values()) / K)) if K else 0.0

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
