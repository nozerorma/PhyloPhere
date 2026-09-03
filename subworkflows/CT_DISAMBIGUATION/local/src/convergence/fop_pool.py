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
* ``core   = P(>=2 domains carry a clean change)``  (inclusion-exclusion, the
  same ``_p_at_least_2`` as path_scores.py, over the pooled ``c_d``);
* ``independence`` / ``mrca_diversity`` / ``derived_agreement`` /
  ``conservation_gate`` are pooled across hypotheses by the same PSS weights
  (they are carried on the null's PositionAxes record for exactly this);
* recombine with the path_scores.py algebra, verbatim from fop_pool.R::

      diversity_mult = 0.75 + 0.25 * mrca_diversity_pooled
      replication    = independence_pooled * core_pooled
      strength       = diversity_mult * derived_agreement_pooled
      asr_pooled     = replication * strength * conservation_gate_pooled

Single hypothesis (or non-FOP) -> the lone ``asr_path_score`` passes through.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple


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
) -> Dict[str, float]:
    """Domain-pool one position's FOP hypothesis records.

    Args:
        hyp_records: list of dicts, one per hypothesis, each with keys
            ``hyp`` (e.g. "H3"), ``asr_path_score`` (float), and
            ``pair_scores`` ({domain_id:int -> s(p, site):float}).
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

    # Per-hypothesis weight = min pair PSS across its domains (a hypothesis is
    # only as strong as its weakest contrast — mirrors the FOP ranking key and
    # fop_pool.R).
    def hyp_weight(h: str) -> Optional[float]:
        if not pss_by_hyp_domain:
            return None
        ws = [v for (hh, _d), v in pss_by_hyp_domain.items()
              if hh == h and v is not None and v == v]
        return min(ws) if ws else None

    hw = {r["hyp"]: hyp_weight(r["hyp"]) for r in hyps}
    row_w = [hw.get(r["hyp"]) for r in hyps]

    # c_d = PSS-weighted mean of s(p, site) over the pairs domain d contributed.
    domains = sorted({d for r in hyps for d in (r.get("pair_scores") or {}).keys()})
    c_vec: List[float] = []
    for d in domains:
        xs, ws = [], []
        for r in hyps:
            s = (r.get("pair_scores") or {}).get(d)
            if s is None:
                continue
            xs.append(float(s))
            ws.append(hw.get(r["hyp"]))
        m = _wmean(xs, ws)
        if m is not None:
            c_vec.append(m)

    core_pooled = p_at_least_2(c_vec)

    def _pool(field: str, default: float) -> float:
        v = _wmean([_num(r.get(field)) for r in hyps], row_w)
        return v if v is not None else default

    indep = _pool("independence", 1.0)
    div = _pool("mrca_diversity", 0.0)
    da = _pool("derived_agreement", 1.0)
    cg = _pool("conservation_gate", 1.0)

    diversity_mult = 0.75 + 0.25 * div
    asr_pooled = max(0.0, min(1.0, indep * core_pooled * diversity_mult * da * cg))

    return {
        "asr_path_score": asr_pooled,
        "core": core_pooled,
        "independence": indep,
        "mrca_diversity": div,
        "derived_agreement": da,
        "conservation_gate": cg,
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
