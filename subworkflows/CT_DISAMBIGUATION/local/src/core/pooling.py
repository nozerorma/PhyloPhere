"""What a domain-pooled record says per phenotype side: one rule for the observed and the null.

`fop_pool.pool_domains` returns {"top": agg, "bottom": agg, "n_hypotheses": M}. A side takes part when at
least one domain changed on it; a position where no side does is reported as a single `side="none"` row.
The observed rows (`ConvergenceResult`) and the null rows (`PositionAxes`) are built from these summaries,
so the participation rule, the float casts and the agreement ratio exist once.
"""
from typing import Any, Dict, List, Optional

SIDES = ("top", "bottom")


def pooled_sides(pooled: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One summary per participating side of a pool_domains result ([] when no domain changed anywhere).

    Keys: side, asr_path_score (float), derived_agreement (agree_num / agree_den, None without a changed
    domain), participating_hyps (comma-joined or None), domain_scores / domain_anc / domain_der /
    domain_der_support / domain_anc_support (dict or None), and convergence_type when the pool carries it.
    """
    out = []
    for side in SIDES:
        d = (pooled or {}).get(side) or {}
        if int(d.get("n_participating", 0) or 0) <= 0:
            continue
        den = int(d.get("agree_den", 0) or 0)
        summary = {
            "side": side,
            "asr_path_score": float(d.get("asr_path_score", 0.0) or 0.0),
            "derived_agreement": (int(d.get("agree_num", 0) or 0) / den) if den else None,
            "participating_hyps": (",".join(d.get("participating_hyps") or []) or None),
        }
        for key in ("domain_scores", "domain_anc", "domain_der", "domain_der_support", "domain_anc_support"):
            summary[key] = dict(d.get(key) or {}) or None
        if "convergence_type" in d:
            summary["convergence_type"] = d["convergence_type"]
        out.append(summary)
    return out
