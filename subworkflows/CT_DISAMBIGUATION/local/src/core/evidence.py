# evidence.py — Per-domain evidence rows of a position and the selection of the top positions to explain.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
Evidence of a position: what each domain of each hypothesis saw, and the choice of the positions worth looking at.

Imported by: explain_positions.py (run by CT_EVIDENCE, subworkflows/CT_DISAMBIGUATION/ct_evidence.nf)
Inputs: unpooled scorer rows; position_scores.tsv rows (Gene, Position, side, CAAS_score, p.emp)
Outputs: EVIDENCE_COLUMNS rows (dicts of text) and the ranked top positions; the files are written by the caller

`evidence_rows` turns the rows of the scorer before the hypotheses of a position are pooled
(`analyze_gene_disambiguation(..., keep_unpooled=True)`) into one table row per entry and domain: the MRCA of the domain's
pair with its reconstructed state and posterior, the species and residues at the tips of both sides, the score of the
domain on each side and the LCA of the domains that share a residue. The master holds the pooled numbers; this is how
they came about. `select_top_positions` picks the positions to explain from position_scores.tsv.
"""
import math
import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

EVIDENCE_COLUMNS = [
    "gene", "msa_pos", "caap_group", "hypothesis", "tag", "caas", "domain",
    "mrca_node", "mrca_state", "mrca_posterior",
    "top_species", "bottom_species", "top_tip_aa", "bottom_tip_aa",
    "top_domain_score", "bottom_domain_score", "top_pair_lca", "bottom_pair_lca",
]


def _text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def _by_text(mapping: Optional[Dict[Any, Any]]) -> Dict[str, Any]:
    """The keys of a per-domain mapping as text: the scorer keys domains by int, a round trip through JSON by str."""
    return {str(k): v for k, v in (mapping or {}).items()}


def _hypothesis_number(label: Any) -> int:
    m = re.search(r"\d+", str(label or ""))
    return int(m.group()) if m else 0


def _lca_of_domain(pair_lca: Iterable[Any], domain: str) -> str:
    """The LCA entries ([a, b, lca_node, contribution]) that involve the domain, as 'domain-other:lca:contribution'."""
    parts = []
    for a, b, lca, contrib in pair_lca or []:
        if lca is None or domain not in (str(a), str(b)):
            continue
        other = b if str(a) == domain else a
        parts.append(f"{domain}-{other}:{lca}:{float(contrib):.4f}")
    return "|".join(parts)


def evidence_rows(unpooled: Iterable[Any]) -> List[Dict[str, str]]:
    """One row per unpooled scorer row and domain (EVIDENCE_COLUMNS), ordered by position, caap_group, hypothesis
    number and domain. Every domain of the design is listed, also the ones that did not change (score and tips empty)."""
    rows: List[Dict[str, str]] = []
    for r in unpooled:
        sides = getattr(r, "sides", None) or {}
        meta = _by_text(sides.get("domain_meta"))
        pairs = {str(p["pair_id"]): p for p in (getattr(r, "pair_details", None) or [])
                 if isinstance(p, dict) and p.get("pair_id") is not None}
        side = {s: sides.get(s) or {} for s in ("top", "bottom")}
        scores = {s: _by_text(side[s].get("domain_scores")) for s in side}
        domains = set(meta) | set(pairs) | set(scores["top"]) | set(scores["bottom"])
        for d in sorted(domains, key=lambda x: (0, int(x)) if x.isdigit() else (1, x)):
            m, p = meta.get(d) or {}, pairs.get(d) or {}
            mrca = m.get("mrca_id") if m.get("mrca_id") is not None else p.get("node_id")
            state = m.get("state") if m.get("state") is not None else p.get("focal_state")
            posterior = m.get("posterior") if m.get("posterior") is not None else p.get("focal_prob")
            rows.append({
                "gene": _text(getattr(r, "gene", "")),
                "msa_pos": _text(getattr(r, "position", "")),
                "caap_group": _text(getattr(r, "caap_group", "US") or "US"),
                "hypothesis": _text(getattr(r, "hypothesis", "")),
                "tag": _text(getattr(r, "tag", "")),
                "caas": _text(getattr(r, "caas", "")),
                "domain": d,
                "mrca_node": _text(mrca),
                "mrca_state": _text(state),
                "mrca_posterior": _text(float(posterior)) if posterior not in (None, "") else "",
                "top_species": _text(p.get("top_species")),
                "bottom_species": _text(p.get("bottom_species")),
                "top_tip_aa": _text(p.get("top_tip_mode") or p.get("top_tip_residue")),
                "bottom_tip_aa": _text(p.get("bottom_tip_mode") or p.get("bottom_tip_residue")),
                "top_domain_score": _text(scores["top"].get(d)),
                "bottom_domain_score": _text(scores["bottom"].get(d)),
                "top_pair_lca": _lca_of_domain(side["top"].get("pair_lca"), d),
                "bottom_pair_lca": _lca_of_domain(side["bottom"].get("pair_lca"), d),
            })
    rows.sort(key=lambda x: (int(x["msa_pos"]), x["caap_group"], _hypothesis_number(x["hypothesis"]),
                             (0, int(x["domain"])) if x["domain"].isdigit() else (1, x["domain"])))
    return rows


def _number(value: Any) -> Optional[float]:
    if value in (None, "", "NA"):
        return None
    x = float(value)
    return None if math.isnan(x) else x


def select_top_positions(score_rows: Iterable[Dict[str, Any]], n: int) -> List[Tuple[str, str, Dict[str, Optional[float]]]]:
    """The n best positions of position_scores.tsv rows (mappings with Gene, Position, side, CAAS_score and p.emp).

    A position is ranked by the best CAAS_score over its sides; ties go to the smaller p.emp (NA last), then to the gene
    name and the position, so the choice is the same on every run. A position with no score is not a candidate.
    Returns (gene, position, {"score": ..., "p_emp": ...}) in rank order; n <= 0 gives none."""
    if n <= 0:
        return []
    best: Dict[Tuple[str, str], Dict[str, Optional[float]]] = {}
    for row in score_rows:
        score = _number(row.get("CAAS_score"))
        if score is None:
            continue
        key = (str(row["Gene"]), str(row["Position"]))
        p = _number(row.get("p.emp"))
        cur = best.setdefault(key, {"score": score, "p_emp": p})
        cur["score"] = max(cur["score"], score)
        if p is not None and (cur["p_emp"] is None or p < cur["p_emp"]):
            cur["p_emp"] = p
    ranked = sorted(best.items(), key=lambda kv: (-kv[1]["score"], math.inf if kv[1]["p_emp"] is None else kv[1]["p_emp"],
                                                  kv[0][0], int(kv[0][1])))
    return [(g, p, info) for (g, p), info in ranked[:n]]
