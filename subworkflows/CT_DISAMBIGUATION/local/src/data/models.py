# models.py — Dataclasses of the CAAS entries and the convergence results of the disambiguation.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/data/

"""
Data models for convergence-type disambiguation outputs.

- CAASPosition: one discovery row (position, substitution, scheme, discovering hypothesis).
- ConvergenceResult: one scored row (position, scheme and, per side, the pooled ASR score and diagnostics).

Imported by: src/convergence/disambiguate_single.py, src/core/observed.py, src/utils/gene_wrapper.py
"""

__all__ = [
    "CAASPosition",
    "ConvergenceResult",
    "BiochemResults",
]

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class CAASPosition:
    """Container for CAAS position metadata."""

    position: int          # 0-based MSA column index
    position_one_based: int
    tag: str
    caas: str
    trait1_aa: List[str] = field(
        default_factory=list
    )  # High phenotype amino acids (trait=1)
    trait0_aa: List[str] = field(
        default_factory=list
    )  # Low phenotype amino acids (trait=0)
    is_focus: bool = False
    caap_group: str = "US"
    amino_encoded: str = ""
    is_conserved_meta: bool = False
    conserved_pair: str = ""
    # Discovering hypothesis of this row: the traitfile token (e.g. "traitfile_H5.tab"
    # or "H5") when the design has several hypotheses, empty for a single trait file.
    # The disambiguation scores the row against that hypothesis's contrast pairs
    # only, never against the union across hypotheses.
    trait: str = ""


@dataclass
class ConvergenceResult:
    """
    Container for ASR-driven convergence analysis results.
    """

    # Core identification
    gene: str
    position: int
    tag: str
    caas: str

    # State information (from CAAS metadata or ASR)
    ancestral: str  # displayed in the progress log
    derived: str  # displayed in the progress log

    # Pattern classification
    convergence_type: str

    # Tip-level pattern analysis
    trait1_aa: List[str] = field(default_factory=list)
    trait0_aa: List[str] = field(default_factory=list)
    tip_pattern_comment: Optional[str] = None
    pair_details: Optional[List[dict]] = None
    caap_group: str = "US"
    amino_encoded: str = ""
    # Discovering hypothesis ("H<n>"; None for a single-contrast run). One
    # ConvergenceResult is emitted per (position, scheme, hypothesis).
    hypothesis: Optional[str] = None
    # Hypotheses pooled for this position and scheme that drove >= 1 changed
    # domain on THIS SIDE, comma-joined: the only hypothesis-provenance field
    # downstream (scoring_compute.R reads it by name). Per side, so the top and
    # bottom rows of a position can differ; never nulled by a multi-hypothesis
    # pool (unlike `hypothesis`).
    participating_hypotheses: Optional[str] = None
    # Number of hypotheses (M, from fop_pool.pool_domains) pooled for this
    # (position, scheme); identical on both emitted side rows.
    n_hypotheses: Optional[int] = None
    # Cross-hypothesis support tallies for the fields above that otherwise take
    # the value of the first hypothesis row in file order (see
    # `_emit_pooled_side_rows`). They sit beside those fields and do not replace them.
    tag_support: str = ""
    caas_support: str = ""
    amino_encoded_support: str = ""

    # Node mapping and state tracking
    node_mapping: Optional[Dict[str, int]] = None
    asr_ancestral_state: Optional[str] = None
    asr_descendant_states: Optional[List[str]] = None
    node_state_details: Optional[Dict[str, Any]] = None
    node_posteriors: Optional[Dict[str, Any]] = None

    # Root/MRCA states
    root_state: Optional[str] = None
    mrca_state: Optional[str] = None
    focal_states: Optional[Dict[str, Optional[str]]] = None
    node_state_summary: Optional[Dict[str, Optional[str]]] = None
    state_source: str = "unknown"

    # Scoring and quality
    score: Optional[Any] = None
    position_one_based: Optional[int] = None

    # Change tracking
    is_focus: bool = False
    # Direction key (top / bottom / none). A position that changes on both sides
    # is emitted as TWO ConvergenceResult rows keyed (gene, position, side), each
    # carrying that direction's own asr_path_score / derived_agreement /
    # convergence_type. A position with no participating pair is one row with
    # side == "none".
    side: str = "none"

    # CAAS convergence score of the side: the pooled mean over the Voronoi domains
    # (computed in src/convergence/path_scores.py, pooled in src/convergence/fop_pool.py).
    asr_path_score: Optional[float] = None
    # Diagnostic only: agree_num / agree_den (largest same-encoded-residue group
    # over the domains changed in >= 1 hypothesis / count of those domains).
    derived_agreement: Optional[float] = None
    # True when a changed domain's derived residue was tied for the highest support and settled by a convention
    # (smallest residue), so derived_agreement and convergence_type may rest on that choice.
    agreement_ambiguous: Optional[bool] = None
    # Per-domain pooled score s̄_d for the emitted side.
    domain_scores: Optional[Dict[int, float]] = None
    # Raw (un-encoded) ancestral and per-side derived residues per changed domain,
    # modal over the pooled hypotheses. Flattened to domain_<d>_anc_aa / _top_aa / _bot_aa.
    domain_anc_aa: Optional[Dict[int, Any]] = None
    domain_der_top_aa: Optional[Dict[int, str]] = None
    domain_der_bot_aa: Optional[Dict[int, str]] = None
    # Cross-hypothesis support tallies (e.g. "I:1,V:1") for the modal residues
    # above. They sit beside those fields: domain_<d>_top_aa / _bot_aa / _anc_aa
    # hold the modal winner only.
    domain_der_support_top_aa: Optional[Dict[Any, str]] = None
    domain_der_support_bot_aa: Optional[Dict[Any, str]] = None
    domain_anc_support_aa: Optional[Dict[Any, str]] = None
    # All domains: {d: {"mrca_id", "state", "posterior"}}, the node, ancestral
    # state and posterior of each domain. Domains without a reconstruction have
    # state None, posterior 0.0.
    domain_meta: Optional[Dict[int, Dict[str, Any]]] = None
    # Union across pooled hypotheses of (domain_a, domain_b, lca_node_id,
    # contrib) for the same-residue domain pairs that drove this side's
    # `asr_path_score` (see path_scores.score_domains_side); duplicate (a, b, lca)
    # triples across hypotheses are averaged on `contrib`.
    pair_lca: Optional[List[Tuple[Any, Any, int, float]]] = None


# Alias of ConvergenceResult, exported as BiochemResults
BiochemResults = ConvergenceResult
