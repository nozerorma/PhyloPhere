"""Data models for convergence-type disambiguation outputs.

Provided classes
----------------
- CAASPosition: tip/position metadata and significance flags.
- ConvergenceResult: consolidated ASR/convergence results and diagnostics.

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-07
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
    # Discovering hypothesis for this row. FOP runs carry the source traitfile
    # token (e.g. "…/traitfile_H5.tab"); a single non-FOP traitfile leaves this
    # empty. Disambiguation uses it to score the row against that hypothesis's
    # contrast pairs only, never the union across hypotheses.
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
    ancestral: str  # Used for display/legacy compatibility
    derived: str  # Used for display/legacy compatibility

    # Pattern classification
    convergence_type: str

    # Tip-level pattern analysis
    trait1_aa: List[str] = field(default_factory=list)
    trait0_aa: List[str] = field(default_factory=list)
    tip_pattern_comment: Optional[str] = None
    pair_details: Optional[List[dict]] = None
    caap_group: str = "US"
    amino_encoded: str = ""
    # Discovering hypothesis (FOP: "H<n>"; single-contrast run: None). One
    # ConvergenceResult is emitted per (position, scheme, hypothesis).
    hypothesis: Optional[str] = None
    # Hypotheses (across the harvest) that drove >= 1 changed domain on THIS
    # SIDE, comma-joined -- the sole hypothesis-provenance field downstream
    # consumes (SCORING's pos_scores aggregation reads it by name). Per-side,
    # so top/bottom rows for the same position can legitimately differ; never
    # nulled by a multi-hypothesis pool (unlike `hypothesis`).
    participating_hypotheses: Optional[str] = None
    # Harvest size (M, from fop_pool.pool_domains) for this (position, scheme)
    # pool -- position/scheme level, identical on both emitted side rows.
    n_hypotheses: Optional[int] = None
    # Cross-hypothesis support tallies for fields that otherwise silently pass
    # through an arbitrary first-hypothesis-in-file-order row (see
    # `_emit_pooled_side_rows`). Kept alongside the status-quo passthrough
    # fields, not replacing them.
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
    # First-class direction key (top / bottom / none). T4b retired the
    # change_top/change_bottom/change_side triplet: a "both" position is always
    # emitted as TWO ConvergenceResult rows keyed (gene, position, side), each
    # carrying that direction's own asr_path_score / derived_agreement /
    # convergence_type. A position with no participating pair is one row with
    # side == "none".
    side: str = "none"

    # CAAS convergence score on the Voronoi domain (scoring_v2 core v3; computed
    # in src/convergence/path_scores.py + pooled in src/convergence/fop_pool.py).
    # The per-side pooled domain mean (formerly duplicated onto a `core` field --
    # retired since it was always bit-identical to this one).
    asr_path_score: Optional[float] = None
    # Diagnostic only: agree_num / agree_den (largest same-encoded-residue group
    # over the domains changed in >= 1 hypothesis / count of those domains).
    derived_agreement: Optional[float] = None
    # Per-domain pooled score s̄_d for the emitted side (was pair_path_scores).
    domain_scores: Optional[Dict[int, float]] = None
    # Raw (un-encoded) ancestral + per-side derived residues per changed domain,
    # modal over the harvest. Flattened to domain_<d>_anc_aa / _top_aa / _bot_aa.
    domain_anc_aa: Optional[Dict[int, Any]] = None
    domain_der_top_aa: Optional[Dict[int, str]] = None
    domain_der_bot_aa: Optional[Dict[int, str]] = None
    # Cross-hypothesis support tallies (e.g. "I:1,V:1") for the modal residues
    # above -- siblings, not replacements: domain_N_top_aa/bot_aa/anc_aa stay
    # the modal winner CT_POSTPROC's residue_descriptors.py reads directly.
    domain_der_support_top_aa: Optional[Dict[Any, str]] = None
    domain_der_support_bot_aa: Optional[Dict[Any, str]] = None
    domain_anc_support_aa: Optional[Dict[Any, str]] = None
    # All K domains: {d: {"mrca_id", "state", "posterior"}} — carries the node /
    # ancestral state / posterior per domain (was the mrca_<i>_node/state/posterior
    # block sourced from node_mapping). Domains without a reconstruction have
    # state None, posterior 0.0.
    domain_meta: Optional[Dict[int, Dict[str, Any]]] = None
    # Union across pooled hypotheses of (domain_a, domain_b, lca_node_id,
    # contrib) for the same-residue domain pairs that drove this side's
    # `asr_path_score` (see path_scores.score_domains_side); duplicate (a, b, lca) triples across
    # hypotheses are averaged on `contrib`. Debug-tree visualization only.
    pair_lca: Optional[List[Tuple[Any, Any, int, float]]] = None


# Backward-compatible alias
BiochemResults = ConvergenceResult
