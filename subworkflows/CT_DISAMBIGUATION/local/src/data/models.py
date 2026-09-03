"""Data models for convergence-type disambiguation outputs.

Provided classes
----------------
- CAASPosition: tip/position metadata and significance flags.
- ConvergenceResult: consolidated ASR/convergence results and diagnostics.
- ContrastDefinition: species/contrast definitions and tip residue holders.

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
    "ContrastDefinition",
]

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


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
    recovery_boot: Optional[float] = None
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
    pair_transition_summary: Optional[List[dict]] = None
    caap_group: str = "US"
    amino_encoded: str = ""
    is_conserved_meta: bool = False
    conserved_pair: str = ""
    # Discovering hypothesis (FOP: "H<n>"; single-contrast run: None). One
    # ConvergenceResult is emitted per (position, scheme, hypothesis).
    hypothesis: Optional[str] = None

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

    # Derived state analysis
    derived_similarity: Optional[Dict[str, Any]] = None

    # Scoring and quality
    score: Optional[Any] = None
    position_one_based: Optional[int] = None
    pvalue: Optional[float] = None
    recovery_boot: Optional[float] = None

    # Change tracking
    is_focus: bool = False
    change_top: str = "no_change"
    change_bottom: str = "no_change"
    change_side: str = "none"

    # ASR path score (unified replacement for binary ASR gate + convergence +
    # parallel; computed in src/convergence/path_scores.py)
    asr_path_score: Optional[float] = None
    independence: Optional[float] = None
    mrca_diversity: Optional[float] = None
    derived_agreement: Optional[float] = None
    conservation_gate: Optional[float] = None
    core: Optional[float] = None
    pair_path_scores: Optional[Dict[int, float]] = None
    pair_path_contaminated: Optional[Dict[int, bool]] = None
    # Per-conserved-pair conservation-to-root + MRCA node id, analogous to
    # pair_path_scores/pair_path_contaminated. Used by the FOP domain-pooler to
    # rebuild conservation_gate from the DISTINCT conserved pairs shared across
    # hypotheses (dedup by node) rather than averaging per-hypothesis gates.
    conserved_pair_path_scores: Optional[Dict[int, float]] = None
    conserved_pair_path_nodes: Optional[Dict[int, Any]] = None
    # Raw (un-encoded) ancestral + per-side derived residues per changed pair.
    # Plumbed to the flat TSV as mrca_<i>_anc_aa / _top_aa / _bot_aa so the FOP
    # domain-pooler can recompute derived_agreement HARVEST-WIDE and PER SCHEME
    # (POINT 3) instead of pooling per-hypothesis derived_agreement.
    pair_ancestral_aa: Optional[Dict[int, Any]] = None
    pair_derived_top_aa: Optional[Dict[int, str]] = None
    pair_derived_bot_aa: Optional[Dict[int, str]] = None


@dataclass
class ContrastDefinition:
    """Container for species contrast definitions."""

    pair_id: str
    top_taxa: List[str]
    bottom_taxa: List[str]
    top_species: List[str]
    bottom_species: List[str]
    all_taxa: List[str]
    node_id: Optional[int] = None
    top_tip_residues: Optional[List[dict]] = None
    bottom_tip_residues: Optional[List[dict]] = None
    top_tip_mode: Optional[str] = None
    bottom_tip_mode: Optional[str] = None
    top_tip_residue: Optional[str] = None
    bottom_tip_residue: Optional[str] = None
    focal_state: Optional[str] = None
    mrca_contrast: Optional[str] = None
    mrca_modal_aa: Optional[str] = None


# Backward-compatible alias
BiochemResults = ConvergenceResult
