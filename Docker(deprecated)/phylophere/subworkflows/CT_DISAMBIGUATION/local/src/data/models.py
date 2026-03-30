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

    position: int
    position_zero_based: int
    position_one_based: int
    tag: str
    caas: str
    trait1_aa: List[str] = field(
        default_factory=list
    )  # High phenotype amino acids (trait=1)
    trait0_aa: List[str] = field(
        default_factory=list
    )  # Low phenotype amino acids (trait=0)
    pvalue: Optional[float] = None
    pvalue_boot: Optional[float] = None
    is_significant: bool = False
    is_focus: bool = False
    caap_group: str = "US"
    amino_encoded: str = ""
    is_conserved_meta: bool = False
    conserved_pair: str = ""
    sig_hyp: Optional[bool] = None
    sig_perm: Optional[bool] = None
    sig_both: Optional[bool] = None


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
    is_significant: bool

    # State information (from CAAS metadata or ASR)
    ancestral: str  # Used for display/legacy compatibility
    derived: str  # Used for display/legacy compatibility

    # Pattern classification
    pattern_type: str
    convergence_description: str

    # Fields with defaults (must come after non-default fields)
    convergence_mode: str = "focal_clade"  # 'mrca' or 'focal_clade'

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
    sig_hyp: Optional[bool] = None
    sig_perm: Optional[bool] = None
    sig_both: Optional[bool] = None

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
    position_zero_based: Optional[int] = None
    position_one_based: Optional[int] = None
    caas_pvalue: Optional[float] = None
    pvalue_boot: Optional[float] = None
    low_confidence_nodes: Optional[List[str]] = None

    # Change tracking
    is_focus: bool = False
    top_change_type: str = "none"
    bottom_change_type: str = "none"
    change_side: str = "none"

    # Ambiguity and conserved-pair flags (dynamic multi-pair)
    ambiguous: bool = False
    asr_is_conserved: bool = False


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
