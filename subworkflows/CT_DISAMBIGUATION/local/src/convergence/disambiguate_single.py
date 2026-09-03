#!/usr/bin/env python3
"""Convergence-type analysis module for single-gene CAAS disambiguation.

Handles ASR-driven convergence detection including:
- Tip-level residue collection and modal analysis
- Phylogenetic tree traversal and MRCA identification
- Convergence pattern classification (tip-level and node-level)
- Conserved-pair validation from metadata + ASR

Author: Refactored from test_nutm2a_real_caas.py
Date: 2025-11-24
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import logging
from collections import Counter, namedtuple

# Add src to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

# Core imports
from src.convergence.convergence import (
    extract_node_states_from_node_level,
    build_alignment_lookup,
    collect_tip_residues,
    extract_tip_residue,
    format_amino_display,
    classify_change_and_parallelism,
)
from src.asr.tree_parser import get_mrca
from src.data.models import CAASPosition, ConvergenceResult
from src.data.loaders import list_gene_caas_entries, parse_trait_pairs
from src.biochem.grouping import get_grouping_scheme
from src.convergence.path_scores import build_node_index, compute_asr_path_score

logger = logging.getLogger(__name__)


# Lightweight per-position record emitted by the axes-only (permulation) path.
# The permulation null keeps ONLY these five fields (Gene/cycle are added by the
# worker), so there is no reason to assemble the full 40-field ConvergenceResult.
# It exposes the same attribute names the perm worker reads via getattr, so it is
# a drop-in for a ConvergenceResult on that path.
PositionAxes = namedtuple(
    "PositionAxes",
    ["position", "caap_group", "asr_path_score", "change_top", "change_bottom",
     "hypothesis", "pair_scores", "independence", "mrca_diversity",
     "derived_agreement", "conservation_gate", "core",
     "conserved_pair_scores", "conserved_pair_nodes",
     "pair_ancestral", "pair_derived_top", "pair_derived_bot"],
)
# All fields after `change_bottom` are optional. Single-contrast perm replay
# leaves `hypothesis` None and the FOP axis fields None. For the FOP null they
# carry the full compute_asr_path_score output so the aggregation domain-pools
# with the exact same algebra scoring_compute.R §2b / fop_pool.R use on the
# observed side.
PositionAxes.__new__.__defaults__ = (
    None, None, None, None, None, None, None, None, None, None, None, None,
)


def _build_per_node_dist(
    posterior_data: Optional[Dict[int, Dict[int, Dict[str, float]]]],
    paml_site: Optional[int],
) -> Dict[int, Dict[str, float]]:
    """Collect ``node_id -> {aa: posterior}`` for one focal site.

    Mirrors exactly what the full scorer stores in ``node_posteriors["per_node"]``
    (same site key, same ``dict(sorted(...))`` ordering, same skip of empty nodes),
    so both the observed path and the perm replay feed
    :func:`compute_asr_path_score` the same input. The ``dict(sorted(...))`` order
    is retained because :func:`modal_encoded` breaks argmax ties by iteration order.
    """
    per_node_dist: Dict[int, Dict[str, float]] = {}
    if not posterior_data or paml_site is None:
        return per_node_dist
    for node_id, node_sites in posterior_data.items():
        site_probs = node_sites.get(paml_site)
        if not site_probs:  # missing/empty site → node contributes nothing
            continue
        per_node_dist[int(node_id)] = dict(sorted(site_probs.items()))
    return per_node_dist


def _position_axes(
    caas_pos: CAASPosition,
    tree_data,
    posterior_data: Optional[Dict[int, Dict[int, Dict[str, float]]]],
    node_index: Optional[Dict[int, Any]],
    pair_details_list: Optional[List[Dict[str, Any]]],
    per_site_dist_cache: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    walk_cache: Optional[Dict[Any, Any]] = None,
) -> Dict[str, Any]:
    """Reduced ASR-path-score kernel shared by the full scorer and the perm replay.

    Builds ``per_node_dist`` for the focal site directly from the posterior map and
    runs :func:`compute_asr_path_score`. This is the ONLY numeric computation the
    permulation null keeps (``asr_path_score`` + its five axes). Extracting it into
    one helper guarantees the observed (full) path and the axes-only perm path score
    each position through *identical* code — bit-for-bit — instead of two drifting
    copies.

    ``per_site_dist_cache`` (perm replay only): ``per_node_dist`` depends solely on
    ``(posterior_data, site)`` — it is invariant to the grouping scheme (applied
    later inside :func:`compute_asr_path_score`) and to the permuted phenotype. So a
    site that recurs across schemes within a cycle, or across the N permulation
    cycles for a gene, is built once and reused. When supplied the map is keyed by
    ``position_one_based``; the observed path passes ``None`` (each site scored once).
    """
    paml_site = caas_pos.position_one_based

    if per_site_dist_cache is not None:
        per_node_dist = per_site_dist_cache.get(paml_site)
        if per_node_dist is None:
            per_node_dist = _build_per_node_dist(posterior_data, paml_site)
            per_site_dist_cache[paml_site] = per_node_dist
    else:
        per_node_dist = _build_per_node_dist(posterior_data, paml_site)

    if node_index is None:  # per-gene invariant; hoisted by the caller when available
        node_index = build_node_index(getattr(tree_data, "root", None))

    return compute_asr_path_score(
        pair_details=pair_details_list,
        per_node_dist=per_node_dist,
        node_index=node_index,
        scheme=getattr(caas_pos, "caap_group", "US") or "US",
        is_conserved_meta=bool(getattr(caas_pos, "is_conserved_meta", False)),
        conserved_pair=str(getattr(caas_pos, "conserved_pair", "") or "").strip(),
        walk_cache=walk_cache,
        site_key=paml_site,
    )


def analyze_caas_position_disambiguation(
    gene: str,
    caas_pos: CAASPosition,
    tree_data,
    tip_level_pattern: Optional[dict] = None,
    posterior_data: Optional[dict] = None,
    tip_diagnostics: Optional[Dict[str, Any]] = None,
    posterior_threshold: float = 0.7,
    convergence_mode: str = "focal_clade",
    node_index: Optional[Dict[int, Any]] = None,
    build_node_posteriors: bool = False,
    per_site_dist_cache: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    hypothesis: Optional[str] = None,
    walk_cache: Optional[Dict[Any, Any]] = None,
) -> ConvergenceResult:
    """
    Perform complete convergence/disambiguation analysis for a CAAS position.

    ``build_node_posteriors`` (default False): populate the heavy
    ``node_posteriors["per_node"]`` map — the modal AA + full 20-AA distribution
    for *every* tree node at the focal site. This field is NOT in the master CSV,
    the aggregation DB, or the per-gene JSON summary, and the tree plots reload
    node posteriors from the separately-dumped ASR object (see
    ``gene_wrapper.export_posteriors_to_jsonl``), overwriting whatever is here.
    It is therefore redundant, so it is skipped by default; pass True only if an
    in-memory consumer genuinely needs it.

    Args:
        gene: Gene name
        caas_pos: CAAS position information
        tree_data: Tree structure data
        tip_level_pattern: Pre-computed tip-level pattern analysis
        posterior_data: ASR posterior probabilities
        posterior_threshold: Posterior probability threshold for accepting node states

    Returns:
        ConvergenceResult object with complete analysis
    """
    logger.info(
        f"Analyzing convergence for {gene} position {caas_pos.position_one_based}"
    )
    node_posteriors: Dict[str, Any] = {}  # Ensure node_posteriors is always defined

    # Initialize analysis variables
    ancestral = "?"
    derived = "?"
    tip_diagnostics = tip_diagnostics or {}
    state_source = "unknown"
    tip_pattern_comment = caas_pos.caas or ""

    # Extract change/parallelism classification from pre-computed result
    cp_result = {}
    if tip_level_pattern and isinstance(tip_level_pattern, dict):
        cp_result = tip_level_pattern
    elif tip_diagnostics.get("pair_details"):
        cp_result = classify_change_and_parallelism(
            tip_diagnostics["pair_details"],
            convergence_mode=convergence_mode,
            grouping_scheme=getattr(caas_pos, "caap_group", None),
        )

    change_top = cp_result.get("change_top", "no_change")
    change_bottom = cp_result.get("change_bottom", "no_change")
    change_side = cp_result.get("change_side", "none")
    convergence_type = cp_result.get("convergence_type", "no_change")

    # Build per-pair transition status map for annotations
    pair_status_map: Dict[str, Dict[str, str]] = {}
    for summary_entry in tip_diagnostics.get("pair_transition_summary") or []:
        pair_id_raw = summary_entry.get("pair_id")
        pair_id = str(pair_id_raw) if pair_id_raw is not None else None
        transitions = summary_entry.get("transitions") or {}
        statuses: Dict[str, str] = {}
        for side in ("top", "bottom"):
            status = (transitions.get(side) or {}).get("status")
            if status and status != "unknown":
                statuses[side] = status
        if pair_id and statuses:
            pair_status_map[pair_id] = statuses

    # Perform node-level convergence analysis using ASR node mapping
    if posterior_data is None:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position}: missing posterior data."
        )

    node_state_info = None
    node_state_details = None
    node_role_mapping = tip_diagnostics.get("node_mapping")

    if not node_role_mapping:
        # tip_diagnostics doesn't have node_mapping - this happens when ASR is available
        # but tip-level analysis hasn't been run yet. We need to build it here.
        node_role_mapping = {}
        if tree_data and hasattr(tree_data, "root") and tree_data.root:
            node_role_mapping["root"] = tree_data.root.node_id

        # We can't determine focal nodes without tip analysis, so skip ASR node states
        if not node_role_mapping or len(node_role_mapping) < 4:
            logger.debug(
                f"Skipping ASR node states for {gene} position {caas_pos.position}: "
                f"focal node mapping not available (need tip-level analysis)"
            )
            node_role_mapping = None

    try:
        paml_site = caas_pos.position_one_based  # 1-based index for PAML
        node_state_info = extract_node_states_from_node_level(
            posterior_data,
            node_role_mapping,
            paml_site or -1,
            gene,
            posterior_threshold=posterior_threshold,
            tree_node_lookup=getattr(tree_data, "node_mapping", None),
        )

        node_posteriors: Dict[str, Any] = {"roles": {}, "per_node": {}}
        if node_state_info:
            state_source = "asr"
            node_state_details = {
                "root": node_state_info.root,
                "root_prob": node_state_info.root_prob,
                "mrca_contrast": node_state_info.mrca_contrast,
                "mrca_contrast_prob": node_state_info.mrca_contrast_prob,
                "focal_states": node_state_info.focal_states,
                "focal_probs": node_state_info.focal_probs,
            }

            # Debug: Log focal data structure
            logger.debug(
                f"focal_states={node_state_info.focal_states}, "
                f"focal_probs={node_state_info.focal_probs}, "
                f"len(focal_states)={len(node_state_info.focal_states)}"
            )

            # Dynamic focal states extraction (also store as individual keys)
            for idx in range(1, len(node_state_info.focal_states) + 1):
                state = (
                    node_state_info.focal_states[idx - 1]
                    if idx - 1 < len(node_state_info.focal_states)
                    else None
                )
                prob = (
                    node_state_info.focal_probs[idx - 1]
                    if idx - 1 < len(node_state_info.focal_probs)
                    else None
                )
                node_state_details[f"focal_{idx}"] = state
                node_state_details[f"focal_{idx}_prob"] = prob

            # Store root and mrca_contrast
            for role, data in (
                ("root", (node_state_info.root, node_state_info.root_prob)),
                (
                    "mrca_contrast",
                    (node_state_info.mrca_contrast, node_state_info.mrca_contrast_prob),
                ),
            ):
                aa, prob = data
                if aa:
                    node_posteriors.setdefault("roles", {})[role] = {
                        "aa": aa,
                        "prob": prob,
                    }

            # Store dynamic focal nodes
            for idx, (aa, prob) in enumerate(
                zip(node_state_info.focal_states, node_state_info.focal_probs), 1
            ):
                if aa:
                    node_posteriors.setdefault("roles", {})[f"focal_{idx}"] = {
                        "aa": aa,
                        "prob": prob,
                    }
            # Capture per-node annotations for downstream debugging/plots. Redundant
            # in the normal path (not serialized anywhere; plots reload posteriors
            # from the ASR dump), so only built when explicitly requested.
            if build_node_posteriors:
                per_node_states: Dict[int, Dict[str, Any]] = {}
                for node_id, node_sites in posterior_data.items():
                    site_probs = node_sites.get(paml_site)
                    if not site_probs:
                        continue
                    try:
                        modal_aa, modal_prob = max(
                            site_probs.items(), key=lambda x: x[1]
                        )
                    except ValueError:
                        continue
                    per_node_states[int(node_id)] = {
                        "aa": modal_aa,
                        "prob": modal_prob,
                        "distribution": dict(sorted(site_probs.items())),
                    }
                if per_node_states:
                    node_posteriors["per_node"] = per_node_states

    except Exception as e:
        logger.warning(f"Node-level analysis failed: {e}")

    if state_source != "asr" or not node_state_info:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position}; "
            "cannot analyze without posterior-supported nodes."
        )

    # Determine ancestral/derived states using ASR node states
    if not node_state_info or not node_state_info.mrca_contrast:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position}: "
            "missing MRCA contrast state."
        )

    ancestral = node_state_info.mrca_contrast

    asr_descendants: List[str] = [
        state
        for state in node_state_info.focal_states
        if state and state not in {"-", "?", "X"}
    ]
    if not asr_descendants:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position}: "
            "missing focal lineage states."
        )

    trait1_list = caas_pos.trait1_aa or []
    trait0_list = caas_pos.trait0_aa or []

    def _normalize_node_role_mapping(
        mapping: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        if not mapping:
            return {}
        norm = dict(mapping)
        focal_nodes = mapping.get("focal_nodes")
        if isinstance(focal_nodes, (list, tuple)):
            for idx, node_id in enumerate(focal_nodes, 1):
                if isinstance(node_id, int):
                    norm[f"focal_{idx}"] = node_id
        return norm

    derived_states = sorted(
        {
            state
            for state in asr_descendants
            if state and state not in {"-", "?", "X"} and state != ancestral
        }
    )
    if not derived_states:
        derived = ancestral
    elif len(derived_states) == 1:
        derived = derived_states[0]
    else:
        derived = "/".join(derived_states)

    pair_details_list: List[Dict[str, Any]] = tip_diagnostics.get("pair_details") or []

    if trait1_list or trait0_list:
        top_desc = format_amino_display(trait1_list)
        bottom_desc = format_amino_display(trait0_list)
        tip_pattern_comment = (
            f"{caas_pos.caas or ''} (trait1: {top_desc}, trait0: {bottom_desc})"
        )

    node_summary = {
        "root": node_state_info.root if node_state_info else None,
        "mrca_contrast": node_state_info.mrca_contrast if node_state_info else None,
    }

    # Add dynamic focal states to summary
    if node_state_info:
        for idx, state in enumerate(node_state_info.focal_states, 1):
            node_summary[f"focal_{idx}"] = state

    # Conserved-pair logic driven by metadata row
    is_cons_meta = bool(getattr(caas_pos, "is_conserved_meta", False))
    conserved_pair = str(getattr(caas_pos, "conserved_pair", "") or "").strip()

    # ── ASR path score ───────────────────────────────────────────────────────
    # Unified replacement for the legacy binary ASR gate + convergence + parallel
    # scores. Walks each pair's MRCA up to the root using the per-node posteriors
    # captured above, scoring how isolated each tip change is from the deeper
    # background (signed isolation) or how deep the conserved state extends
    # (unsigned conservation). See src/convergence/path_scores.py.
    asr_path_score = 0.0
    independence = 1.0
    mrca_diversity = 0.0
    derived_agreement = 1.0
    conservation_gate = 1.0
    core = 0.0
    pair_path_scores: Dict[int, float] = {}
    pair_path_contaminated: Dict[int, bool] = {}
    conserved_pair_path_scores: Dict[int, float] = {}
    conserved_pair_path_nodes: Dict[int, Any] = {}
    pair_ancestral_aa: Dict[int, Any] = {}
    pair_derived_top_aa: Dict[int, str] = {}
    pair_derived_bot_aa: Dict[int, str] = {}
    try:
        # Single code path with the axes-only perm replay: _position_axes rebuilds
        # per_node_dist from posterior_data exactly as node_posteriors["per_node"]
        # was built above, so the score here stays bit-identical while guaranteeing
        # the perm null scores each position through the same helper.
        path_result = _position_axes(
            caas_pos, tree_data, posterior_data, node_index, pair_details_list,
            per_site_dist_cache=per_site_dist_cache, walk_cache=walk_cache,
        )
        asr_path_score = path_result["asr_path_score"]
        independence = path_result.get("independence", 1.0)
        mrca_diversity = path_result["mrca_diversity"]
        derived_agreement = path_result["derived_agreement"]
        conservation_gate = path_result["conservation_gate"]
        core = path_result.get("core", 0.0)
        pair_path_scores = path_result["pair_scores"]
        pair_path_contaminated = path_result["pair_contaminated"]
        conserved_pair_path_scores = path_result.get("conserved_pair_scores", {}) or {}
        conserved_pair_path_nodes = path_result.get("conserved_pair_nodes", {}) or {}
        pair_ancestral_aa = path_result.get("pair_ancestral", {}) or {}
        pair_derived_top_aa = path_result.get("pair_derived_top", {}) or {}
        pair_derived_bot_aa = path_result.get("pair_derived_bot", {}) or {}
    except Exception as e:  # never let path scoring break disambiguation
        logger.warning(
            f"ASR path scoring failed for {gene}:{caas_pos.position}: {e}"
        )

    return ConvergenceResult(
        gene=gene,
        position=caas_pos.position,
        tag=caas_pos.tag,
        caas=caas_pos.caas,
        position_one_based=caas_pos.position_one_based,
        ancestral=ancestral,
        derived=derived,
        convergence_type=convergence_type,
        trait1_aa=trait1_list,
        trait0_aa=trait0_list,
        tip_pattern_comment=tip_pattern_comment,
        pair_details=tip_diagnostics.get("pair_details"),
        pair_transition_summary=tip_diagnostics.get("pair_transition_summary"),
        node_mapping=tip_diagnostics.get("node_mapping"),
        asr_ancestral_state=node_state_info.mrca_contrast if node_state_info else None,
        asr_descendant_states=asr_descendants if asr_descendants else None,
        node_state_details=node_state_details,
        node_posteriors=node_posteriors if node_posteriors else None,
        root_state=node_state_info.root if node_state_info else None,
        mrca_state=node_state_info.mrca_contrast if node_state_info else None,
        focal_states=(
            {
                f"focal_{idx}": state
                for idx, state in enumerate(node_state_info.focal_states, 1)
            }
            if node_state_info
            else None
        ),
        node_state_summary=node_summary,
        state_source=state_source,
        derived_similarity=None,
        change_top=change_top,
        change_bottom=change_bottom,
        change_side=change_side,
        caap_group=getattr(caas_pos, "caap_group", "US"),
        amino_encoded=getattr(caas_pos, "amino_encoded", ""),
        is_conserved_meta=is_cons_meta,
        conserved_pair=conserved_pair,
        hypothesis=hypothesis,
        asr_path_score=asr_path_score,
        independence=independence,
        mrca_diversity=mrca_diversity,
        derived_agreement=derived_agreement,
        conservation_gate=conservation_gate,
        core=core,
        pair_path_scores=pair_path_scores or None,
        pair_path_contaminated=pair_path_contaminated or None,
        conserved_pair_path_scores=conserved_pair_path_scores or None,
        conserved_pair_path_nodes=conserved_pair_path_nodes or None,
        pair_ancestral_aa=pair_ancestral_aa or None,
        pair_derived_top_aa=pair_derived_top_aa or None,
        pair_derived_bot_aa=pair_derived_bot_aa or None,
        score=None,
        pvalue_boot=getattr(caas_pos, "pvalue_boot", None),
    )


def analyze_gene_disambiguation(
    gene: str,
    alignment_data,
    tree_data,
    caas_positions: List[int],
    caas_entries: Optional[List[CAASPosition]] = None,
    caas_metadata_path: Optional[Path] = None,
    trait_file_path: Optional[Path] = None,
    taxid_mapping: Optional[Dict[str, str]] = None,
    posterior_data: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    posterior_threshold: float = 0.7,
    diagnostics_dir: Optional[Path] = None,
    convergence_mode: str = "focal_state",
    asr_mode: str = "precomputed",
    axes_only: bool = False,
    per_site_dist_cache: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    build_node_posteriors: bool = False,
) -> Tuple[List[ConvergenceResult], Dict[str, Any]]:
    """
    Perform complete convergence/disambiguation analysis for a gene's CAAS positions.

    Args:
        gene: Gene name
        alignment_data: Alignment and lookup data
        tree_data: Tree structure data
        caas_positions: List of CAAS positions to analyze
        caas_metadata_path: Optional path to CAAS metadata
        trait_file_path: Optional path to trait file
        taxid_mapping: Optional species to taxid mapping
        posterior_data: Optional ASR posterior data
        posterior_threshold: Posterior probability threshold for node state extraction
        axes_only: Reduced-kernel mode for the permulation replay. When True, each
            position still builds pair_details + the (cheap) change classification,
            but SKIPS the full per-position scorer (its redundant per-node
            posterior-map rebuild + the 40-field ConvergenceResult the perm null
            discards). It emits a
            lightweight :class:`PositionAxes` per position carrying only
            asr_path_score + change_top/change_bottom. The asr_path_score is scored
            by the same :func:`_position_axes` helper as the full path, so it is
            bit-for-bit identical.

    Returns:
        Tuple of (results list, diagnostics dict). In axes_only mode the list holds
        :class:`PositionAxes` records instead of :class:`ConvergenceResult`.
    """
    logger.info(
        f"Starting convergence disambiguation for {gene} ({len(caas_positions)} positions)"
    )
    logger.debug(
        f"Using posterior threshold {posterior_threshold:.3f} for node state extraction"
    )

    if per_site_dist_cache is None:
        per_site_dist_cache = {}

    # Memoises the labeling-invariant MRCA->root walks inside compute_asr_path_score
    # (see path_scores._changed_side_walk). Shared across every position, scheme and
    # FOP hypothesis of this gene: a candidate pair recurring across hypotheses /
    # cycles is walked once. Per-node posteriors for a site are content-identical on
    # every rebuild, so the cache stays valid for the observed path too. Bounded by
    # (#sites x #schemes x #distinct pair MRCAs x #residues) entries per gene.
    gene_walk_cache: Dict[Any, Any] = {}

    results: List[ConvergenceResult] = []
    diagnostics: Dict[str, Any] = {
        "skipped_positions": 0,
        "skip_reasons": Counter(),
        "tip_dump_file": None,
    }
    # Prepare a streaming JSONL file for tip diagnostics if requested
    tip_file_handle = None
    tip_file_path = None
    if diagnostics_dir:
        tip_dir = Path(diagnostics_dir) / "tip_details"
        tip_dir.mkdir(parents=True, exist_ok=True)
        tip_file_path = tip_dir / f"{gene.lower()}_tip_details.jsonl"
        # We open lazily when the first record is ready so we avoid creating empty files

    # Load row-wise CAAS metadata entries
    if caas_entries is None:
        caas_entries = []
        if caas_metadata_path:
            caas_entries = list_gene_caas_entries(Path(caas_metadata_path), gene)
        if not caas_entries:
            caas_entries = [
                CAASPosition(
                    position=pos,
                    position_one_based=pos + 1,
                    tag=f"POS{pos}",
                    caas="",
                    trait1_aa=[],
                    trait0_aa=[],
                )
                for pos in caas_positions
            ]

    # ── Trait pairs, grouped by contrast ──────────────────────────────────────
    # parse_trait_pairs returns {contrast -> pairs}: one contrast per FOP
    # hypothesis (traitfile_H<n>.tab -> n), or a single contrast for a plain
    # non-FOP traitfile. Each CAAS metadata row belongs to the hypothesis that
    # discovered it (its `trait` field); it is disambiguated against THAT
    # hypothesis's pairs only. Unioning every hypothesis's pairs (the old
    # `flattened_pairs`) polluted every multi-pair axis of the ASR path score —
    # LAC merge points, independence, mrca_diversity — with contrasts from
    # unrelated hypotheses, and discarded the per-hypothesis Dunn independence
    # the FOP harvest enforces.
    trait_pairs_all: Dict[int, List[Tuple[str, str]]] = {}
    if trait_file_path:
        trait_pairs_all = parse_trait_pairs(Path(trait_file_path))

    def _dedup_pairs(pairs: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
        seen: set = set()
        out: List[Tuple[str, str]] = []
        for pair in pairs:
            key = tuple(pair)
            if key not in seen:
                seen.add(key)
                out.append(pair)
        return out

    contrast_pairs_by_key: Dict[int, List[Tuple[str, str]]] = {
        k: _dedup_pairs(v) for k, v in trait_pairs_all.items()
    }
    _multi_contrast = len(contrast_pairs_by_key) > 1
    _single_contrast_key = (
        next(iter(contrast_pairs_by_key)) if len(contrast_pairs_by_key) == 1 else None
    )
    # Safety net only: a multi-contrast run whose metadata row carries no
    # resolvable hypothesis tag falls back to the union (with a warning) rather
    # than silently dropping the position.
    _flattened_fallback = _dedup_pairs(
        [p for pairs in trait_pairs_all.values() for p in pairs]
    )
    any_pairs = bool(_flattened_fallback)

    def _resolve_contrast(entry) -> Tuple[Optional[int], Optional[str]]:
        """(contrast key, hypothesis label) for one CAAS metadata row.

        FOP rows carry `trait` containing "H<n>"; map to contrast <n>. A single
        non-FOP traitfile ignores `trait` and uses its lone contrast, emitting no
        hypothesis label so downstream output is byte-identical to before.
        """
        raw = str(getattr(entry, "trait", "") or "").strip()
        m = re.search(r"H(\d+)", raw)
        if m:
            k = int(m.group(1))
            return (k if k in contrast_pairs_by_key else None), f"H{k}"
        if _single_contrast_key is not None and not _multi_contrast:
            return _single_contrast_key, None
        if _single_contrast_key is not None:
            return _single_contrast_key, (raw or None)
        return None, (raw or None)

    # ── Hoisted per-gene/tree invariants (computed ONCE, not per position) ──────
    # The alignment lookup, the tree node index, and each species-pair's MRCA
    # depend only on the alignment + tree, which are fixed for the gene. Rebuilding
    # them inside the per-position loop (and, for the permulation replay, once per
    # position × per cycle) is pure redundant work. Hoisting/memoising them here
    # speeds up both the observed disambiguation and the perm replay with NO change
    # to results.
    hoisted_seq_by_id = hoisted_seq_by_species = None
    if any_pairs and taxid_mapping:
        hoisted_seq_by_id, hoisted_seq_by_species = build_alignment_lookup(
            alignment_data.alignment, alignment_data.taxid_to_species
        )
    hoisted_node_index = build_node_index(getattr(tree_data, "root", None))
    _mrca_cache: Dict[tuple, Any] = {}

    def _get_mrca_cached(taxa: List[str]):
        if not tree_data or not taxa:
            return None
        key = tuple(sorted(taxa))
        if key not in _mrca_cache:
            _mrca_cache[key] = get_mrca(tree_data.root, list(taxa))
        return _mrca_cache[key]

    for idx, caas_pos in enumerate(caas_entries):
        pos = caas_pos.position
        # Contrast (hypothesis) this row belongs to; its pairs are the ONLY ones
        # this row is scored against.
        _ckey, _hyp_label = _resolve_contrast(caas_pos)
        if _ckey is not None:
            contrast_pairs = contrast_pairs_by_key.get(_ckey, [])
        elif _hyp_label and _hyp_label.startswith("H"):
            # Row names a hypothesis whose traitfile is absent → cannot score it.
            contrast_pairs = []
            logger.warning(
                f"{gene} pos {pos}: metadata hypothesis {_hyp_label} has no "
                "matching traitfile; position skipped"
            )
        else:
            if _multi_contrast:
                logger.warning(
                    f"{gene} pos {pos}: multi-hypothesis run but row has no "
                    "resolvable hypothesis tag; falling back to the pooled union"
                )
            contrast_pairs = _flattened_fallback
        try:
            # Skip if no amino acid conversion data
            if not caas_pos.caas:
                logger.debug(f"Skipping position {pos} - no amino acid conversion data")
                diagnostics["skip_reasons"]["no_caasersion"] += 1
                diagnostics["skipped_positions"] += 1
                continue

            if posterior_data is None:
                diagnostics["skip_reasons"]["no_asr"] += 1
                diagnostics["skipped_positions"] += 1
                if asr_mode != "skip":
                    logger.debug(
                        f"Skipping position {pos} - ASR posterior data missing (asr_mode={asr_mode})"
                    )
                continue

            # Perform tip-level convergence analysis across ALL pairs from trait file
            tip_level_pattern = None
            tip_diagnostics: Dict[str, Any] = {}
            try:
                if contrast_pairs and taxid_mapping:
                    seq_by_id, seq_by_species = hoisted_seq_by_id, hoisted_seq_by_species

                    pair_details = []
                    all_taxa: List[str] = []

                    def _modal_state(
                        node_id: Optional[int],
                    ) -> Tuple[Optional[str], Optional[float]]:
                        if node_id is None or posterior_data is None:
                            return None, None
                        site = caas_pos.position_one_based
                        node_dict = posterior_data.get(node_id, {})
                        site_post = node_dict.get(site, {}) if site is not None else {}
                        if not site_post:
                            return None, None
                        aa, prob = max(site_post.items(), key=lambda x: x[1])
                        return aa, prob

                    for pair_idx, (high_species, low_species) in enumerate(
                        contrast_pairs, 1
                    ):
                        top_taxid = str(taxid_mapping.get(high_species, high_species))
                        bottom_taxid = str(taxid_mapping.get(low_species, low_species))

                        all_taxa.extend([top_taxid, bottom_taxid])

                        mrca_node = _get_mrca_cached([top_taxid, bottom_taxid])
                        mrca_state, mrca_prob = _modal_state(
                            mrca_node.node_id if mrca_node else None
                        )

                        top_tip_records = collect_tip_residues(
                            [top_taxid],
                            [high_species],
                            caas_pos.position,
                            seq_by_id,
                            seq_by_species,
                            alignment_data.taxid_to_species,
                        )
                        bottom_tip_records = collect_tip_residues(
                            [bottom_taxid],
                            [low_species],
                            caas_pos.position,
                            seq_by_id,
                            seq_by_species,
                            alignment_data.taxid_to_species,
                        )

                        top_tip = extract_tip_residue(top_tip_records)
                        bottom_tip = extract_tip_residue(bottom_tip_records)

                        pair_details.append(
                            {
                                "pair_id": pair_idx,
                                "node_id": mrca_node.node_id if mrca_node else None,
                                "focal_state": mrca_state,
                                "focal_prob": mrca_prob,
                                "mrca_modal_aa": mrca_state,
                                "top_taxa": [top_taxid],
                                "bottom_taxa": [bottom_taxid],
                                "top_species": [high_species],
                                "bottom_species": [low_species],
                                "top_tip_mode": top_tip,
                                "bottom_tip_mode": bottom_tip,
                                "top_tip_residue": top_tip,
                                "bottom_tip_residue": bottom_tip,
                                "top_tip_residues": top_tip_records,
                                "bottom_tip_residues": bottom_tip_records,
                            }
                        )

                    mrca_node = _get_mrca_cached(all_taxa)
                    # Populate mrca_contrast in each pair using the global contrast MRCA
                    # state (MRCA of all taxa). classify_change_and_parallelism uses this
                    # field as the ancestor when convergence_mode="mrca".
                    global_mrca_state, _ = _modal_state(
                        mrca_node.node_id if mrca_node else None
                    )
                    if global_mrca_state is None:
                        logger.warning(
                            f"Global MRCA state unavailable for {gene} pos {pos}; "
                            "mrca convergence_mode will fall back to focal_state per pair"
                        )
                    for p in pair_details:
                        p["mrca_contrast"] = global_mrca_state if global_mrca_state is not None else p.get("focal_state")
                    node_mapping = {
                        "root": (
                            tree_data.root.node_id
                            if tree_data and tree_data.root
                            else None
                        ),
                        "mrca_contrast": mrca_node.node_id if mrca_node else None,
                        "focal_nodes": [p.get("node_id") for p in pair_details],
                    }

                    tip_level_pattern = classify_change_and_parallelism(
                        pair_details,
                        convergence_mode=convergence_mode,
                        grouping_scheme=getattr(caas_pos, "caap_group", None),
                    )
                    tip_diagnostics["pair_details"] = pair_details
                    tip_diagnostics["node_mapping"] = node_mapping
            except Exception as e:
                logger.warning(f"Tip-level analysis failed for position {pos}: {e}")

            if diagnostics_dir and tip_diagnostics.get("pair_details"):
                # Stream tip diagnostic record to JSONL, don't accumulate in memory
                import json

                record = {
                    "gene": gene,
                    "position": caas_pos.position,
                    "position_one_based": caas_pos.position_one_based,
                    "tag": caas_pos.tag,
                    "pair_details": tip_diagnostics.get("pair_details"),
                }
                try:
                    if tip_file_handle is None:
                        if tip_file_path is None:
                            raise ValueError("tip_file_path is None")
                        tip_file_handle = open(tip_file_path, "a", encoding="utf-8")
                        diagnostics["tip_dump_file"] = str(tip_file_path)
                    tip_file_handle.write(json.dumps(record) + "\n")
                    tip_file_handle.flush()
                except Exception as e:
                    logger.warning(
                        f"Failed to write tip diagnostic for {gene} pos {pos}: {e}"
                    )

            # If no valid trait pairs overlap the alignment/taxid set, we cannot
            # build the focal node mapping needed for ASR node-level analysis.
            if not tip_diagnostics.get("pair_details"):
                diagnostics["skip_reasons"]["no_valid_pairs"] += 1
                diagnostics["skipped_positions"] += 1
                logger.warning(
                    f"Skipping position {pos} - no valid trait pairs overlap alignment/taxid mapping"
                )
                continue

            # ── Axes-only reduced kernel (permulation replay) ──────────────────
            # The perm null keeps only asr_path_score + the change partition, so we
            # skip the full per-position scorer — whose cost is the redundant per-node
            # posterior-map rebuild (node_posteriors["per_node"]) plus the 40-field
            # ConvergenceResult assembly the null discards — and call the shared
            # _position_axes helper directly (which builds per_node_dist once, cached).
            # change_top/change_bottom come straight from the (cheap) classification
            # computed just above, so they are the EXACT categories the full path
            # would have emitted — not a placeholder.
            if axes_only:
                cp = tip_level_pattern or {}
                try:
                    path_result = _position_axes(
                        caas_pos,
                        tree_data,
                        posterior_data,
                        hoisted_node_index,
                        tip_diagnostics.get("pair_details"),
                        per_site_dist_cache=per_site_dist_cache,
                        walk_cache=gene_walk_cache,
                    )
                    axes_score = path_result.get("asr_path_score", 0.0)
                    axes_pair_scores = path_result.get("pair_scores", None)
                    axes_cons_scores = path_result.get("conserved_pair_scores", None) or None
                    axes_cons_nodes = path_result.get("conserved_pair_nodes", None) or None
                    axes_anc = path_result.get("pair_ancestral", None) or None
                    axes_der_top = path_result.get("pair_derived_top", None) or None
                    axes_der_bot = path_result.get("pair_derived_bot", None) or None
                    axes_extra = {
                        k: path_result.get(k)
                        for k in ("independence", "mrca_diversity",
                                  "derived_agreement", "conservation_gate", "core")
                    }
                except Exception as e:  # never let path scoring break the replay
                    logger.warning(
                        f"ASR path scoring failed for {gene}:{caas_pos.position}: {e}"
                    )
                    axes_score = 0.0
                    axes_pair_scores = None
                    axes_cons_scores = None
                    axes_cons_nodes = None
                    axes_anc = axes_der_top = axes_der_bot = None
                    axes_extra = {}
                results.append(
                    PositionAxes(
                        position=caas_pos.position,
                        caap_group=getattr(caas_pos, "caap_group", "US"),
                        asr_path_score=axes_score,
                        change_top=cp.get("change_top", "no_change"),
                        change_bottom=cp.get("change_bottom", "no_change"),
                        hypothesis=_hyp_label,
                        pair_scores=axes_pair_scores,
                        independence=axes_extra.get("independence"),
                        mrca_diversity=axes_extra.get("mrca_diversity"),
                        derived_agreement=axes_extra.get("derived_agreement"),
                        conservation_gate=axes_extra.get("conservation_gate"),
                        core=axes_extra.get("core"),
                        conserved_pair_scores=axes_cons_scores,
                        conserved_pair_nodes=axes_cons_nodes,
                        pair_ancestral=axes_anc,
                        pair_derived_top=axes_der_top,
                        pair_derived_bot=axes_der_bot,
                    )
                )
                logger.debug(
                    f"✓ Axes-only position {pos}: asr_path_score={axes_score}"
                )
                continue

            # Perform convergence/disambiguation analysis
            result = analyze_caas_position_disambiguation(
                gene,
                caas_pos,
                tree_data,
                tip_level_pattern,
                posterior_data,
                tip_diagnostics,
                posterior_threshold=posterior_threshold,
                convergence_mode=convergence_mode,
                node_index=hoisted_node_index,
                build_node_posteriors=build_node_posteriors,
                per_site_dist_cache=per_site_dist_cache,
                hypothesis=_hyp_label,
                walk_cache=gene_walk_cache,
            )

            results.append(result)
            logger.info(
                f"✓ Analyzed position {pos}: {result.convergence_type} pattern, "
                f"{result.ancestral}→{result.derived}"
            )

        except Exception as e:
            logger.exception(f"Failed to analyze position {pos}: {e}")
            diagnostics["skipped_positions"] += 1
            diagnostics["skip_reasons"][str(e).split(":")[0]] += 1
            continue

    if tip_file_handle is not None:
        try:
            tip_file_handle.close()
        except Exception:
            pass
        logger.info(f"Tip details written to {diagnostics.get('tip_dump_file')}")

    logger.info(
        f"✓ Completed convergence disambiguation: {len(results)}/{len(caas_entries)} metadata rows"
    )
    diagnostics["skip_reasons"] = dict(diagnostics["skip_reasons"])
    return results, diagnostics


# Backward-compatible aliases
analyze_caas_position_biochemistry = analyze_caas_position_disambiguation
analyze_gene_biochemistry = analyze_gene_disambiguation
