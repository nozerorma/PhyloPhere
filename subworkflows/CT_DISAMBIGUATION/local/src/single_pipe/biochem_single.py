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

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import logging
from collections import Counter

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
    classify_tip_level_pattern
)

from src.biochem.state_inference import compute_change_side
from src.asr.tree_parser import get_mrca
from src.data.models import CAASPosition, BiochemResults
from src.data.loaders import list_gene_caas_entries, parse_trait_pairs

logger = logging.getLogger(__name__)


def analyze_caas_position_biochemistry(
    gene: str,
    caas_pos: CAASPosition,
    tree_data,
    tip_level_pattern: Optional[dict] = None,
    posterior_data: Optional[dict] = None,
    tip_diagnostics: Optional[Dict[str, Any]] = None,
    posterior_threshold: float = 0.7,
    convergence_mode: str = 'focal_clade'
) -> BiochemResults:
    """
    Perform complete biochemical analysis for a CAAS position.

    Args:
        gene: Gene name
        caas_pos: CAAS position information
        tree_data: Tree structure data
        tip_level_pattern: Pre-computed tip-level pattern analysis
        posterior_data: ASR posterior probabilities
        posterior_threshold: Posterior probability threshold for accepting node states

    Returns:
        BiochemResults object with complete analysis
    """
    logger.info(f"Analyzing convergence for {gene} position {caas_pos.position_one_based}")
    node_posteriors: Dict[str, Any] = {}  # Ensure node_posteriors is always defined

    # Initialize analysis variables
    ancestral = '?';
    derived = '?';
    pattern_type = 'unknown'
    convergence_desc = 'No analysis available'
    tip_diagnostics = tip_diagnostics or {}
    state_source = 'unknown'
    tip_pattern_comment = caas_pos.caas or ''
    # Determine pattern type from tip-level analysis
    if tip_level_pattern and isinstance(tip_level_pattern, dict):
        pattern_type = tip_level_pattern.get('pattern', 'unknown')
        convergence_desc = tip_level_pattern.get('description', 'Unknown pattern')
    elif tip_diagnostics.get('pair_details'):
        inferred_pattern = classify_tip_level_pattern(
            tip_diagnostics['pair_details'],
            convergence_mode=convergence_mode,
            grouping_scheme=getattr(caas_pos, 'caap_group', None)
        )
        pattern_type = inferred_pattern.get('pattern', 'unknown')
        convergence_desc = inferred_pattern.get('description', 'Unknown pattern')

    # Track which trait groups changed (top/bottom) and how
    top_change_type = 'none'
    bottom_change_type = 'none'
    if tip_level_pattern:
        if tip_level_pattern.get('top_convergent'):
            top_change_type = 'convergent'
        elif tip_level_pattern.get('top_divergent'):
            top_change_type = 'divergent'
        elif tip_level_pattern.get('top_insufficient'):
            top_change_type = 'insufficient'

        if tip_level_pattern.get('bottom_convergent'):
            bottom_change_type = 'convergent'
        elif tip_level_pattern.get('bottom_divergent'):
            bottom_change_type = 'divergent'
        elif tip_level_pattern.get('bottom_insufficient'):
            bottom_change_type = 'insufficient'
            
    # Build per-pair transition status map for annotations
    pair_status_map: Dict[str, Dict[str, str]] = {}
    for summary_entry in tip_diagnostics.get('pair_transition_summary') or []:
        pair_id_raw = summary_entry.get('pair_id')
        pair_id = str(pair_id_raw) if pair_id_raw is not None else None
        transitions = summary_entry.get('transitions') or {}
        statuses: Dict[str, str] = {}
        for side in ('top', 'bottom'):
            status = (transitions.get(side) or {}).get('status')
            if status and status != 'unknown':
                statuses[side] = status
        if pair_id and statuses:
            pair_status_map[pair_id] = statuses

    # Perform node-level convergence analysis using ASR node mapping
    if posterior_data is None:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position_zero_based}: missing posterior data."
        )

    node_state_info = None
    node_state_details = None
    node_role_mapping = tip_diagnostics.get('node_mapping')

    if not node_role_mapping:
        # tip_diagnostics doesn't have node_mapping - this happens when ASR is available
        # but tip-level analysis hasn't been run yet. We need to build it here.
        node_role_mapping = {}
        if tree_data and hasattr(tree_data, 'root') and tree_data.root:
            node_role_mapping['root'] = tree_data.root.node_id

        # We can't determine focal nodes without tip analysis, so skip ASR node states
        if not node_role_mapping or len(node_role_mapping) < 4:
            logger.debug(
                f"Skipping ASR node states for {gene} position {caas_pos.position_zero_based}: "
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
            tree_node_lookup=getattr(tree_data, 'node_mapping', None)
        )

        node_posteriors: Dict[str, Any] = {'roles': {}, 'per_node': {}}
        if node_state_info:
            state_source = 'asr'
            node_state_details = {
                'root': node_state_info.root,
                'root_prob': node_state_info.root_prob,
                'mrca_contrast': node_state_info.mrca_contrast,
                'mrca_contrast_prob': node_state_info.mrca_contrast_prob,
                'focal_states': node_state_info.focal_states,
                'focal_probs': node_state_info.focal_probs,
            }

            # Debug: Log focal data structure
            logger.debug(
                f"focal_states={node_state_info.focal_states}, "
                f"focal_probs={node_state_info.focal_probs}, "
                f"len(focal_states)={len(node_state_info.focal_states)}"
            )

            # Dynamic focal states extraction (also store as individual keys)
            for idx in range(1, len(node_state_info.focal_states) + 1):
                state = node_state_info.focal_states[idx - 1] if idx - 1 < len(node_state_info.focal_states) else None
                prob = node_state_info.focal_probs[idx - 1] if idx - 1 < len(node_state_info.focal_probs) else None
                node_state_details[f'focal_{idx}'] = state
                node_state_details[f'focal_{idx}_prob'] = prob

            # Store root and mrca_contrast
            for role, data in (
                ('root', (node_state_info.root, node_state_info.root_prob)),
                ('mrca_contrast', (node_state_info.mrca_contrast, node_state_info.mrca_contrast_prob)),
            ):
                aa, prob = data
                if aa:
                    node_posteriors.setdefault('roles', {})[role] = {'aa': aa, 'prob': prob}

            # Store dynamic focal nodes
            for idx, (aa, prob) in enumerate(zip(node_state_info.focal_states, node_state_info.focal_probs), 1):
                if aa:
                    node_posteriors.setdefault('roles', {})[f'focal_{idx}'] = {'aa': aa, 'prob': prob}
            if node_state_info.low_confidence_nodes:
                logger.warning(
                    f"Low-confidence node posteriors for {gene}:{caas_pos.position_zero_based} "
                    f"at nodes {node_state_info.low_confidence_nodes}"
                )
                node_state_details['low_confidence_nodes'] = node_state_info.low_confidence_nodes
                node_posteriors['low_confidence_nodes'] = node_state_info.low_confidence_nodes

            # Capture per-node annotations for downstream debugging/plots
            per_node_states: Dict[int, Dict[str, Any]] = {}
            for node_id, node_sites in posterior_data.items():
                site_probs = node_sites.get(paml_site)
                if not site_probs:
                    continue
                try:
                    modal_aa, modal_prob = max(site_probs.items(), key=lambda x: x[1])
                except ValueError:
                    continue
                per_node_states[int(node_id)] = {
                    'aa': modal_aa,
                    'prob': modal_prob,
                    'distribution': dict(sorted(site_probs.items()))
                }
            if per_node_states:
                node_posteriors['per_node'] = per_node_states

    except Exception as e:
        logger.warning(f"Node-level analysis failed: {e}")

    if state_source != 'asr' or not node_state_info:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position_zero_based}; "
            "cannot analyze without posterior-supported nodes."
        )

    # Determine ancestral/derived states using ASR node states
    if not node_state_info or not node_state_info.mrca_contrast:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position_zero_based}: "
            "missing MRCA contrast state."
        )

    ancestral = node_state_info.mrca_contrast

    asr_descendants: List[str] = [
        state for state in node_state_info.focal_states
        if state and state not in {'-', '?', 'X'}
    ]
    if not asr_descendants:
        raise ValueError(
            f"ASR node states unavailable for {gene} position {caas_pos.position_zero_based}: "
            "missing focal lineage states."
        )

    trait1_list = caas_pos.trait1_aa or []
    trait0_list = caas_pos.trait0_aa or []

    def _normalize_node_role_mapping(mapping: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        if not mapping:
            return {}
        norm = dict(mapping)
        focal_nodes = mapping.get('focal_nodes')
        if isinstance(focal_nodes, (list, tuple)):
            for idx, node_id in enumerate(focal_nodes, 1):
                if isinstance(node_id, int):
                    norm[f'focal_{idx}'] = node_id
        return norm

    derived_states = sorted({
        state for state in asr_descendants
        if state and state not in {'-', '?', 'X'} and state != ancestral
    })
    if not derived_states:
        derived = ancestral
    elif len(derived_states) == 1:
        derived = derived_states[0]
    else:
        derived = '/'.join(derived_states)

    pair_details_list: List[Dict[str, Any]] = tip_diagnostics.get('pair_details') or []

    # Compute change side using extracted function (now pattern-aware)
    change_side = compute_change_side(top_change_type, bottom_change_type, pattern_type)

    if trait1_list or trait0_list:
        top_desc = format_amino_display(trait1_list)
        bottom_desc = format_amino_display(trait0_list)
        tip_pattern_comment = f"{caas_pos.caas or ''} (trait1: {top_desc}, trait0: {bottom_desc})"

    node_summary = {
        'root': node_state_info.root if node_state_info else None,
        'mrca_contrast': node_state_info.mrca_contrast if node_state_info else None,
    }

    # Add dynamic focal states to summary
    if node_state_info:
        for idx, state in enumerate(node_state_info.focal_states, 1):
            node_summary[f'focal_{idx}'] = state

    # Conserved-pair logic driven by metadata row
    is_cons_meta = bool(getattr(caas_pos, 'is_conserved_meta', False))
    conserved_pair = str(getattr(caas_pos, 'conserved_pair', '') or '').strip()
    asr_cons = False
    if is_cons_meta and conserved_pair:
        try:
            pair_idx = int(conserved_pair)
        except Exception:
            pair_idx = None
        if pair_idx and 1 <= pair_idx <= len(pair_details_list):
            pair = pair_details_list[pair_idx - 1] or {}
            mrca = pair.get('focal_state') or pair.get('mrca_modal_aa')
            top_tip = pair.get('top_tip_mode') or pair.get('top_tip_residue')
            bottom_tip = pair.get('bottom_tip_mode') or pair.get('bottom_tip_residue')
            asr_cons = bool(mrca and top_tip and bottom_tip and mrca == top_tip == bottom_tip)

    return BiochemResults(
        gene=gene,
        position=caas_pos.position,
        tag=caas_pos.tag,
        caas=caas_pos.caas,
        is_significant=caas_pos.is_significant,
        position_zero_based=caas_pos.position_zero_based,
        position_one_based=caas_pos.position_one_based,
        ancestral=ancestral,
        derived=derived,
        pattern_type=pattern_type,
        convergence_description=convergence_desc,
        convergence_mode=convergence_mode,
        trait1_aa=trait1_list,
        trait0_aa=trait0_list,
        tip_pattern_comment=tip_pattern_comment,
        pair_details=tip_diagnostics.get('pair_details'),
        pair_transition_summary=tip_diagnostics.get('pair_transition_summary'),
        node_mapping=tip_diagnostics.get('node_mapping'),
        asr_ancestral_state=node_state_info.mrca_contrast if node_state_info else None,
        asr_descendant_states=asr_descendants if asr_descendants else None,
        node_state_details=node_state_details,
        node_posteriors=node_posteriors if node_posteriors else None,
        root_state=node_state_info.root if node_state_info else None,
        mrca_state=node_state_info.mrca_contrast if node_state_info else None,
        focal_states={
            f'focal_{idx}': state
            for idx, state in enumerate(node_state_info.focal_states, 1)
        } if node_state_info else None,
        node_state_summary=node_summary,
        low_confidence_nodes=node_state_info.low_confidence_nodes if node_state_info else None,
        state_source=state_source,
        derived_similarity=None,
        top_change_type=top_change_type,
        bottom_change_type=bottom_change_type,
        change_side=change_side,
        caap_group=getattr(caas_pos, 'caap_group', 'US'),
        amino_encoded=getattr(caas_pos, 'amino_encoded', ''),
        is_conserved_meta=is_cons_meta,
        conserved_pair=conserved_pair,
        sig_hyp=getattr(caas_pos, 'sig_hyp', None),
        sig_perm=getattr(caas_pos, 'sig_perm', None),
        sig_both=getattr(caas_pos, 'sig_both', None),
        asr_is_conserved=asr_cons,
        score=None,
        caas_pvalue=caas_pos.pvalue,
        pvalue_boot=getattr(caas_pos, 'pvalue_boot', None)
    )


def analyze_gene_biochemistry(
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
    convergence_mode: str = 'focal_clade',
    asr_mode: str = 'precomputed',
    include_non_significant: bool = False,
) -> Tuple[List[BiochemResults], Dict[str, Any]]:
    """
    Perform complete biochemical analysis for a gene's CAAS positions.

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

    Returns:
        Tuple of (results list, diagnostics dict)
    """
    logger.info(f"Starting biochemical analysis for {gene} ({len(caas_positions)} positions)")
    logger.debug(f"Using posterior threshold {posterior_threshold:.3f} for node state extraction")

    results: List[BiochemResults] = []
    diagnostics: Dict[str, Any] = {
        'skipped_positions': 0,
        'skip_reasons': Counter(),
        'low_confidence_positions': 0,
        'tip_dump_file': None
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
                    position_zero_based=pos,
                    position_one_based=pos + 1,
                    tag=f'POS{pos}',
                    caas='',
                    trait1_aa=[],
                    trait0_aa=[]
                )
                for pos in caas_positions
            ]

    # Load all trait pairs once for uniform processing
    trait_pairs_all: List[Tuple[str, str]] = []
    flattened_pairs: List[Tuple[str, str]] = []
    if trait_file_path:
        trait_pairs_all = parse_trait_pairs(Path(trait_file_path))
        seen_pairs = set()
        for pair in trait_pairs_all:
            pair_tuple = tuple(pair)
            if pair_tuple not in seen_pairs:
                flattened_pairs.append(pair)
                seen_pairs.add(pair_tuple)

    for idx, caas_pos in enumerate(caas_entries):
        pos = caas_pos.position_zero_based
        try:
            # Skip non-significant positions unless explicitly included
            if not include_non_significant and not getattr(caas_pos, 'is_significant', False):
                diagnostics['skip_reasons']['non_significant'] += 1
                diagnostics['skipped_positions'] += 1
                logger.debug(f"Skipping position {pos} (non-significant and include_non_significant=False)")
                continue
            # Skip if no amino acid conversion data
            if not caas_pos.caas:
                logger.debug(f"Skipping position {pos} - no amino acid conversion data")
                diagnostics['skip_reasons']['no_caasersion'] += 1
                diagnostics['skipped_positions'] += 1
                continue

            if posterior_data is None:
                diagnostics['skip_reasons']['no_asr'] += 1
                diagnostics['skipped_positions'] += 1
                if asr_mode != 'skip':
                    logger.debug(
                        f"Skipping position {pos} - ASR posterior data missing (asr_mode={asr_mode})"
                    )
                continue

            # Perform tip-level convergence analysis across ALL pairs from trait file
            tip_level_pattern = None
            tip_diagnostics: Dict[str, Any] = {}
            try:
                if flattened_pairs and taxid_mapping:
                    seq_by_id, seq_by_species = build_alignment_lookup(
                        alignment_data.alignment,
                        alignment_data.taxid_to_species
                    )

                    pair_details = []
                    all_taxa: List[str] = []

                    def _modal_state(node_id: Optional[int]) -> Tuple[Optional[str], Optional[float]]:
                        if node_id is None or posterior_data is None:
                            return None, None
                        site = caas_pos.position_one_based
                        node_dict = posterior_data.get(node_id, {})
                        site_post = node_dict.get(site, {}) if site is not None else {}
                        if not site_post:
                            return None, None
                        aa, prob = max(site_post.items(), key=lambda x: x[1])
                        return aa, prob

                    for pair_idx, (high_species, low_species) in enumerate(flattened_pairs, 1):
                        top_taxid = str(taxid_mapping.get(high_species, high_species))
                        bottom_taxid = str(taxid_mapping.get(low_species, low_species))

                        all_taxa.extend([top_taxid, bottom_taxid])

                        mrca_node = get_mrca(tree_data.root, [top_taxid, bottom_taxid]) if tree_data else None
                        mrca_state, mrca_prob = _modal_state(mrca_node.node_id if mrca_node else None)

                        top_tip_records = collect_tip_residues(
                            [top_taxid], [high_species],
                            caas_pos.position_zero_based,
                            seq_by_id, seq_by_species,
                            alignment_data.taxid_to_species
                        )
                        bottom_tip_records = collect_tip_residues(
                            [bottom_taxid], [low_species],
                            caas_pos.position_zero_based,
                            seq_by_id, seq_by_species,
                            alignment_data.taxid_to_species
                        )

                        top_tip = extract_tip_residue(top_tip_records)
                        bottom_tip = extract_tip_residue(bottom_tip_records)

                        pair_details.append({
                            'pair_id': pair_idx,
                            'node_id': mrca_node.node_id if mrca_node else None,
                            'focal_state': mrca_state,
                            'focal_prob': mrca_prob,
                            'mrca_modal_aa': mrca_state,
                            'top_taxa': [top_taxid],
                            'bottom_taxa': [bottom_taxid],
                            'top_species': [high_species],
                            'bottom_species': [low_species],
                            'top_tip_mode': top_tip,
                            'bottom_tip_mode': bottom_tip,
                            'top_tip_residue': top_tip,
                            'bottom_tip_residue': bottom_tip,
                            'top_tip_residues': top_tip_records,
                            'bottom_tip_residues': bottom_tip_records
                        })

                    mrca_node = get_mrca(tree_data.root, all_taxa) if tree_data and all_taxa else None
                    node_mapping = {
                        'root': tree_data.root.node_id if tree_data and tree_data.root else None,
                        'mrca_contrast': mrca_node.node_id if mrca_node else None,
                        'focal_nodes': [p.get('node_id') for p in pair_details]
                    }

                    tip_level_pattern = classify_tip_level_pattern(
                        pair_details,
                        convergence_mode=convergence_mode,
                        grouping_scheme=getattr(caas_pos, 'caap_group', None)
                    )
                    tip_diagnostics['pair_details'] = pair_details
                    tip_diagnostics['node_mapping'] = node_mapping
            except Exception as e:
                logger.warning(f"Tip-level analysis failed for position {pos}: {e}")

            if diagnostics_dir and tip_diagnostics.get('pair_details'):
                # Stream tip diagnostic record to JSONL, don't accumulate in memory
                import json
                record = {
                    'gene': gene,
                    'position': caas_pos.position,
                    'position_one_based': caas_pos.position_one_based,
                    'tag': caas_pos.tag,
                    'pair_details': tip_diagnostics.get('pair_details'),
                }
                try:
                    if tip_file_handle is None:
                        if tip_file_path is None:
                            raise ValueError("tip_file_path is None")
                        tip_file_handle = open(tip_file_path, 'a', encoding='utf-8')
                        diagnostics['tip_dump_file'] = str(tip_file_path)
                    tip_file_handle.write(json.dumps(record) + '\n')
                    tip_file_handle.flush()
                except Exception as e:
                    logger.warning(f"Failed to write tip diagnostic for {gene} pos {pos}: {e}")

            # If no valid trait pairs overlap the alignment/taxid set, we cannot
            # build the focal node mapping needed for ASR node-level analysis.
            if not tip_diagnostics.get('pair_details'):
                diagnostics['skip_reasons']['no_valid_pairs'] += 1
                diagnostics['skipped_positions'] += 1
                logger.warning(
                    f"Skipping position {pos} - no valid trait pairs overlap alignment/taxid mapping"
                )
                continue

            # Perform biochemical analysis
            result = analyze_caas_position_biochemistry(
                gene,
                caas_pos,
                tree_data,
                tip_level_pattern,
                posterior_data,
                tip_diagnostics,
                posterior_threshold=posterior_threshold,
                convergence_mode=convergence_mode
            )

            results.append(result)
            if getattr(result, 'low_confidence_nodes', None):
                diagnostics['low_confidence_positions'] += 1
            logger.info(f"✓ Analyzed position {pos}: {result.pattern_type} pattern, "
                       f"{result.ancestral}→{result.derived}")

        except Exception as e:
            logger.exception(f"Failed to analyze position {pos}: {e}")
            diagnostics['skipped_positions'] += 1
            diagnostics['skip_reasons'][str(e).split(':')[0]] += 1
            continue

    if tip_file_handle is not None:
        try:
            tip_file_handle.close()
        except Exception:
            pass
        logger.info(f"Tip details written to {diagnostics.get('tip_dump_file')}")

    logger.info(f"✓ Completed biochemical analysis: {len(results)}/{len(caas_entries)} metadata rows")
    diagnostics['skip_reasons'] = dict(diagnostics['skip_reasons'])
    return results, diagnostics
