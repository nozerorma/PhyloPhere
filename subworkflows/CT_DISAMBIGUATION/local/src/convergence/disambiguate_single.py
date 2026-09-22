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

import dataclasses
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import logging
from collections import Counter, defaultdict, namedtuple

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
)
from src.asr.tree_parser import get_mrca, build_name_taxid_index
from src.data.models import CAASPosition, ConvergenceResult
from src.data.loaders import list_gene_caas_entries, parse_trait_pairs
from src.biochem.grouping import get_grouping_scheme
from src.convergence.path_scores import build_node_index, compute_domain_scores
from src.convergence.fop_pool import pool_domains
from src.convergence.support_fmt import fmt_support

logger = logging.getLogger(__name__)


# Lightweight per-position record emitted by the axes-only (permulation) path.
# The permulation null keeps ONLY these five fields (Gene/cycle are added by the
# worker), so there is no reason to assemble the full 40-field ConvergenceResult.
# It exposes the same attribute names the perm worker reads via getattr, so it is
# a drop-in for a ConvergenceResult on that path.
PositionAxes = namedtuple(
    "PositionAxes",
    ["position", "caap_group", "asr_path_score",
     "side", "hypothesis", "domain_scores", "derived_agreement",
     "sides"],
)
# All fields after `side` are optional. `sides` is the raw
# ``compute_domain_scores`` return ``{"top", "bottom", "domain_meta"}`` for the
# FOP null's treeless pooler (``fop_pool.pool_domains``); populated in the
# axes-only replay, ``None`` on the full observed path.
PositionAxes.__new__.__defaults__ = ("none", None, None, None, None)


def _derive_convergent_call(all_rows: List["ConvergenceResult"], attr: str) -> str:
    """Pool one side's per-hypothesis ``caas``/``amino_encoded`` fg/bg strings into a
    single "<derived>/<ancestral>" call: exclude each pair index where fg==bg
    (conserved -- caap_id.py's own conserved-pair check is this same per-index
    group-equality test, so the two are equivalent; no metadata lookup needed),
    then union the surviving divergent residues per side across every pooled
    hypothesis.

    Replaces the old arbitrary-first-hypothesis passthrough for these two
    columns: the union of divergent residues, across every hypothesis in the
    pool, is a genuine summary rather than one hypothesis's raw pattern.
    """
    fg_chars: set = set()
    bg_chars: set = set()
    for r in all_rows or []:
        val = str(getattr(r, attr, "") or "")
        parts = val.split("/")
        if len(parts) != 2:
            continue
        fg, bg = parts
        if len(fg) != len(bg):
            continue
        for f, b in zip(fg, bg):
            if f == b:
                continue
            fg_chars.add(f)
            bg_chars.add(b)
    if not fg_chars and not bg_chars:
        return ""
    return "".join(sorted(fg_chars)) + "/" + "".join(sorted(bg_chars))


def _pooled_pair_lca(hyp_rows: List[Dict[str, Any]], side: str) -> Optional[List[Tuple[Any, Any, int, float]]]:
    """Union across pooled hypotheses of one side's ``(a, b, lca, contrib)``
    triples, dedup'd by ``(a, b, lca)`` (average ``contrib`` across duplicates —
    the same domain pair can recur across hypotheses with a slightly different
    ``contrib`` if the domains' ``mrca_id`` differs by hypothesis)."""
    acc: Dict[Tuple[Any, Any, int], List[float]] = defaultdict(list)
    for hr in hyp_rows:
        for a, b, lca, contrib in ((hr.get("sides") or {}).get(side) or {}).get("pair_lca", []) or []:
            if lca is None:
                continue
            acc[(a, b, lca)].append(contrib)
    if not acc:
        return None
    return [(a, b, lca, sum(cs) / len(cs)) for (a, b, lca), cs in acc.items()]


def _emit_pooled_side_rows(
    base: ConvergenceResult,
    hyp_rows: List[Dict[str, Any]],
    hyp_pairs_pss: Optional[Dict[Tuple[str, Any], float]] = None,
    all_rows: Optional[List[ConvergenceResult]] = None,
) -> List[ConvergenceResult]:
    """Pool ``M >= 1`` per-hypothesis ``compute_domain_scores`` records for one
    ``(Gene, Position, scheme)`` and emit <= 2 per-side ConvergenceResult rows.

    ``hyp_rows`` = ``[{"hyp": str, "sides": <compute_domain_scores return>}]``.
    Always calls :func:`fop_pool.pool_domains` (``M == 1`` degenerates to the
    plain PSS-weighted mean over the K domains). A position with no changed
    domain on either side across the harvest collapses to one ``side="none"``
    row. When the harvest carries more than one distinct hypothesis the emitted
    rows drop the ``hypothesis`` label (a genuine FOP pool); a lone hypothesis
    keeps it. ``all_rows`` (the raw per-hypothesis ``ConvergenceResult`` list,
    same length/order as ``hyp_rows``) is used only for the position-level
    ``tag_support``/``caas``/``amino_encoded`` derivation -- these live on the
    result itself, not inside ``sides``, so they cannot be recovered from
    ``hyp_rows`` alone.
    """
    pooled = pool_domains(hyp_rows, hyp_pairs_pss)
    meta = None
    for hr in hyp_rows:
        meta = (hr.get("sides") or {}).get("domain_meta") or meta
    hyp_labels = {hr.get("hyp") for hr in hyp_rows if hr.get("hyp")}
    hyp_out = None if len(hyp_labels) > 1 else base.hypothesis

    def _tally(attr: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for r in (all_rows or []):
            v = getattr(r, attr, None)
            if v:
                counts[str(v)] = counts.get(str(v), 0) + 1
        return counts

    tag_support = fmt_support(_tally("tag"))
    caas_support = fmt_support(_tally("caas"))
    amino_encoded_support = fmt_support(_tally("amino_encoded"))
    # `caas`/`amino_encoded` are no longer an arbitrary first-hypothesis
    # passthrough: derive them as the union of divergent (non-conserved)
    # residues across every pooled hypothesis (see _derive_convergent_call).
    # Fall back to the first-hypothesis passthrough only if no divergent
    # residue survived exclusion (e.g. every row's harvest fully agrees with
    # its own background) -- never emit a blank caas/amino_encoded.
    derived_caas = _derive_convergent_call(all_rows or [], "caas") or base.caas
    derived_amino_encoded = (
        _derive_convergent_call(all_rows or [], "amino_encoded") or base.amino_encoded
    )

    # Harvest size (M) for this (position, scheme) pool.
    n_hypotheses = int(pooled.get("n_hypotheses", 0) or 0)

    sides = [s for s in ("top", "bottom")
             if int((pooled.get(s) or {}).get("n_participating", 0) or 0) > 0]
    if not sides:
        return [dataclasses.replace(
            base, side="none", hypothesis=hyp_out, participating_hypotheses=None,
            asr_path_score=0.0, convergence_type="no_change",
            domain_scores=None, domain_anc_aa=None,
            domain_der_top_aa=None, domain_der_bot_aa=None,
            domain_der_support_top_aa=None, domain_der_support_bot_aa=None,
            domain_anc_support_aa=None, pair_lca=None,
            domain_meta=(dict(meta) if meta else None),
            caas=derived_caas, amino_encoded=derived_amino_encoded,
            tag_support=tag_support, caas_support=caas_support,
            amino_encoded_support=amino_encoded_support,
            n_hypotheses=n_hypotheses,
        )]

    out: List[ConvergenceResult] = []
    for s in sides:
        d = pooled[s]
        den = int(d.get("agree_den", 0) or 0)
        da = (int(d.get("agree_num", 0) or 0) / den) if den else None
        out.append(dataclasses.replace(
            base, side=s, hypothesis=hyp_out,
            participating_hypotheses=(",".join(d.get("participating_hyps") or []) or None),
            asr_path_score=float(d.get("asr_path_score", 0.0) or 0.0),
            derived_agreement=da,
            convergence_type=d.get("convergence_type", base.convergence_type),
            domain_scores=(dict(d.get("domain_scores") or {}) or None),
            domain_anc_aa=(dict(d.get("domain_anc") or {}) or None),
            domain_der_top_aa=(dict(d.get("domain_der") or {}) or None) if s == "top" else None,
            domain_der_bot_aa=(dict(d.get("domain_der") or {}) or None) if s == "bottom" else None,
            domain_der_support_top_aa=(dict(d.get("domain_der_support") or {}) or None) if s == "top" else None,
            domain_der_support_bot_aa=(dict(d.get("domain_der_support") or {}) or None) if s == "bottom" else None,
            domain_anc_support_aa=(dict(d.get("domain_anc_support") or {}) or None),
            pair_lca=_pooled_pair_lca(hyp_rows, s),
            domain_meta=(dict(meta) if meta else None),
            caas=derived_caas, amino_encoded=derived_amino_encoded,
            tag_support=tag_support, caas_support=caas_support,
            amino_encoded_support=amino_encoded_support,
            n_hypotheses=n_hypotheses,
        ))
    return out


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
) -> Dict[str, Any]:
    """Reduced domain-score kernel shared by the full scorer and the perm replay.

    The return is always ``compute_domain_scores``'s
    ``{"top": {...}, "bottom": {...}, "domain_meta": {...}}``.

    Builds ``per_node_dist`` for the focal site directly from the posterior map
    and runs :func:`compute_domain_scores`, so the observed path and the
    axes-only perm path score each position through identical code.

    ``per_site_dist_cache`` (perm replay only): ``per_node_dist`` depends solely
    on ``(posterior_data, site)`` — invariant to the grouping scheme and to the
    permuted phenotype — so a recurring site is built once and reused. Keyed by
    ``position_one_based``; the observed path passes ``None``.
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

    return compute_domain_scores(
        pair_details=pair_details_list,
        per_node_dist=per_node_dist,
        node_index=node_index,
        scheme=getattr(caas_pos, "caap_group", "US") or "US",
    )


def analyze_caas_position_disambiguation(
    gene: str,
    caas_pos: CAASPosition,
    tree_data,
    posterior_data: Optional[dict] = None,
    tip_diagnostics: Optional[Dict[str, Any]] = None,
    posterior_threshold: float = 0.7,
    node_index: Optional[Dict[int, Any]] = None,
    build_node_posteriors: bool = False,
    per_site_dist_cache: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    hypothesis: Optional[str] = None,
) -> "List[ConvergenceResult]":
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
        posterior_data: ASR posterior probabilities
        posterior_threshold: Posterior probability threshold for accepting node states

    Returns:
        List of per-side ConvergenceResult rows (1-2; a "both" position yields
        one row per participating side, a no-change position one ``side="none"``
        row).
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

    # ``side`` / ``convergence_type`` are no longer produced by a separate
    # categorical classifier — ``side`` is set per row in
    # :func:`_split_result_by_side` and ``convergence_type`` comes from
    # ``compute_asr_path_score``'s ``_convergence_type``. This is the placeholder
    # on ``base_result`` for the no-change collapse.
    convergence_type = "no_change"

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

    # ── CAAS convergence score (core v3, on the Voronoi domain) ───────────────
    # ``compute_domain_scores`` returns {"top", "bottom", "domain_meta"} for this
    # (Gene, Position, scheme, hypothesis). The per-side pooling + the <=2-row
    # split happen once per (position, scheme) in analyze_gene_disambiguation via
    # :func:`_emit_pooled_side_rows` (M == 1 degenerates to the plain mean). Here
    # we only stash the raw record on ``base_result.sides``.
    domain_split: Optional[Dict[str, Any]] = None
    try:
        domain_split = _position_axes(
            caas_pos, tree_data, posterior_data, node_index, pair_details_list,
            per_site_dist_cache=per_site_dist_cache,
        )
    except Exception as e:  # never let path scoring break disambiguation
        domain_split = None
        logger.warning(
            f"ASR path scoring failed for {gene}:{caas_pos.position}: {e}"
        )

    base_result = ConvergenceResult(
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
        side="none",  # overwritten per-side by _emit_pooled_side_rows
        caap_group=getattr(caas_pos, "caap_group", "US"),
        amino_encoded=getattr(caas_pos, "amino_encoded", ""),
        hypothesis=hypothesis,
        asr_path_score=None,
        derived_agreement=None,
        domain_scores=None,
        domain_anc_aa=None,
        domain_der_top_aa=None,
        domain_der_bot_aa=None,
        domain_meta=(dict(domain_split.get("domain_meta"))
                     if domain_split and domain_split.get("domain_meta") else None),
        score=None,
    )

    # Raw compute_domain_scores record for this (position, scheme, hypothesis);
    # analyze_gene_disambiguation groups these by (position, scheme) and pools
    # them with _emit_pooled_side_rows.
    base_result.sides = domain_split
    return [base_result]


def analyze_gene_disambiguation(
    gene: str,
    alignment_data,
    tree_data,
    caas_positions: List[int],
    caas_entries: Optional[List[CAASPosition]] = None,
    caas_metadata_path: Optional[Path] = None,
    trait_file_path: Optional[Path] = None,
    trait_pairs: Optional[Dict[int, List[Tuple[str, str]]]] = None,
    taxid_mapping: Optional[Dict[str, str]] = None,
    posterior_data: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    posterior_threshold: float = 0.7,
    diagnostics_dir: Optional[Path] = None,
    asr_mode: str = "precomputed",
    axes_only: bool = False,
    per_site_dist_cache: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None,
    build_node_posteriors: bool = False,
    hyp_pairs_pss: Optional[Dict[Tuple[str, int], float]] = None,
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
        trait_pairs: Optional pre-parsed {contrast: [(high_species, low_species), ...]},
            the exact return shape of parse_trait_pairs. Takes precedence over
            trait_file_path when given, skipping the file write+reparse round-trip
            the permulation-null replay would otherwise do per cycle.
        taxid_mapping: Optional species to taxid mapping
        posterior_data: Optional ASR posterior data
        posterior_threshold: Posterior probability threshold for node state extraction
        axes_only: Reduced-kernel mode for the permulation replay. When True, each
            position still builds pair_details, but SKIPS the full per-position
            scorer (its redundant per-node posterior-map rebuild + the 40-field
            ConvergenceResult the perm null discards). It emits a lightweight
            :class:`PositionAxes` per position carrying the per-side
            ``compute_asr_path_score`` return (``.sides``). The asr_path_score is
            scored by the same :func:`_position_axes` helper as the full path, so
            it is bit-for-bit identical.

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
    # LCA merge points, independence, mrca_diversity — with contrasts from
    # unrelated hypotheses, and discarded the per-hypothesis Dunn independence
    # the FOP harvest enforces.
    trait_pairs_all: Dict[int, List[Tuple[str, str]]] = {}
    if trait_pairs is not None:
        trait_pairs_all = trait_pairs
    elif trait_file_path:
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
    # get_mrca's tip lookups (find_node_by_name/find_node_by_taxid) are an
    # unindexed O(tree size) recursive search each; confirmed via real-data
    # cProfile as ~42% of total replay wall time (see
    # docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md, Tier 3) since the null
    # replay calls this once per cycle over a tree that never changes. Building
    # the index is itself a single O(tree size) pass, so even rebuilding it on
    # every call here (unavoidable without threading it in from the per-gene
    # caller) turns many O(tree size) lookups into one.
    _name_index, _taxid_index = (
        build_name_taxid_index(tree_data.root) if getattr(tree_data, "root", None) else ({}, {})
    )
    _mrca_cache: Dict[tuple, Any] = {}

    def _get_mrca_cached(taxa: List[str]):
        if not tree_data or not taxa:
            return None
        key = tuple(sorted(taxa))
        if key not in _mrca_cache:
            _mrca_cache[key] = get_mrca(
                tree_data.root, list(taxa),
                name_index=_name_index, taxid_index=_taxid_index,
            )
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
                    node_mapping = {
                        "root": (
                            tree_data.root.node_id
                            if tree_data and tree_data.root
                            else None
                        ),
                        "mrca_contrast": mrca_node.node_id if mrca_node else None,
                        "focal_nodes": [p.get("node_id") for p in pair_details],
                    }

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
            # The perm null keeps only asr_path_score + the per-side split, so we
            # skip the full per-position scorer — whose cost is the redundant per-node
            # posterior-map rebuild (node_posteriors["per_node"]) plus the 40-field
            # ConvergenceResult assembly the null discards — and call the shared
            # _position_axes helper directly (which builds per_node_dist once, cached).
            # The `.sides` field carries the whole {"top": row, "bottom": row}
            # split; gene_wrapper._expand_sides derives per-side change labels from
            # it. The scalar fields here collapse to the stronger side as a legacy
            # fallback.
            if axes_only:
                # NOTE (scoring_v2 core v3): the permulation-null replay
                # (_perms_worker) is rewired in V3-3. This branch now only stashes
                # the raw compute_domain_scores record on ``.sides``; the pooled
                # scalars are recomputed downstream by pool_domains.
                axes_sides = None
                axes_score = 0.0
                try:
                    axes_sides = _position_axes(
                        caas_pos,
                        tree_data,
                        posterior_data,
                        hoisted_node_index,
                        tip_diagnostics.get("pair_details"),
                        per_site_dist_cache=per_site_dist_cache,
                    )
                except Exception as e:  # never let path scoring break the replay
                    logger.warning(
                        f"ASR path scoring failed for {gene}:{caas_pos.position}: {e}"
                    )
                    axes_sides = None
                results.append(
                    PositionAxes(
                        position=caas_pos.position,
                        caap_group=getattr(caas_pos, "caap_group", "US"),
                        asr_path_score=axes_score,
                        side="none",
                        hypothesis=_hyp_label,
                        domain_scores=None,
                        derived_agreement=None,
                        sides=axes_sides,
                    )
                )
                logger.debug(f"✓ Axes-only position {pos}")
                continue

            # Perform convergence/disambiguation analysis
            row_list = analyze_caas_position_disambiguation(
                gene,
                caas_pos,
                tree_data,
                posterior_data,
                tip_diagnostics,
                posterior_threshold=posterior_threshold,
                node_index=hoisted_node_index,
                build_node_posteriors=build_node_posteriors,
                per_site_dist_cache=per_site_dist_cache,
                hypothesis=_hyp_label,
            )

            # One base row per (position, scheme, hypothesis); the per-side
            # split + pooling happens after the loop.
            results.extend(row_list)
            logger.info(
                f"✓ Analyzed position {pos}: {row_list[0].ancestral}→{row_list[0].derived}"
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

    if not axes_only and results:
        # core v3: group the per-hypothesis base rows by (position, scheme) and
        # pool them onto <=2 per-side rows with the treeless mean-of-means pooler
        # (M == 1 degenerates to the plain PSS-weighted mean over the K domains).
        # PSS weights {(hyp, domain): pss} come from hyp_pairs_pss (None -> equal
        # weight). Scoring receives rows already domain-pooled here.
        try:
            by_group: Dict[Tuple[Any, str], List[ConvergenceResult]] = {}
            for r in results:
                by_group.setdefault(
                    (getattr(r, "position", None), getattr(r, "caap_group", "US")), []
                ).append(r)
            pooled_rows: List[ConvergenceResult] = []
            for _rows in by_group.values():
                hyp_rows = [
                    {"hyp": (getattr(r, "hypothesis", None) or f"_H_UNRESOLVED_{i}"),
                     "sides": getattr(r, "sides", None) or {}}
                    for i, r in enumerate(_rows)
                ]
                pooled_rows.extend(
                    _emit_pooled_side_rows(_rows[0], hyp_rows, hyp_pairs_pss, all_rows=_rows)
                )
            results = pooled_rows
        except Exception as e:  # never let pooling break disambiguation
            logger.warning(f"[{gene}] domain pooling failed: {e}")

    logger.info(
        f"✓ Completed convergence disambiguation: {len(results)}/{len(caas_entries)} metadata rows"
    )
    diagnostics["skip_reasons"] = dict(diagnostics["skip_reasons"])
    return results, diagnostics


# Backward-compatible aliases
analyze_caas_position_biochemistry = analyze_caas_position_disambiguation
analyze_gene_biochemistry = analyze_gene_disambiguation
