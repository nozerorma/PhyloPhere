#!/usr/bin/env python3
"""
Aggregated JSON export for multi-gene convergence analysis.

Creates compact, reporter-friendly JSON objects derived from
per-position dictionaries produced by the CAAS/ASR pipeline.

Design goals:
- Keep JSON structure stable for downstream plotting/UI.
- Avoid recomputing heavy biology here; assume upstream results already
  contain core fields (convergence_type, amino_encoded, pair_details, etc.).
- Be defensive with missing optional keys.

Author: ASR Integration
Date: 2025-12-06
Updated: 2025-12-09 (robust pair extraction)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def extract_pair_info(
    result_dict: Dict[str, Any], num_pairs: int
) -> List[Dict[str, Any]]:
    """
    Extract pair-level MRCA and tip information.

    Source for tip residues:
        1) pair_details (if present)

    Returns a list with length == num_pairs.
    """
    pairs: List[Dict[str, Any]] = []

    pair_details = result_dict.get("pair_details") or []

    for idx in range(1, num_pairs + 1):
        mrca_state = result_dict.get(f"domain_{idx}_state")

        pair_info: Dict[str, Any] = {
            "pair_id": idx,
            "mrca_node": result_dict.get(f"domain_{idx}_node"),
            "mrca_state": mrca_state,
            "mrca_posterior": result_dict.get(f"domain_{idx}_posterior"),
        }

        top_tip: Optional[str] = None
        bottom_tip: Optional[str] = None

        # Preferred source: structured pair_details
        if isinstance(pair_details, list) and idx - 1 < len(pair_details):
            pair = pair_details[idx - 1] or {}
            if isinstance(pair, dict):
                top_tip = pair.get("top_tip_residue") or pair.get("top_tip_mode")
                bottom_tip = pair.get("bottom_tip_residue") or pair.get(
                    "bottom_tip_mode"
                )

                top_tip_residues = pair.get("top_tip_residues")
                bottom_tip_residues = pair.get("bottom_tip_residues")

                if top_tip_residues:
                    trs = (
                        top_tip_residues
                        if isinstance(top_tip_residues, list)
                        else [top_tip_residues]
                    )
                    pair_info["top_tip_residues"] = [
                        {
                            "species": r.get("species"),
                            "taxid": r.get("taxid"),
                            "residue": r.get("residue"),
                        }
                        for r in trs
                        if isinstance(r, dict) and r.get("residue")
                    ]

                if bottom_tip_residues:
                    brs = (
                        bottom_tip_residues
                        if isinstance(bottom_tip_residues, list)
                        else [bottom_tip_residues]
                    )
                    pair_info["bottom_tip_residues"] = [
                        {
                            "species": r.get("species"),
                            "taxid": r.get("taxid"),
                            "residue": r.get("residue"),
                        }
                        for r in brs
                        if isinstance(r, dict) and r.get("residue")
                    ]

        if top_tip or bottom_tip:
            pair_info["tip_states"] = {"top": top_tip, "bottom": bottom_tip}

            if mrca_state and top_tip:
                pair_info["top_change"] = {
                    "from": mrca_state,
                    "to": top_tip,
                    "changed": mrca_state != top_tip,
                }

            if mrca_state and bottom_tip:
                pair_info["bottom_change"] = {
                    "from": mrca_state,
                    "to": bottom_tip,
                    "changed": mrca_state != bottom_tip,
                }

            top_changed = bool(mrca_state and top_tip and (mrca_state != top_tip))
            bottom_changed = bool(
                mrca_state and bottom_tip and (mrca_state != bottom_tip)
            )
            pair_info["pair_change"] = top_changed or bottom_changed

            pair_info["conserved"] = bool(
                top_tip and bottom_tip and (top_tip == bottom_tip)
            )

        pairs.append(pair_info)

    return pairs


def extract_convergence_summary(
    result_dict: Dict[str, Any], num_pairs: int
) -> Dict[str, Any]:
    """
    Extract a simplified convergence summary for one CAAS position.
    """
    summary: Dict[str, Any] = {
        "gene": result_dict.get("gene"),
        "msa_pos": result_dict.get("msa_pos"),
        "tag": result_dict.get("tag"),
        "convergence_type": result_dict.get("convergence_type"),
        "caas": result_dict.get("caas", ""),
        "caap_group": result_dict.get("caap_group", "US"),
        "amino_encoded": result_dict.get("amino_encoded", ""),
        # First-class direction key (top / bottom / none). T4b retired the
        # change_top/change_bottom/change_side triplet.
        "side": result_dict.get("side", "none"),
    }

    summary["pairs"] = extract_pair_info(result_dict, num_pairs)

    # Optional visualization payloads
    if result_dict.get("node_mapping"):
        summary["node_mapping"] = result_dict["node_mapping"]
    if result_dict.get("node_state_details"):
        summary["node_state_details"] = result_dict["node_state_details"]

    if result_dict.get("multi_hypothesis"):
        summary["multi_hypothesis"] = result_dict["multi_hypothesis"]

    return summary


