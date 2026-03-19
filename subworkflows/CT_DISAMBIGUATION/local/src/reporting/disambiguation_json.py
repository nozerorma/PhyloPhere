#!/usr/bin/env python3
"""
Aggregated JSON export for multi-gene convergence analysis.

Creates compact, reporter-friendly JSON objects derived from
per-position dictionaries produced by the CAAS/ASR pipeline.

Design goals:
- Keep JSON structure stable for downstream plotting/UI.
- Avoid recomputing heavy biology here; assume upstream results already
  contain core fields (pattern_type, amino_encoded, pair_details, etc.).
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
        mrca_state = result_dict.get(f"mrca_{idx}_state")

        pair_info: Dict[str, Any] = {
            "pair_id": idx,
            "mrca_node": result_dict.get(f"mrca_{idx}_node"),
            "mrca_state": mrca_state,
            "mrca_posterior": result_dict.get(f"mrca_{idx}_posterior"),
            "fallback_depth": result_dict.get(f"mrca_{idx}_fallback_depth", 0),
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
        "position": result_dict.get("position"),
        "position_zero_based": result_dict.get("msa_pos"),
        "tag": result_dict.get("tag"),
        "pattern_type": result_dict.get("pattern_type"),
        "caas": result_dict.get("caas", ""),
        "caap_group": result_dict.get("caap_group", "US"),
        "amino_encoded": result_dict.get("amino_encoded", ""),
        "is_conserved_meta": bool(result_dict.get("is_conserved_meta", False)),
        "conserved_pair": result_dict.get("conserved_pair", ""),
        "sig_hyp": result_dict.get("sig_hyp"),
        "sig_perm": result_dict.get("sig_perm"),
        "sig_both": result_dict.get("sig_both"),
        "change_top": result_dict.get("change_top", "no_change"),
        "change_bottom": result_dict.get("change_bottom", "no_change"),
        "change_side": result_dict.get("change_side", "none"),
        "parallel_top": result_dict.get("parallel_top"),
        "parallel_bottom": result_dict.get("parallel_bottom"),
        "parallel_type": result_dict.get("parallel_type", "none"),
        "is_significant": result_dict.get("is_significant"),
        "ambiguous": result_dict.get("ambiguous"),
        "asr_is_conserved": result_dict.get("asr_is_conserved"),
        "asr_root_conserved": result_dict.get("asr_root_conserved"),
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


def detect_max_pairs(results: List[Dict[str, Any]]) -> int:
    """Detect maximum number of focal pairs across all results (minimum 1)."""
    max_pairs = 0
    for result in results:
        idx = 1
        while f"mrca_{idx}_node" in result:
            max_pairs = max(max_pairs, idx)
            idx += 1
    return max(max_pairs, 1)


def export_aggregated_convergence_json(
    caas_results: List[Dict[str, Any]], output_path: Path
) -> Path:
    """
    Export aggregated convergence summaries for all genes to one JSON file.
    """
    logger.info(f"Exporting aggregated convergence JSON: {output_path}")

    num_pairs = detect_max_pairs(caas_results)
    logger.info(f"Detected {num_pairs} pairs across all results")

    by_gene: Dict[str, List[Dict[str, Any]]] = {}
    for result in caas_results:
        gene = result.get("gene", "unknown")
        by_gene.setdefault(gene, []).append(
            extract_convergence_summary(result, num_pairs)
        )

    for gene in by_gene:
        by_gene[gene].sort(key=lambda x: (x.get("position") is None, x.get("position")))

    output = {
        "metadata": {
            "num_genes": len(by_gene),
            "total_positions": len(caas_results),
            "num_pairs": num_pairs,
            "schema_version": "2025-12-09_simplified",
        },
        "genes": by_gene,
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    logger.info(f"✓ Exported aggregated convergence JSON: {output_path}")
    logger.info(
        f"  Genes: {len(by_gene)}, Positions: {len(caas_results)}, Pairs: {num_pairs}"
    )
    return output_path


def export_gene_summaries_json(
    caas_results: List[Dict[str, Any]], output_dir: Path
) -> List[Path]:
    """
    Export one JSON file per gene containing simplified convergence summaries.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    num_pairs = detect_max_pairs(caas_results)

    by_gene: Dict[str, List[Dict[str, Any]]] = {}
    for result in caas_results:
        gene = result.get("gene", "unknown")
        by_gene.setdefault(gene, []).append(
            extract_convergence_summary(result, num_pairs)
        )

    written_paths: List[Path] = []
    for gene, summaries in sorted(by_gene.items()):
        summaries.sort(key=lambda x: (x.get("position") is None, x.get("position")))

        output_file = output_dir / f"{gene.lower()}_convergence_summary.json"
        gene_output = {
            "metadata": {
                "gene": gene,
                "num_positions": len(summaries),
                "num_pairs": num_pairs,
                "schema_version": "2025-12-09_gene_level",
            },
            "positions": summaries,
        }

        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(gene_output, f, indent=2, ensure_ascii=False, default=str)

        written_paths.append(output_file)
        logger.debug(f"  {gene}: {len(summaries)} positions → {output_file}")

    logger.info(f"✓ Exported {len(written_paths)} per-gene JSON summaries")
    return written_paths
