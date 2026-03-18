#!/usr/bin/env python3
"""Disambiguation CSV Writers.

Writes CAAS convergence results to CSV files with dynamic schema supporting
variable numbers of pairs (1-N).

**REDESIGNED FOR DYNAMIC PAIRS (2025-12-05)**:
CSV columns are generated dynamically based on the maximum number of pairs
found across all results.

Author: ASR Integration
Date: 2025-12-03
Updated: 2025-12-05 (Dynamic pair support)
"""

import logging
import csv
from pathlib import Path
from typing import List, Dict, Optional
import sqlite3
from typing import Tuple
import json as _json

from src.utils.disambiguation_db import fetch_alignment_for_gene
from src.utils.gene_wrapper import convert_convergence_result_to_dict

logger = logging.getLogger(__name__)


def write_caas_convergence_csvs(
    results: List[Dict], output_dir: Path, max_pairs: Optional[int] = None
) -> List[Path]:
    """
    Write master, ambiguous, and no_change debug CAAS convergence CSVs.

    Args:
        results: List of CAAS result dictionaries
        output_dir: Output directory for CSV files
        max_pairs: Maximum number of pairs (auto-detected if None)

    Returns:
        List of paths to written files
    """
    logger.info(f"Writing CAAS convergence CSVs to {output_dir}")

    master_filename = output_dir / "caas_convergence_master.csv"
    no_change_filename = output_dir / "no_change_debug.csv"

    # Build comments for debugging ASR/fallback issues
    for result in results:
        result["comments"] = _build_comments(result)

    _write_csv(results, master_filename, max_pairs=max_pairs)

    # Export no_change cases for debugging
    no_change_results = [r for r in results if r.get("pattern_type") == "no_change"]
    if no_change_results:
        _write_csv(
            no_change_results,
            no_change_filename,
            max_pairs=max_pairs,
        )
        logger.info(
            f"  Wrote {len(no_change_results)} no_change debug cases to {no_change_filename}"
        )

    logger.info(f"  Wrote {len(results)} results to {master_filename}")
    return (
        [master_filename, no_change_filename]
        if no_change_results
        else [master_filename]
    )


def _build_comments(result: Dict) -> str:
    """
    Build debug comments string from result data.

    Captures fallback-to-parent cases (ASR states with posterior=0)
    and other edge cases for debugging.

    Args:
        result: Single CAAS result dictionary

    Returns:
        Comments string (empty if no issues detected)
    """
    comments = []

    # Check for missing ASR states (posterior=0 indicates fallback)
    if "pair_details" in result and result["pair_details"]:
        for pair in result["pair_details"]:
            pair_id = pair.get("pair_id", "?")
            posterior = pair.get("focal_posterior", pair.get("mrca_posterior"))
            state = pair.get("focal_state", pair.get("mrca_modal_aa"))

            if posterior is not None and posterior == 0:
                comments.append(f"pair{pair_id}_fallback:{state}")

    # Check for low ASR posteriors (< 0.5)
    idx = 1
    while f"mrca_{idx}_posterior" in result:
        posterior = result.get(f"mrca_{idx}_posterior")
        state = result.get(f"mrca_{idx}_state")
        if posterior is not None and posterior < 0.5:
            comments.append(f"pair{idx}_low_posterior:{posterior:.3f}")
        idx += 1

    # Check for insufficient changes
    if result.get("pattern_type") == "no_change":
        top_count = result.get("top_change_count", 0)
        bottom_count = result.get("bottom_change_count", 0)
        if top_count > 0 or bottom_count > 0:
            comments.append(
                f"insufficient_changes(top={top_count},bottom={bottom_count})"
            )

    return ";".join(comments) if comments else ""


def _detect_max_pairs(results: List[Dict]) -> int:
    """
    Detect maximum number of focal pairs across all results.

    Args:
        results: List of CAAS result dictionaries

    Returns:
        Maximum number of pairs found (minimum 1)
    """
    max_pairs = 0
    for result in results:
        # Count mrca_N_node columns present in result
        pair_count = 0
        idx = 1
        while f"mrca_{idx}_node" in result:
            pair_count += 1
            idx += 1
        max_pairs = max(max_pairs, pair_count)

    # Minimum of 1 for backward compatibility
    return max(max_pairs, 1)


def _generate_dynamic_fields(max_pairs: int) -> List[str]:
    """
    Generate field list with dynamic focal node columns.

    **DYNAMIC PAIR SUPPORT (2025-12-05)**:
    ASR-related fields moved to END to maintain core CSV structure.
    Only extends schema when ASR data available.

    Args:
        max_pairs: Maximum number of pairs to generate columns for

    Returns:
        List of field names for CSV header
    """
    fields = [
        # Core identification (stable structure)
        "gene",
        "msa_pos",
        "tag",
        "caas",
        "is_significant",
        "pvalue",
        "pvalue_boot",
        # Pattern classification
        "pattern_type",
        "convergence_description",
        "convergence_mode",
        # Metadata-driven convergence context
        "caap_group",
        "amino_encoded",
        "is_conserved_meta",
        "conserved_pair",
        "sig_hyp",
        "sig_perm",
        "sig_both",
        # Change tracking
        "top_change_type",
        "bottom_change_type",
        "change_side",
        "low_confidence_nodes",
        # Conserved-pair validation flag
        "asr_is_conserved",
        "asr_root_conserved",
        # Debug and comments
        "comments",
        # ASR fields (AT END - only present when ASR available)
        "all_mrca_state",
        "all_mrca_posterior",
        "all_mrca_node",
    ]

    # Add dynamic focal node columns for N pairs (ASR-related, at end)
    for idx in range(1, max_pairs + 1):
        fields.extend(
            [
                f"mrca_{idx}_node",
                f"mrca_{idx}_state",
                f"mrca_{idx}_posterior",
            ]
        )

    return fields


def _write_csv(
    results: List[Dict],
    filename: Path,
    max_pairs: Optional[int] = None,
) -> None:
    """
    Write results to CSV with dynamic field generation.

    Args:
        results: List of result dictionaries
        filename: Output file path
    """
    if not results:
        logger.warning(f"No results to write to {filename}; writing header only")

    # Detect maximum number of pairs in results if not supplied
    max_pairs = max_pairs or _detect_max_pairs(results)

    # Generate dynamic field list
    fields = _generate_dynamic_fields(max_pairs)

    # Serialize list fields to strings for CSV
    def serialize_value(val):
        if isinstance(val, (list, tuple)):
            return ",".join(str(v) for v in val)
        elif isinstance(val, bool):
            return str(val)
        elif val is None:
            return ""
        return str(val)

    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()

        # Serialize list values
        for result in results:
            serialized = {k: serialize_value(result.get(k)) for k in fields}
            writer.writerow(serialized)

    logger.debug(f"Wrote {len(results)} rows with {len(fields)} fields to {filename}")


def export_from_db(
    db_path: Path, output_dir: Path, max_pairs: Optional[int] = None
) -> Tuple[List[Path], Path]:
    """
    Export CAAS convergence master CSV, no_change debug CSV, and per-gene JSONs directly from the aggregation SQLite DB.

    Streams rows from DB to avoid loading all results into memory.

    Returns:
        (list_of_caas_files, summary_json)
    """
    logger.info(f"Exporting CAAS convergence outputs from DB: {db_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    master_filename = output_dir / "caas_convergence_master.csv"
    no_change_filename = output_dir / "no_change_debug.csv"
    json_dir = output_dir / "json_summaries"
    json_dir.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(str(db_path))
    try:
        # First pass: detect max pairs if not supplied (use stored pair_count for efficiency)
        if max_pairs is None:
            try:
                cur = conn.cursor()
                cur.execute("SELECT MAX(pair_count) FROM results")
                row = cur.fetchone()
                max_pairs = int(row[0]) if row and row[0] else 1
            except Exception:
                max_pairs = 1

        # Prepare CSV writers (streaming)
        master_fields = _generate_dynamic_fields(max_pairs)

        def serialize_value(val):
            if isinstance(val, (list, tuple)):
                return ",".join(str(v) for v in val)
            elif isinstance(val, bool):
                return str(val)
            elif val is None:
                return ""
            return str(val)

        master_f = open(master_filename, "w", newline="")
        master_writer = csv.DictWriter(
            master_f, fieldnames=master_fields, extrasaction="ignore"
        )
        master_writer.writeheader()

        no_change_f = open(no_change_filename, "w", newline="")
        no_change_writer = csv.DictWriter(
            no_change_f, fieldnames=master_fields, extrasaction="ignore"
        )
        no_change_writer.writeheader()

        # Iterate rows ordered by gene, msa_pos, id and write rows one-by-one.
        # This preserves Tag-level hypotheses even when they share the same msa_pos.
        cur = conn.cursor()

        current_gene = None
        gene_file = None
        per_gene_counts = {}
        total_positions = 0

        cur.execute(
            "SELECT id, gene, msa_pos, result_json FROM results ORDER BY gene, msa_pos, id"
        )
        for _, gene, msa_pos, result_json in cur.fetchall():
            if not result_json:
                continue
            try:
                result = _json.loads(result_json)
            except Exception:
                continue

            align_data = fetch_alignment_for_gene(conn, gene) or {}
            alignment = align_data.get("alignment")
            taxid_to_species = align_data.get("taxid_to_species")
            seq_by_id = align_data.get("seq_by_id")
            seq_by_species = align_data.get("seq_by_species")
            alignment_extras = align_data.get("alignment_extras")
            posterior_dump_jsonl = None
            if alignment_extras:
                posterior_dump_jsonl = alignment_extras.get("posterior_dump_jsonl")

            caas_dict = convert_convergence_result_to_dict(
                result,
                multi_hypothesis=None,
                alignment=alignment,
                seq_by_id=seq_by_id,
                seq_by_species=seq_by_species,
                trait_pairs=None,
                taxid_to_species=taxid_to_species,
            )
            try:
                caas_dict["comments"] = _build_comments(caas_dict)
            except Exception:
                caas_dict["comments"] = ""

            master_writer.writerow(
                {k: serialize_value(caas_dict.get(k)) for k in master_fields}
            )
            total_positions += 1
            per_gene_counts[gene] = per_gene_counts.get(gene, 0) + 1

            if caas_dict.get("pattern_type") == "no_change":
                no_change_writer.writerow(
                    {k: serialize_value(caas_dict.get(k)) for k in master_fields}
                )

            if current_gene != gene:
                if gene_file is not None:
                    gene_file.close()
                gene_file = open(
                    json_dir / f"{gene.lower()}_convergence_positions.jsonl",
                    "a",
                    encoding="utf-8",
                )
                current_gene = gene
            try:
                from src.reporting.disambiguation_json import (
                    extract_convergence_summary,
                )

                summary = extract_convergence_summary(caas_dict, max_pairs)
            except Exception:
                summary = {"gene": gene, "msa_pos": msa_pos}
            if posterior_dump_jsonl:
                summary["posterior_dump_jsonl"] = posterior_dump_jsonl
            if gene_file is not None:
                gene_file.write(_json.dumps(summary, ensure_ascii=False) + "\n")
            else:
                logger.error(f"gene_file is None for gene: {gene}")

        # flush final gene file
        if gene_file is not None:
            gene_file.close()

        # close files
        master_f.close()
        no_change_f.close()
        # Write aggregated summary JSON (compact)
        summary_path = output_dir / "caas_convergence_summary.json"
        summary_obj = {
            "metadata": {
                "num_genes": len(per_gene_counts),
                "total_positions": total_positions,
                "num_pairs": max_pairs,
                "schema_version": "2025-12-08_db_export",
            },
            "by_gene_counts": per_gene_counts,
        }
        with open(summary_path, "w", encoding="utf-8") as f:
            _json.dump(summary_obj, f, indent=2, ensure_ascii=False)

    finally:
        conn.close()

    caas_files = [master_filename]
    if Path(no_change_filename).exists():
        caas_files.append(no_change_filename)

    logger.info(
        f"Exported CAAS master CSV: {master_filename}; JSON summary: {summary_path}"
    )
    return caas_files, summary_path
