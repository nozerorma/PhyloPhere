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


def _generate_dynamic_fields(max_pairs: int) -> List[str]:
    """
    Generate field list with dynamic per-domain columns.

    Args:
        max_pairs: Maximum number of pairs to generate columns for

    Returns:
        List of field names for CSV header
    """
    fields = [
        # Core identification (stable structure)
        "gene",
        "msa_pos",
        "caas",
        # Pattern classification
        "convergence_type",
        # Metadata-driven convergence context
        "caap_group",
        "amino_encoded",
        # `caas`/`amino_encoded` above are the union of divergent (non-conserved)
        # residues across every pooled hypothesis, not one hypothesis's raw
        # pattern (see disambiguate_single._derive_convergent_call). These
        # tallies carry the full per-hypothesis breakdown; `tag_support` also
        # stands in for the row identifier no longer carried as its own column.
        "tag_support",
        "caas_support",
        "amino_encoded_support",
        # Hypotheses that drove >= 1 changed domain on THIS SIDE -- the sole
        # hypothesis-provenance column, side-aware (top/bottom rows for the
        # same position can legitimately differ); SCORING's pos_scores
        # aggregation consumes it by name.
        "participating_hypotheses",
        # Harvest size (M) for this (position, scheme) pool, not per-side.
        "n_hypotheses",
        # First-class direction key (top / bottom / none). T4b retired the
        # change_top/change_bottom/change_side triplet.
        "side",
        # CAAS convergence score on the Voronoi domain (scoring_v2 core v3)
        "asr_path_score",
        "derived_agreement",
    ]

    # Per-domain columns for the K fixed Voronoi domains (at end).
    for idx in range(1, max_pairs + 1):
        fields.extend(
            [
                f"domain_{idx}_posterior",
                f"domain_{idx}_score",
                # Raw derived/ancestral residues (modal over the harvest): feed
                # the FOP harvest-wide per-scheme derived_agreement rebuild.
                # _top_aa / _bot_aa empty when that side did not change.
                f"domain_{idx}_anc_aa",
                f"domain_{idx}_top_aa",
                f"domain_{idx}_bot_aa",
                # Cross-hypothesis support tallies for the modal residues above.
                f"domain_{idx}_anc_aa_support",
                f"domain_{idx}_top_aa_support",
                f"domain_{idx}_bot_aa_support",
            ]
        )

    return fields


def export_from_db(
    db_path: Path, output_dir: Path, max_pairs: Optional[int] = None
) -> Tuple[List[Path], Path]:
    """
    Export the decoration outputs (no_change debug CSV, per-gene JSONs, summary JSON) from the aggregation SQLite DB.

    The master CSV is not written here: core.master writes it from the workers' rows.
    Streams rows from DB to avoid loading all results into memory.

    Returns:
        (list_of_caas_files, summary_json)
    """
    logger.info(f"Exporting CAAS convergence outputs from DB: {db_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    diag_dir = output_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)
    no_change_filename = diag_dir / "no_change_debug.csv"
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

        # The no_change rows keep the master's column schema
        master_fields = _generate_dynamic_fields(max_pairs)

        from src.core.master import serialize_value

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

            total_positions += 1
            per_gene_counts[gene] = per_gene_counts.get(gene, 0) + 1

            if (caas_dict.get("side") or "none") == "none":
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

    caas_files = []
    if Path(no_change_filename).exists():
        caas_files.append(no_change_filename)

    logger.info(f"Exported no_change debug CSV and per-gene JSONs; JSON summary: {summary_path}")
    return caas_files, summary_path
