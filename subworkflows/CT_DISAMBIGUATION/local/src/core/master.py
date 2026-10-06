# master.py — Schema and writer of caas_convergence_master.csv for the observed labeling.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
The master CSV of the observed labeling: rows are serialized once, in the worker, and written once.

Imported by: contract_main.py, observed_b0_main.py, src/core/contract.py, src/core/observed.py
Inputs: record dicts built from the per-position results (see src/utils/gene_wrapper.py)
Outputs: caas_convergence_master.csv, or a per-gene .master.csv.gz shard when the path ends in .gz

caas_convergence_master.csv is the table the rest of the pipeline reads. Its rows are built from the records the workers
hold and written here ordered by gene, then msa_pos (rows of one position keep the order they were produced in).
"""
import csv
import gzip
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


def serialize_value(val: Any) -> str:
    """CSV cell of a record value: lists comma-joined, None empty, everything else str()."""
    if isinstance(val, (list, tuple)):
        return ",".join(str(v) for v in val)
    if val is None:
        return ""
    return str(val)


def master_fields(max_pairs: int) -> List[str]:
    """The columns of the master CSV: the fixed ones, then eight per Voronoi domain for `max_pairs` domains
    (the largest pair id of the observed design). Every batch and every route writes the same schema."""
    fields = [
        # Identification
        "gene",
        "msa_pos",
        "caas",
        # Pattern classification
        "convergence_type",
        # Grouping scheme of the row
        "caap_group",
        "amino_encoded",
        # `caas`/`amino_encoded` above are the union of divergent (non-conserved)
        # residues across every pooled hypothesis, not one hypothesis's raw
        # pattern (see disambiguate_single._derive_convergent_call). These
        # tallies carry the full per-hypothesis breakdown. There is no separate
        # row-identifier column: CT_ACCUMULATION reads `tag_support` as the row tag.
        "tag_support",
        "caas_support",
        "amino_encoded_support",
        # Hypotheses that drove >= 1 changed domain on THIS SIDE: the only
        # hypothesis-provenance column, side-aware (top/bottom rows for the
        # same position can legitimately differ). SCORING reads it by name
        # (scoring_compute.R).
        "participating_hypotheses",
        # Number of hypotheses (M) pooled for this (position, scheme), not per side.
        "n_hypotheses",
        # Direction key of the row (top / bottom / none).
        "side",
        # CAAS convergence score on the Voronoi domain
        "asr_path_score",
        "derived_agreement",
        # True when a changed domain's derived residue was tied and settled by convention (smallest residue).
        "agreement_ambiguous",
    ]

    # Per-domain columns for the K fixed Voronoi domains (at end).
    for idx in range(1, max_pairs + 1):
        fields.extend(
            [
                f"domain_{idx}_posterior",
                f"domain_{idx}_score",
                # Ancestral and derived residues of the domain (modal over the
                # pooled hypotheses); _top_aa / _bot_aa are empty when that side
                # did not change.
                f"domain_{idx}_anc_aa",
                f"domain_{idx}_top_aa",
                f"domain_{idx}_bot_aa",
                # Cross-hypothesis support tallies of the modal residues above.
                f"domain_{idx}_anc_aa_support",
                f"domain_{idx}_top_aa_support",
                f"domain_{idx}_bot_aa_support",
            ]
        )

    return fields


def master_row(record: Dict[str, Any], fields: Sequence[str]) -> Dict[str, str]:
    return {k: serialize_value(record.get(k)) for k in fields}


def write_master_csv(
    rows: Iterable[Tuple[str, Any, Dict[str, str]]], path: Path, fields: Sequence[str]
) -> int:
    """Write (gene, msa_pos, row) triples sorted by gene then msa_pos (stable: production order within
    a position is kept). Returns the number of rows written; an empty input writes just the header.
    A path ending in .gz is written compressed."""
    ordered = sorted(rows, key=lambda t: (t[0], -1 if t[1] is None else t[1]))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for _gene, _pos, row in ordered:
            writer.writerow(row)
    return len(ordered)
