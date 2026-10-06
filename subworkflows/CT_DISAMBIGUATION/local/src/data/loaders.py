# loaders.py — Small loaders and parsers of discovery-row fields shared by the observed and null sides.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/data/

"""
Data loaders for the CAAS scorer.

Imported by: explain_positions.py, observed_b0_main.py, src/core/observed.py, src/utils/gene_wrapper.py
Inputs: an Ensembl genes TSV/CSV with a `gene` column; discovery-row cell values
Outputs: in-memory sets and normalized values

Small helpers shared by the observed and the null side:

- **as_bool**: the pipeline's boolean spelling of a flag value
- **load_ensembl_genes**: Ensembl gene names from a TSV/CSV file with a `gene` column
- **_parse_conserved_pair**: the `conserved_pair` cell of a discovery row, as the comma-joined text the master holds

Trait/contrast definitions are read by `src.core.labelings` (one parser for the observed and the null labelings); the CAAS
entries come from the discovery rows (`src.core.observed.observed_entries`) or from the permulation exports.
"""

import csv
import logging
from pathlib import Path
from typing import Any, Set

logger = logging.getLogger(__name__)


def _parse_conserved_pair(raw: str) -> str:
    """Normalise the conserved_pair field from the CT output.

    caas_id.py (CT) writes the field as ``"{count}:{pair_id1},{pair_id2},..."``,
    e.g. ``"1:3"`` or ``"2:1,4"``. A value without the count prefix (a plain pair
    id such as ``"1"``) is returned unchanged.

    Returns a comma-separated string of pair ids, or ``""`` when no conserved
    pairs are present (``"0:"`` or empty input).
    """
    raw = raw.strip()
    if not raw:
        return ""
    if ":" in raw:
        # "{count}:{pairs}": drop the count prefix
        pairs_part = raw.split(":", 1)[1]
        # "0:" gives "" (no conserved pairs)
        return pairs_part.strip()
    return raw  # already a plain id or comma-separated ids


# ── Metadata parsing ──────────────────────────────────────────────────────────


def as_bool(v: Any) -> bool:
    """A metadata cell as a boolean: True/1/yes/y in any case; None and everything else False."""
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    return str(v).strip().lower() in {"true", "1", "yes", "y"}


# ── Ensembl genes ─────────────────────────────────────────────────────────────


def load_ensembl_genes(ensembl_genes_file: Path) -> Set[str]:
    """Load Ensembl gene names from a TSV/CSV file (expects a 'gene' column)."""
    if not ensembl_genes_file.exists():
        raise FileNotFoundError(f"Ensembl genes file not found: {ensembl_genes_file}")

    with ensembl_genes_file.open(newline="") as handle:
        sample = handle.read(2048)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        # reader.fieldnames is None without a header: check it before the membership test
        if not reader.fieldnames or "gene" not in reader.fieldnames:
            raise ValueError("Ensembl genes file must contain a 'gene' column")
        genes = {row["gene"].strip() for row in reader if row.get("gene")}
    if not genes:
        raise ValueError("No genes found in Ensembl genes file")
    return genes
