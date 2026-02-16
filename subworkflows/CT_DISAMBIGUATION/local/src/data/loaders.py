"""
Data Loaders for CAAS Analysis
===============================

Unified helpers for loading CAAS metadata and species pair definitions.
Trait parsing flows through a single parser (`parse_trait_pairs`) ensuring
consistent ordering and row handling across aggregation and single-gene pipelines.

Workflow
--------
1. **read_caas_metadata_table**: Load CAAS metadata from file, normalize columns, filter by gene if specified
2. **list_gene_caas_positions**: Extract all CAAS positions for a gene from metadata
3. **get_caas_position_info**: Retrieve metadata for a single CAAS position
4. **build_caas_positions_map**: Build CAASPosition objects for selected positions
5. **parse_trait_pairs**: Parse trait file into species pairs (no contrast column)
6. **load_ensembl_genes**: Load Ensembl gene names from file

Position Indexing Convention
------------------------------
Historical CAAS exports (e.g., `meta_caas.output`) record positions as **zero-based MSA
indices** in the `GenePos` column (e.g., `NUTM2A_648` → column 648). ASR consumers and
biology workflows expect **1-based** positions. All loaders return both forms:

    - position_zero_based: MSA column index (0-indexed)
    - position_one_based: Biological position (1-indexed)

This dual indexing prevents off-by-one errors in cross-consumer pipelines.

Data Contracts & Validation
------------------------------
1. **CAAS Metadata** (read_caas_metadata_table, get_caas_position_info):
     - Required columns: Tag, AminoConv, isSignificant
     - Optional columns: Pvalue, Pvalue.boot
     - GenePos normalization: Accepts 'GenePos' or 'Gene_Pos'; returns 'GenePos'
     - FileNotFoundError raised if file missing; ValueError for missing columns

2. **Trait Pairs** (parse_trait_pairs):
     - Tab-separated, no header: species, trait, pair columns
     - Deterministic ordering: Pairs sorted by numeric pair_id (ascending)
     - Incomplete rows (missing fields/invalid values) silently skipped with debug logging
     - Returns: [(high_species, low_species), ...]

Logging Strategy
------------------
- INFO: File loading summary, total counts (e.g., "Loaded X positions", "Y pairs")
- DEBUG: Per-item details (species→taxid lookup, pair details)
- WARNING: Missing files, invalid values, incomplete rows, tree validation failures
- Exception logging: Full traceback on critical failures (e.g., file read errors)

Usage Examples
---------------
::


        # Load CAAS metadata for a gene
        meta_file = Path('data/meta_caas.output')
        caas_positions = build_caas_positions_map(meta_file, gene='NUTM2A', positions=[85, 141])
        for pos_idx, caas_pos in caas_positions.items():
                print(f"{caas_pos.tag}: {caas_pos.trait1_aa} (significant: {caas_pos.is_significant})")

        # Load trait pairs with validation
        trait_file = Path('data/trait_pairs.txt')
        trait_pairs = parse_trait_pairs(trait_file)
        print(f"Loaded {len(trait_pairs)} pairs")

References
-----------
- CAASPosition dataclass: src.data.models.CAASPosition
- Amino acid normalization: src.convergence.convergence.normalize_amino_list()

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12
"""

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

import pandas as pd

from src.convergence.convergence import normalize_amino_list

from .models import CAASPosition

logger = logging.getLogger(__name__)


def _parse_gene_pos_token(gene_pos: str) -> Tuple[Optional[str], Optional[int]]:
    """Return (gene_name, zero_based_index) from a GenePos token like "NUTM2A_648".

    Returns (None, None) on malformed tokens.
    This is the single authoritative parser for GenePos-style tokens.
    """
    if not gene_pos:
        return None, None
    token = str(gene_pos).strip()
    if "_" not in token:
        return None, None
    gene, pos_str = token.rsplit("_", 1)
    try:
        pos_val = int(pos_str)
    except ValueError:
        return gene, None
    return gene, pos_val


#-- Functions for CAAS Metadata Loading and Parsing --#

def read_caas_metadata_table(
    metadata_file: Path, gene_name: Optional[str] = None
) -> pd.DataFrame:
    """
    Load CAAS metadata (Tag, AminoConv, isSignificant, optional Pvalue.boot).

    - Tries comma-separated first, falls back to tab-separated.
    - Accepts `GenePos` or `Gene_Pos` and normalizes to `GenePos`.
    - Returns a filtered DataFrame if `gene_name` is provided.
    - Positions in `GenePos` remain zero-based; callers may derive +1.
    """
    if not metadata_file.exists():
        raise FileNotFoundError(f"CAAS metadata file not found: {metadata_file}")

    try:
        # Try comma-separated first (typical for .output files), then tab-separated
        df = pd.read_csv(metadata_file, sep=",")
        # Check if parsing worked - if we got a single column, try tab-separated
        if len(df.columns) == 1:
            df = pd.read_csv(metadata_file, sep="\t")
    except Exception as e:
        raise ValueError(f"Error reading CAAS metadata file: {e}") from e

    # Validate required columns (flexible with Gene_Pos vs GenePos)
    required_cols = ["Tag", "AminoConv"]
    gene_pos_col = "GenePos" if "GenePos" in df.columns else "Gene_Pos"

    if gene_pos_col not in df.columns:
        raise ValueError(
            f"CAAS metadata file missing gene position column (GenePos or Gene_Pos). "
            f"Available columns: {list(df.columns)}"
        )

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"CAAS metadata file missing required columns: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    # Normalize column name for consistency
    if gene_pos_col == "Gene_Pos":
        df = df.rename(columns={"Gene_Pos": "GenePos"})

    # Harmonize significance columns to a single canonical boolean
    if "isSignificant" not in df.columns:
        if "sig_both" in df.columns:
            df["isSignificant"] = df["sig_both"]
        elif "sig_hyp" in df.columns:
            df["isSignificant"] = df["sig_hyp"]
        else:
            df["isSignificant"] = False

    # Ensure optional fields exist with stable names used downstream
    for col, default in {
        "CAAP_Group": "US",
        "AminoEncoded": "",
        "ConservedPair": "",
        "IsConserved": False,
        "sig_hyp": None,
        "sig_perm": None,
        "sig_both": None,
    }.items():
        if col not in df.columns:
            df[col] = default

    # Filter by gene if specified
    if gene_name:
        if "Gene" in df.columns:
            df = df[df["Gene"].astype(str) == str(gene_name)].copy()
        else:
            df = df[df["GenePos"].str.startswith(f"{gene_name}_")].copy()

    return df


def list_gene_caas_positions(caas_metadata_path: Path, gene: str) -> List[int]:
    """
    Load all CAAS positions available for a specific gene from metadata file.

    Args:
        caas_metadata_path: Path to CAAS metadata (.output) file
        gene: Gene name (e.g., "NUTM2A")

    Returns:
        List of zero-based MSA positions for the gene
    """
    logger.info(f"Loading CAAS positions for {gene} from {caas_metadata_path}")
    df = read_caas_metadata_table(caas_metadata_path, gene)
    logger.debug(f"Metadata has {len(df)} entries for {gene}")

    positions = []

    for gene_pos in df['GenePos']:
        # Use canonical parser for GenePos token
        gene_name_parsed, pos = _parse_gene_pos_token(gene_pos)
        if gene_name_parsed == gene and pos is not None:
            positions.append(pos)

    logger.info(f"Found {len(positions)} positions for {gene}: {sorted(positions[:5])}...")
    return sorted(set(positions))


def list_gene_caas_entries(caas_metadata_path: Path, gene: str) -> List[CAASPosition]:
    """
    Load CAAS metadata rows for a gene as independent CAASPosition entries.

    Each metadata row is treated as an independent hypothesis row
    (e.g., same GenePos with different CAAP_Group values).
    """
    logger.info(f"Loading CAAS entries for {gene} from {caas_metadata_path}")
    df = read_caas_metadata_table(caas_metadata_path, gene)

    entries: List[CAASPosition] = []
    for _, row in df.iterrows():
        gene_pos = str(row.get("GenePos", ""))
        _, pos0 = _parse_gene_pos_token(gene_pos)
        if pos0 is None:
            pos_raw = row.get("Position")
            if pd.notna(pos_raw):
                try:
                    pos0 = int(pos_raw) - 1
                except Exception:
                    pos0 = None
        if pos0 is None:
            continue

        caas = str(row.get("AminoConv", "") or "")
        parts = caas.split("/") if caas else []
        trait1 = normalize_amino_list(list(parts[0])) if len(parts) == 2 else []
        trait0 = normalize_amino_list(list(parts[1])) if len(parts) == 2 else []

        def _b(v: Any) -> bool:
            if isinstance(v, bool):
                return v
            if v is None:
                return False
            return str(v).strip().lower() in {"true", "1", "yes", "y"}

        entry = CAASPosition(
            position=pos0,
            position_zero_based=pos0,
            position_one_based=pos0 + 1,
            tag=str(row.get("Tag", f"POS{pos0}")),
            caas=caas,
            trait1_aa=trait1,
            trait0_aa=trait0,
            pvalue=float(row["Pvalue"]) if "Pvalue" in row.index and pd.notna(row["Pvalue"]) else None,
            pvalue_boot=float(row["Pvalue.boot"]) if "Pvalue.boot" in row.index and pd.notna(row["Pvalue.boot"]) else None,
            is_significant=_b(row.get("isSignificant")),
            caap_group=str(row.get("CAAP_Group", "US") or "US"),
            amino_encoded=str(row.get("AminoEncoded", "") or ""),
            is_conserved_meta=_b(row.get("IsConserved")),
            conserved_pair=str(row.get("ConservedPair", "") or ""),
            sig_hyp=_b(row.get("sig_hyp")) if pd.notna(row.get("sig_hyp")) else None,
            sig_perm=_b(row.get("sig_perm")) if pd.notna(row.get("sig_perm")) else None,
            sig_both=_b(row.get("sig_both")) if pd.notna(row.get("sig_both")) else None,
        )
        entries.append(entry)

    logger.info("Loaded %d metadata rows for %s", len(entries), gene)
    return entries


def get_caas_position_info(
    metadata_file: Path, gene_name: str, position: int
) -> Optional[Dict[str, Any]]:
    """
    Get metadata for a single CAAS position (zero-based index).

    Args:
        metadata_file: Path to CAAS metadata file
        gene_name: Gene name (e.g., "NUTM2A")
        position: Position as recorded in metadata (zero-based).

    Returns:
        Dict with tag, caas, pvalue, significance, and both
        zero/one-based positions, or None if not found.

    Example:
        >>> info = get_caas_position_info(meta_file, "NUTM2A", 85)
        >>> print(f"Significant: {info['is_significant']}")
        Significant: False
    """
    df = read_caas_metadata_table(metadata_file, gene_name)

    # Construct gene_pos identifier
    gene_pos = f"{gene_name}_{position}"

    # Find matching row
    matching = df[df["GenePos"] == gene_pos]

    if len(matching) == 0:
        logger.warning("No CAAS metadata found for %s", gene_pos)
        return None

    if len(matching) > 1:
        logger.warning("Multiple CAAS entries found for %s, using first", gene_pos)

    row = matching.iloc[0]

    # Extract core information (always present)
    info = {
        "gene_pos": str(row["GenePos"]),
        "tag": str(row["Tag"]),
        "caas": str(row["AminoConv"]),
        "is_significant": str(row["isSignificant"]).upper() == "TRUE",
    }

    _, zero_based_pos = _parse_gene_pos_token(str(info["gene_pos"]))
    info["position_zero_based"] = zero_based_pos
    info["position_one_based"] = (
        zero_based_pos + 1 if zero_based_pos is not None else None
    )

    # Add Pvalue (should always be present)
    if "Pvalue" in row.index:
        try:
            info["pvalue"] = float(row["Pvalue"])
        except (ValueError, TypeError):
            logger.warning("Could not parse Pvalue for %s: %s", gene_pos, row["Pvalue"])
            info["pvalue"] = None

    # Add Pvalue.boot if present (optional column)
    if "Pvalue.boot" in row.index:
        try:
            pvalue_boot = row["Pvalue.boot"]
            if pd.notna(pvalue_boot):
                info["pvalue_boot"] = float(pvalue_boot)
            else:
                info["pvalue_boot"] = None
        except (ValueError, TypeError):
            logger.warning(
                "Could not parse Pvalue.boot for %s: %s", gene_pos, row["Pvalue.boot"]
            )
            info["pvalue_boot"] = None

    logger.debug(
        "Found CAAS info for %s: tag=%s, significant=%s",
        gene_pos,
        info["tag"],
        info["is_significant"],
    )

    return info


def build_caas_positions_map(
    caas_metadata_path, gene: str, positions: List[int]
) -> Dict[int, CAASPosition]:
    """
    Load CAAS metadata for selected positions and return CAASPosition objects.

    Parses caas into normalized trait1/trait0 amino lists when present.
    Missing positions are skipped with a warning.
    """
    # Ensure we have a Path object
    caas_metadata_path = Path(caas_metadata_path)

    caas_positions: Dict[int, CAASPosition] = {}

    if not caas_metadata_path.exists():
        logger.warning("CAAS metadata file not found: %s", caas_metadata_path)
        return caas_positions

    for pos in positions:
        try:
            info = get_caas_position_info(caas_metadata_path, gene_name=gene, position=pos)
            if info:
                caas_pos = CAASPosition(
                    position=pos,
                    position_zero_based=info.get("position_zero_based", pos),
                    position_one_based=info.get("position_one_based", pos + 1),
                    tag=info.get("tag", f"POS{pos}"),
                    caas=info.get("caas", ""),
                    trait1_aa=[],
                    trait0_aa=[],
                    pvalue=info.get("pvalue"),
                    pvalue_boot=info.get("pvalue_boot"),
                    is_significant=info.get("is_significant", False),
                )

                # Parse amino acid conversion string
                parts = caas_pos.caas.split("/")
                if len(parts) == 2:
                    caas_pos.trait1_aa = normalize_amino_list(
                        list(parts[0])
                    )  # High trait
                    caas_pos.trait0_aa = normalize_amino_list(
                        list(parts[1])
                    )  # Low trait

                caas_positions[pos] = caas_pos

        except Exception as e:
            logger.warning("Failed to load metadata for position %d: %s", pos, e)
            continue

    logger.info("✓ Loaded metadata for %d CAAS positions", len(caas_positions))
    return caas_positions


#-- Functions for Trait Pair and Contrast Definition Loading --#

def parse_trait_pairs(
    trait_file_path: Path,
) -> List[Tuple[str, str]]:
    """
    Parse trait file once and return species pairs.

    - Tab-separated input with NO header and exactly three columns:
      species, trait, pair
    - Returns [(high_species, low_species), ...].
    - Pairs are sorted by numeric pair_id where possible.
    - Ignores rows with missing fields or invalid trait values.
    - No taxid mapping or validation here (pure parsing)
    """
    if not trait_file_path.exists():
        logger.error("Trait file not found: %s", trait_file_path)
        raise FileNotFoundError(f"Trait file not found: {trait_file_path}")

    # Structure: pair_id -> {'high': [species], 'low': [species]}
    by_pair: Dict[str, Dict[str, list]] = defaultdict(lambda: {"high": [], "low": []})

    try:
        with open(trait_file_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")

            for row in reader:
                # Expected row format (no header): species, trait, pair
                if len(row) < 3:
                    continue

                species = str(row[0]).strip()
                trait = str(row[1]).strip()
                pair_id = str(row[2]).strip()

                if not all([species, trait, pair_id]):
                    continue

                # Skip legacy header-like rows if present by mistake
                if species.lower() == "species" and trait.lower() == "trait":
                    continue

                try:
                    trait_val = int(trait)
                except (ValueError, TypeError):
                    logger.debug("Skipping row with invalid trait value: %s", row)
                    continue

                if trait_val == 1:
                    by_pair[pair_id]["high"].append(species)
                elif trait_val == 0:
                    by_pair[pair_id]["low"].append(species)

        # Build final ordered pair list
        pairs_list: List[Tuple[str, str]] = []
        for pair_id in sorted(by_pair.keys(), key=lambda x: int(x) if x.isdigit() else x):
            group = by_pair[pair_id]
            high_species = group["high"]
            low_species = group["low"]

            if high_species and low_species:
                pairs_list.append((high_species[0], low_species[0]))

                if len(high_species) > 1 or len(low_species) > 1:
                    logger.warning(
                        "Pair %s has multiple species per side: high=%s, low=%s. Using first from each.",
                        pair_id,
                        high_species,
                        low_species,
                    )
            else:
                logger.debug(
                    "Skipping incomplete pair %s: high=%s, low=%s",
                    pair_id,
                    high_species,
                    low_species,
                )

        total_pairs = len(pairs_list)
        logger.info(
            "Loaded %d total pairs from %s",
            total_pairs,
            trait_file_path,
        )
        return pairs_list

    except Exception as e:
        logger.error("Failed to load trait pairs: %s", e, exc_info=True)
        return []


#-- Function to Load Ensembl Genes --#

def load_ensembl_genes(ensembl_genes_file: Path) -> Set[str]:
    """Load Ensembl gene names from a TSV/CSV file (expects a 'gene' column)."""
    if not ensembl_genes_file.exists():
        raise FileNotFoundError(f"Ensembl genes file not found: {ensembl_genes_file}")

    with ensembl_genes_file.open(newline='') as handle:
        sample = handle.read(2048)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        # reader.fieldnames may be None (no header) — guard against that before membership test
        if not reader.fieldnames or 'gene' not in reader.fieldnames:
            raise ValueError("Ensembl genes file must contain a 'gene' column")
        genes = {row['gene'].strip() for row in reader if row.get('gene')}
    if not genes:
        raise ValueError("No genes found in Ensembl genes file")
    return genes