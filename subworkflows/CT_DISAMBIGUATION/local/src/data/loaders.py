"""
Data Loaders for CAAS Analysis
===============================

Unified helpers for loading CAAS metadata, trait/contrast definitions, and species pairs.
Trait parsing flows through a single parser (`parse_trait_pairs`) ensuring
consistent ordering and row handling across aggregation and single-gene pipelines.

Workflow
--------
1. **read_caas_metadata_table**: Load CAAS metadata from file, normalize columns, filter by gene if specified
2. **list_gene_caas_positions**: Extract all CAAS positions for a gene from metadata
3. **get_caas_position_info**: Retrieve metadata for a single CAAS position
4. **build_caas_positions_map**: Build CAASPosition objects for selected positions
5. **parse_trait_pairs**: Parse trait file into species pairs grouped by contrast
6. **load_ensembl_genes**: Load Ensembl gene names from file

Position Indexing Convention
------------------------------
Historical CAAS exports (e.g., `meta_caas.output`) record positions as **zero-based MSA
indices** in the `GenePos` column (e.g., `NUTM2A_648` → column 648). ASR consumers and
biology workflows expect **1-based** positions. All loaders return both forms:

    - position: MSA column index (0-indexed)
    - position_one_based: Biological position (1-indexed)

This dual indexing prevents off-by-one errors in cross-consumer pipelines.

Data Contracts & Validation
------------------------------
1. **CAAS Metadata** (read_caas_metadata_table, get_caas_position_info):
     - Required columns: tag, caas, is_significant, GenePos
     - Optional columns: pvalue, pvalue_boot
     - FileNotFoundError raised if file missing; ValueError for missing columns

2. **Trait Pairs** (parse_trait_pairs):
     - Tab-separated: species, contrast, trait, pair columns
     - Deterministic ordering: Pairs sorted by numeric pair_id (ascending)
     - Incomplete rows (missing fields/invalid values) silently skipped with debug logging
     - Returns: {contrast → [(high_species, low_species), ...]}

Logging Strategy
------------------
- INFO: File loading summary, total counts (e.g., "Loaded X positions", "Y pairs")
- DEBUG: Per-item details (position→contrast mapping, species→taxid lookup, pair details)
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
        contrast_defs = parse_trait_pairs(trait_file)
        print(f"Loaded {len(contrast_defs)} contrasts")

References
-----------
- CAASPosition dataclass: src.data.models.CAASPosition
- ContrastDefinition dataclass: src.data.models.ContrastDefinition
- Amino acid normalization: src.utils.amino.normalize_amino_list()

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12
"""

import csv
import io
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

import pandas as pd

from src.utils.amino import normalize_amino_list

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


def _parse_conserved_pair(raw: str) -> str:
    """Normalise the conserved_pair field from the CT output.

    caas_id.py writes the field as ``"{count}:{pair_id1},{pair_id2},..."``,
    e.g. ``"1:3"`` or ``"2:1,4"``.  Legacy data (max_conserved=0 runs or old
    tooling) may already contain just a plain pair id like ``"1"``.

    Returns a comma-separated string of pair ids, or ``""`` when no conserved
    pairs are present (``"0:"`` or empty input).
    """
    raw = raw.strip()
    if not raw:
        return ""
    if ":" in raw:
        # "{count}:{pairs}" — drop the count prefix
        pairs_part = raw.split(":", 1)[1]
        # "0:" → "" (no conserved pairs)
        return pairs_part.strip()
    return raw  # already a plain id or comma-separated ids (legacy)


# -- Functions for CAAS Metadata Loading and Parsing --#


def _stream_gene_filtered(
    metadata_file: Path, sep: str, gene: str
) -> Optional[str]:
    """Return ``header + only the rows for ``gene`` as a text buffer, or None.

    Streams the file line-by-line with a cheap substring prefilter, splitting and
    exact-matching only candidate lines. Matches on the ``Gene`` column if present,
    else on the ``GenePos`` ``"{gene}_"`` prefix — mirroring the post-read gene
    filter in :func:`read_caas_metadata_table`. Returns None (caller falls back to a
    full parse) if neither key column is in the header, so behaviour degrades safely.
    """
    with open(metadata_file, "r") as fh:
        header = fh.readline()
        if not header:
            return None
        cols = header.rstrip("\n").split(sep)
        gene_idx = cols.index("Gene") if "Gene" in cols else None
        gp_idx = cols.index("GenePos") if "GenePos" in cols else None
        if gene_idx is None and gp_idx is None:
            return None  # no key column → let the caller parse the whole file
        prefix = f"{gene}_"
        parts_out = [header]
        for line in fh:
            if gene not in line:  # C-level superset guard: no false negatives
                continue
            parts = line.rstrip("\n").split(sep)
            if gene_idx is not None and gene_idx < len(parts):
                if parts[gene_idx] == gene:
                    parts_out.append(line)
            elif gp_idx is not None and gp_idx < len(parts):
                if parts[gp_idx].startswith(prefix):
                    parts_out.append(line)
    return "".join(parts_out)


def read_caas_metadata_table(
    metadata_file: Path, gene_name: Optional[str] = None
) -> pd.DataFrame:
    """
    Load CAAS metadata (tag, caas, is_significant, optional pvalue_boot).

    - Tries comma-separated first, falls back to tab-separated.
    - Returns a filtered DataFrame if `gene_name` is provided.
    - Positions in `GenePos` remain zero-based; callers may derive +1.
    """
    if not metadata_file.exists():
        raise FileNotFoundError(f"CAAS metadata file not found: {metadata_file}")

    try:
        # Sniff the delimiter from the header line. meta_caas is TSV (write_tsv);
        # legacy .output files are CSV. A comma-first guess is fragile because a TSV
        # field may itself contain commas (e.g. conserved_pair="1,3"), which makes
        # the C parser raise on inconsistent field counts. Matches concatenate.py.
        with open(metadata_file, "r") as _fh:
            _header = _fh.readline()
        sep = "\t" if "\t" in _header else ","

        # Gene-filtered fast path: the permulation cycle metadata is genome-wide
        # (~150k rows, all genes), and every per-gene worker × cycle would otherwise
        # pandas-parse the WHOLE file just to keep one gene's ~hundreds of rows. Stream
        # the header + only the matching lines and hand that small buffer to pandas —
        # same rows, same dtypes-on-the-kept-subset, ~7× faster. The substring
        # prefilter is a pure superset guard (a row with Gene==gene necessarily
        # contains the token, so no false negatives); the exact column check below
        # drops false positives (e.g. "A2M" matching an "A2ML1" line).
        buf = _stream_gene_filtered(metadata_file, sep, gene_name) if gene_name else None
        if buf is not None:
            df = pd.read_csv(io.StringIO(buf), sep=sep)
        else:
            df = pd.read_csv(metadata_file, sep=sep)
    except Exception as e:
        raise ValueError(f"Error reading CAAS metadata file: {e}") from e

    # Validate required columns
    required_cols = ["tag", "caas", "is_significant", "GenePos"]

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"CAAS metadata file missing required columns: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    # Ensure optional fields exist with stable names used downstream
    for col, default in {
        "caap_group": "US",
        "amino_encoded": "",
        "conserved_pair": "",
        "is_conserved_meta": False,
        "sig_hyp": None,
        "sig_perm": None,
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

    for gene_pos in df["GenePos"]:
        # Use canonical parser for GenePos token
        gene_name_parsed, pos = _parse_gene_pos_token(gene_pos)
        if gene_name_parsed == gene and pos is not None:
            positions.append(pos)

    logger.info(
        f"Found {len(positions)} positions for {gene}: {sorted(positions[:5])}..."
    )
    return sorted(set(positions))


def list_gene_caas_entries(caas_metadata_path: Path, gene: str) -> List[CAASPosition]:
    """
    Load CAAS metadata rows for a gene as independent CAASPosition entries.

    Each metadata row is treated as an independent hypothesis row
    (e.g., same GenePos with different caap_group values).
    """
    logger.info(f"Loading CAAS entries for {gene} from {caas_metadata_path}")
    df = read_caas_metadata_table(caas_metadata_path, gene)

    def _b(v: Any) -> bool:
        if isinstance(v, bool):
            return v
        if v is None:
            return False
        return str(v).strip().lower() in {"true", "1", "yes", "y"}

    entries: List[CAASPosition] = []
    # ``to_dict("records")`` once is far cheaper than ``iterrows()`` + per-column
    # ``Series.__getitem__``/``.get`` (each of which builds/searches a Series). The
    # dict values are the same typed cell values, so ``pd.notna`` and membership
    # checks behave identically — bit-for-bit same entries, less pandas overhead.
    for row in df.to_dict("records"):
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

        caas = str(row.get("caas", "") or "")
        parts = caas.split("/") if caas else []
        trait1 = normalize_amino_list(list(parts[0])) if len(parts) == 2 else []
        trait0 = normalize_amino_list(list(parts[1])) if len(parts) == 2 else []

        entry = CAASPosition(
            position=pos0,
            position_one_based=pos0 + 1,
            tag=str(row.get("tag", f"POS{pos0}")),
            caas=caas,
            trait1_aa=trait1,
            trait0_aa=trait0,
            pvalue=(
                float(row["pvalue"])
                if "pvalue" in row and pd.notna(row["pvalue"])
                else None
            ),
            pvalue_boot=(
                float(row["pvalue_boot"])
                if "pvalue_boot" in row and pd.notna(row["pvalue_boot"])
                else None
            ),
            is_significant=_b(row.get("is_significant")),
            caap_group=str(row.get("caap_group", "US") or "US"),
            amino_encoded=str(row.get("amino_encoded", "") or ""),
            is_conserved_meta=_b(row.get("is_conserved_meta")),
            conserved_pair=_parse_conserved_pair(str(row.get("conserved_pair", "") or "")),
            sig_hyp=_b(row.get("sig_hyp")) if pd.notna(row.get("sig_hyp")) else None,
            sig_perm=_b(row.get("sig_perm")) if pd.notna(row.get("sig_perm")) else None,
            pvalue_fdr=(
                float(row["pvalue_fdr"])
                if "pvalue_fdr" in row and pd.notna(row["pvalue_fdr"])
                else None
            ),
            pvalue_boot_fdr=(
                float(row["pvalue_boot_fdr"])
                if "pvalue_boot_fdr" in row and pd.notna(row["pvalue_boot_fdr"])
                else None
            ),
            alpha_fdr=(
                float(row["alpha_fdr"])
                if "alpha_fdr" in row and pd.notna(row["alpha_fdr"])
                else None
            ),
            sig_hyp_fdr=_b(row.get("sig_hyp_fdr")) if pd.notna(row.get("sig_hyp_fdr")) else None,
            sig_perm_fdr=_b(row.get("sig_perm_fdr")) if pd.notna(row.get("sig_perm_fdr")) else None,
            is_significant_fdr=_b(row.get("is_significant_fdr")) if pd.notna(row.get("is_significant_fdr")) else None,
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
        Dict with tag, contrast, caas, pvalue, significance, and both
        zero/one-based positions, or None if not found.

    Example:
        >>> info = get_caas_position_info(meta_file, "NUTM2A", 85)
        >>> print(f"Contrast: {info['contrast']}, Significant: {info['is_significant']}")
        Contrast: 2, Significant: False
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
        "tag": str(row["tag"]),
        "caas": str(row["caas"]),
        "is_significant": str(row["is_significant"]).upper() == "TRUE",
    }

    _, zero_based_pos = _parse_gene_pos_token(str(info["gene_pos"]))
    info["position"] = zero_based_pos
    info["position_one_based"] = (
        zero_based_pos + 1 if zero_based_pos is not None else None
    )

    # Add pvalue (should always be present)
    if "pvalue" in row.index:
        try:
            info["pvalue"] = float(row["pvalue"])
        except (ValueError, TypeError):
            logger.warning("Could not parse pvalue for %s: %s", gene_pos, row["pvalue"])
            info["pvalue"] = None

    # Add pvalue_boot if present (optional column)
    if "pvalue_boot" in row.index:
        try:
            pvalue_boot = row["pvalue_boot"]
            if pd.notna(pvalue_boot):
                info["pvalue_boot"] = float(pvalue_boot)
            else:
                info["pvalue_boot"] = None
        except (ValueError, TypeError):
            logger.warning(
                "Could not parse pvalue_boot for %s: %s", gene_pos, row["pvalue_boot"]
            )
            info["pvalue_boot"] = None

    for fld in ["pvalue_fdr", "pvalue_boot_fdr", "alpha_fdr"]:
        if fld in row.index:
            try:
                val = row[fld]
                info[fld] = float(val) if pd.notna(val) else None
            except (ValueError, TypeError):
                info[fld] = None

    for fld in ["sig_hyp_fdr", "sig_perm_fdr", "is_significant_fdr"]:
        if fld in row.index:
            val = row[fld]
            info[fld] = str(val).upper() == "TRUE" if pd.notna(val) else None

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
            info = get_caas_position_info(
                caas_metadata_path, gene_name=gene, position=pos
            )
            if info:
                caas_pos = CAASPosition(
                    position=pos,
                    position_one_based=info.get("position_one_based", pos + 1),
                    tag=info.get("tag", f"POS{pos}"),
                    caas=info.get("caas", ""),
                    trait1_aa=[],
                    trait0_aa=[],
                    pvalue=info.get("pvalue"),
                    pvalue_boot=info.get("pvalue_boot"),
                    is_significant=info.get("is_significant", False),
                    pvalue_fdr=info.get("pvalue_fdr"),
                    pvalue_boot_fdr=info.get("pvalue_boot_fdr"),
                    sig_hyp_fdr=info.get("sig_hyp_fdr"),
                    sig_perm_fdr=info.get("sig_perm_fdr"),
                    is_significant_fdr=info.get("is_significant_fdr"),
                    alpha_fdr=info.get("alpha_fdr"),
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


# -- Functions for Trait Pair and Contrast Definition Loading --#


def parse_trait_pairs(
    trait_file_path: Path,
) -> Dict[int, List[Tuple[str, str]]]:
    """
    Parse trait file once and return species pairs grouped by contrast.

    Expected tab-separated format (only supported format):
    - No header
    - Exactly 3 columns per row: species, trait, pair
    - contrast defaults to 1 for all rows
    - Returns {contrast: [(high_species, low_species), ...]} with pairs sorted by
      numeric pair_id where possible
    - Ignores rows with missing fields or invalid trait values
    - No taxid mapping or validation here (pure parsing)
    """
    if not trait_file_path.exists():
        logger.error("Trait file not found: %s", trait_file_path)
        raise FileNotFoundError(f"Trait file not found: {trait_file_path}")

    # Structure: contrast -> pair_id -> {'high': [species], 'low': [species]}
    by_contrast_and_pair: Dict[int, Dict[str, Dict[str, list]]] = defaultdict(
        lambda: defaultdict(lambda: {"high": [], "low": []})
    )

    def _to_int(value: str) -> Optional[int]:
        try:
            return int(str(value).strip())
        except (ValueError, TypeError):
            return None

    try:
        with open(trait_file_path, "r", encoding="utf-8-sig") as f:
            reader = csv.reader(f, delimiter="\t")
            rows_seen = 0

            def _consume_row(cells: List[str], row_num: int) -> None:
                if not cells or all(not str(c).strip() for c in cells):
                    return

                if len(cells) != 3:
                    logger.debug(
                        "Skipping malformed trait row %d (expected 3 columns): %s",
                        row_num,
                        cells,
                    )
                    return

                species = cells[0].strip()
                trait_val = _to_int(cells[1])
                pair_id = cells[2].strip()
                contrast_num = 1

                if not species or trait_val is None or not pair_id:
                    logger.debug("Skipping malformed trait row %d: %s", row_num, cells)
                    return

                if trait_val == 1:
                    by_contrast_and_pair[contrast_num][pair_id]["high"].append(species)
                elif trait_val == 0:
                    by_contrast_and_pair[contrast_num][pair_id]["low"].append(species)
                else:
                    logger.debug(
                        "Skipping row %d with non-binary trait value (%s): %s",
                        row_num,
                        trait_val,
                        cells,
                    )

            for row_num, row in enumerate(reader, start=1):
                rows_seen += 1
                _consume_row(row, row_num)

            if rows_seen == 0:
                logger.warning("Trait file is empty: %s", trait_file_path)
                return {}

        # Build final structure: contrast -> list of pairs
        contrast_to_pairs = {}
        for contrast_num, pairs_dict in by_contrast_and_pair.items():
            pairs_list = []
            # Sort by pair_id numerically to maintain order (pair 1, pair 2, pair 3)
            for pair_id in sorted(
                pairs_dict.keys(), key=lambda x: int(x) if x.isdigit() else x
            ):
                group = pairs_dict[pair_id]
                high_species = group["high"]
                low_species = group["low"]

                if high_species and low_species:
                    # Take first species from each side
                    pairs_list.append((high_species[0], low_species[0]))

                    if len(high_species) > 1 or len(low_species) > 1:
                        logger.warning(
                            "Contrast %d, pair %s has multiple species per side: "
                            "high=%s, low=%s. Using first from each.",
                            contrast_num,
                            pair_id,
                            high_species,
                            low_species,
                        )
                else:
                    logger.debug(
                        "Skipping incomplete pair %s in contrast %d: high=%s, low=%s",
                        pair_id,
                        contrast_num,
                        high_species,
                        low_species,
                    )

            contrast_to_pairs[contrast_num] = pairs_list
            logger.debug("Contrast %d: %d pairs loaded", contrast_num, len(pairs_list))

        total_pairs = sum(len(p) for p in contrast_to_pairs.values())
        logger.info(
            "Loaded %d contrasts with %d total pairs from %s",
            len(contrast_to_pairs),
            total_pairs,
            trait_file_path,
        )
        return contrast_to_pairs

    except Exception as e:
        logger.error("Failed to load trait pairs: %s", e, exc_info=True)
        return {}


# -- Function to Load Ensembl Genes --#


def load_ensembl_genes(ensembl_genes_file: Path) -> Set[str]:
    """Load Ensembl gene names from a TSV/CSV file (expects a 'gene' column)."""
    if not ensembl_genes_file.exists():
        raise FileNotFoundError(f"Ensembl genes file not found: {ensembl_genes_file}")

    with ensembl_genes_file.open(newline="") as handle:
        sample = handle.read(2048)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        # reader.fieldnames may be None (no header) — guard against that before membership test
        if not reader.fieldnames or "gene" not in reader.fieldnames:
            raise ValueError("Ensembl genes file must contain a 'gene' column")
        genes = {row["gene"].strip() for row in reader if row.get("gene")}
    if not genes:
        raise ValueError("No genes found in Ensembl genes file")
    return genes
