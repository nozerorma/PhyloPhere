"""
Species name to taxon ID mapping utilities.

This module provides functions to map species names (used in trait files)
to numeric taxon IDs (used in phylogenetic trees).
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict

import pandas as pd
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.SeqRecord import SeqRecord

from src.phylo.tree_utils import prune_tree

logger = logging.getLogger(__name__)

# Per-process deduplication: log each taxid conflict only once per worker to
# avoid flooding stderr when processing thousands of genes with the same conflict.
_WARNED_CONFLICT_TAXIDS: Set[str] = set()


def read_taxid_mapping(taxid_file: Path) -> Dict[str, str]:
    """
    Read taxid mapping file and return species name → taxon ID dictionary.

    The mapping file should have columns: tax_id, species, family, rank, name_class
    Example format:
        tax_id  species                 family          rank     name_class
        9606    Homo_sapiens           Hominidae       species  scientific name
        9598    Pan_troglodytes        Hominidae       species  scientific name

    Args:
        taxid_file: Path to taxid mapping file (tab-separated) - Path object or string

    Returns:
        Dictionary mapping species names to taxon IDs (both as strings)
        Example: {'Homo_sapiens': '9606', 'Pan_troglodytes': '9598'}

    Raises:
        FileNotFoundError: If taxid file doesn't exist
        ValueError: If required columns are missing
    """
    # Convert to Path if input is string
    if isinstance(taxid_file, str):
        taxid_file = Path(taxid_file)

    if not taxid_file.exists():
        raise FileNotFoundError(f"Taxid mapping file not found: {taxid_file}")

    logger.info(f"Reading taxid mapping from {taxid_file}")

    try:
        df = pd.read_csv(taxid_file, sep="\t", dtype=str)
    except Exception as e:
        raise ValueError(f"Error reading taxid file: {e}")

    # Validate required columns
    required_cols = ["tax_id", "species"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Taxid file missing required columns: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    # Create mapping dictionary
    mapping = {}
    for _, row in df.iterrows():
        species_name = str(row["species"]).strip()
        taxon_id = str(row["tax_id"]).strip()

        if species_name and taxon_id:
            mapping[species_name] = taxon_id

    logger.info(f"Loaded {len(mapping)} species → taxon ID mappings")
    logger.debug(f"Sample mappings: {dict(list(mapping.items())[:5])}")

    return mapping


def validate_taxids_in_tree(
    taxids: List[str], tree_tips: List[str], group_name: str = "species"
) -> Tuple[List[str], List[str]]:
    """
    Validate that taxon IDs exist in the phylogenetic tree.

    Args:
        taxids: List of taxon IDs to validate
        tree_tips: List of tip labels from the tree
        group_name: Name of the group for logging (e.g., "TOP", "BOTTOM")

    Returns:
        Tuple of (found_taxids, missing_taxids)

    Example:
        >>> tree_tips = ['9606', '9598', '9544']
        >>> taxids = ['9606', '9598', '9999']
    >>> found, missing = validate_taxids_in_tree(taxids, tree_tips, "TOP")
        >>> print(f"Found: {found}, Missing: {missing}")
        Found: ['9606', '9598'], Missing: ['9999']
    """
    tree_tips_set = set(tree_tips)

    found = [tid for tid in taxids if tid in tree_tips_set]
    missing = [tid for tid in taxids if tid not in tree_tips_set]

    logger.info(
        f"Validated {group_name} taxon IDs: "
        f"{len(found)}/{len(taxids)} found in tree"
    )

    if missing:
        logger.warning(f"Missing {group_name} taxon IDs not in tree: {missing}")

    return found, missing


def match_tree_alignment_by_taxid(
    tree: Tree,
    alignment: MultipleSeqAlignment,
    tax_mapping: Dict[str, str],
) -> Tuple[
    Tree,
    MultipleSeqAlignment,
    Dict[str, str],
    Dict[str, str],
    Dict[str, Dict[str, str]],
]:
    """Match tree and alignment species using tax_id mapping.

    Handles:
    - Synonyms (multiple names → same tax_id)
    - Spelling variations
    - Subspecies/species differences
    - Taxonomy conflicts (creates synthetic tax_ids for phylogenetically distinct
      species sharing same NCBI tax_id)

    Args:
        tree: Input phylogenetic tree
        alignment: Input multiple sequence alignment
        tax_mapping: Dict mapping species_name → tax_id

    Returns:
        Tuple of:
        - Matched tree (pruned, terminals relabeled with tax_ids)
        - Matched alignment (filtered, sequences relabeled with tax_ids)
        - Dict mapping tax_id → original tree name
        - Dict mapping tax_id → original alignment name
        - Dict mapping species_name → {original_taxid, synthetic_taxid, reason}
          (for species assigned synthetic tax_ids due to conflicts)

    Raises:
        ValueError: If no species match between tree and alignment
    """

    logger.info("Matching tree and alignment species using tax_id mapping...")

    tree_species = {tip.name for tip in tree.get_terminals()}
    aln_species = {rec.id for rec in alignment}

    tree_sp_to_taxid: Dict[str, str] = {}
    tree_taxid_to_sp: Dict[str, str] = {}
    tree_unmatched: List[str] = []

    for sp in tree_species:
        if sp in tax_mapping:
            taxid = tax_mapping[sp]
            tree_sp_to_taxid[sp] = taxid
            tree_taxid_to_sp[taxid] = sp
        else:
            tree_unmatched.append(sp)
            logger.warning(f"Tree species not in tax_id mapping: {sp}")

    logger.info(
        "Tree: %d species mapped to tax_ids, %d unmatched",
        len(tree_sp_to_taxid),
        len(tree_unmatched),
    )

    aln_sp_to_taxid: Dict[str, str] = {}
    aln_taxid_to_sp: Dict[str, str] = {}
    aln_taxid_duplicates: Dict[str, List[str]] = defaultdict(list)
    aln_unmatched: List[str] = []

    for sp in aln_species:
        if sp in tax_mapping:
            taxid = tax_mapping[sp]
            aln_sp_to_taxid[sp] = taxid
            aln_taxid_duplicates[taxid].append(sp)
            aln_taxid_to_sp.setdefault(taxid, sp)
        else:
            aln_unmatched.append(sp)
            logger.warning(f"Alignment species not in tax_id mapping: {sp}")

    synthetic_taxids: Dict[str, Dict[str, str]] = {}
    all_existing_taxids: Set[str] = set(tree_sp_to_taxid.values()) | set(
        aln_sp_to_taxid.values()
    )

    for taxid, species_list in aln_taxid_duplicates.items():
        if len(species_list) <= 1:
            continue

        species_list = sorted(species_list)
        # Only emit the full warning block once per taxid conflict per worker process.
        # In a multi-gene run this conflict recurs for every gene sharing these species;
        # repeated logging floods stderr and can cause SLURM to kill the job.
        first_occurrence = taxid not in _WARNED_CONFLICT_TAXIDS
        _WARNED_CONFLICT_TAXIDS.add(taxid)

        if first_occurrence:
            logger.warning(
                "TAXONOMY CONFLICT — tax_id %s shared by: %s. "
                "Assigning synthetic tax_ids to duplicates. "
                "This message is shown once per worker; further occurrences suppressed.",
                taxid,
                ", ".join(species_list),
            )
        else:
            logger.debug(
                "Repeated taxonomy conflict for tax_id %s (%s) — synthetic tax_id reassignment applied silently.",
                taxid,
                ", ".join(species_list),
            )

        kept = species_list[0]
        duplicates = species_list[1:]

        for i, dup_sp in enumerate(duplicates, start=1):
            synthetic_taxid = str(int(taxid) + i)
            max_attempts = 1000
            attempts = 0
            while synthetic_taxid in all_existing_taxids and attempts < max_attempts:
                synthetic_taxid = str(int(synthetic_taxid) + 1)
                attempts += 1

            if attempts >= max_attempts:
                raise RuntimeError(
                    f"Could not find unused synthetic tax_id for {dup_sp} (tried {max_attempts} IDs)"
                )

            synthetic_taxids[dup_sp] = {
                "original_taxid": taxid,
                "synthetic_taxid": synthetic_taxid,
                "reason": f"Duplicate of {kept}",
            }

            all_existing_taxids.add(synthetic_taxid)
            aln_sp_to_taxid[dup_sp] = synthetic_taxid
            aln_taxid_to_sp[synthetic_taxid] = dup_sp

            if dup_sp in tree_sp_to_taxid:
                old_tree_taxid = tree_sp_to_taxid[dup_sp]
                tree_sp_to_taxid[dup_sp] = synthetic_taxid

                other_species = [
                    sp
                    for sp, tid in tree_sp_to_taxid.items()
                    if tid == old_tree_taxid and sp != dup_sp
                ]

                if not other_species:
                    tree_taxid_to_sp.pop(old_tree_taxid, None)
                else:
                    tree_taxid_to_sp[old_tree_taxid] = other_species[0]

                tree_taxid_to_sp[synthetic_taxid] = dup_sp
                logger.debug(
                    "Updated TREE mapping for '%s': %s -> %s",
                    dup_sp,
                    old_tree_taxid,
                    synthetic_taxid,
                )

            logger.debug(
                "Synthetic tax_id %s assigned to '%s' (original: %s, kept: %s)",
                synthetic_taxid,
                dup_sp,
                taxid,
                kept,
            )

        logger.debug("Kept '%s' with original tax_id %s", kept, taxid)

    if synthetic_taxids:
        logger.debug(
            "SUMMARY: %d species assigned synthetic tax_ids due to conflicts: %s",
            len(synthetic_taxids),
            ", ".join(
                f"{sp}={info['synthetic_taxid']}" for sp, info in synthetic_taxids.items()
            ),
        )

    logger.info(
        "Alignment: %d species mapped to tax_ids, %d unmatched",
        len(aln_sp_to_taxid),
        len(aln_unmatched),
    )

    tree_taxids = set(tree_sp_to_taxid.values())
    aln_taxids = set(aln_sp_to_taxid.values())
    common_taxids = tree_taxids & aln_taxids

    if not common_taxids:
        raise ValueError(
            "No common species between tree and alignment!\n"
            f"Tree: {len(tree_species)} species ({len(tree_taxids)} with tax_ids)\n"
            f"Alignment: {len(aln_species)} species ({len(aln_taxids)} with tax_ids)\n"
            "Consider checking tax_id mapping file."
        )

    logger.info(
        "Found %d common tax_ids between tree and alignment", len(common_taxids)
    )

    tree_only_taxids = tree_taxids - aln_taxids
    aln_only_taxids = aln_taxids - tree_taxids

    if tree_only_taxids:
        logger.info("Tax_ids in tree but not alignment: %d", len(tree_only_taxids))
        if len(tree_only_taxids) <= 10:
            for taxid in list(tree_only_taxids)[:10]:
                sp = tree_taxid_to_sp.get(taxid, "?")
                logger.debug("  %s (%s)", taxid, sp)

    if aln_only_taxids:
        logger.info("Tax_ids in alignment but not tree: %d", len(aln_only_taxids))
        if len(aln_only_taxids) <= 10:
            for taxid in list(aln_only_taxids)[:10]:
                sp = aln_taxid_to_sp.get(taxid, "?")
                logger.debug("  %s (%s)", taxid, sp)

    species_to_keep_in_tree = [
        tree_taxid_to_sp[taxid] for taxid in common_taxids if taxid in tree_taxid_to_sp
    ]

    pruned_tree = prune_tree(tree, species_to_keep_in_tree)
    logger.info("Pruned tree to %d species", len(species_to_keep_in_tree))

    for tip in pruned_tree.get_terminals():
        original_name = tip.name
        if original_name in tree_sp_to_taxid:
            tip.name = tree_sp_to_taxid[original_name]
        else:
            logger.error(
                "Tree tip %s not found in mapping (should not happen)", original_name
            )

    filtered_records = []
    seen_taxids: Set[str] = set()

    for rec in alignment:
        original_id = rec.id
        if original_id not in aln_sp_to_taxid:
            continue

        taxid = aln_sp_to_taxid[original_id]
        if taxid not in common_taxids:
            continue

        if taxid in seen_taxids:
            logger.error(
                "🚨 DUPLICATE tax_id %s found for species '%s'!", taxid, original_id
            )
            continue

        new_rec = SeqRecord(
            seq=rec.seq,
            id=taxid,
            name=taxid,
            description=f"[original: {original_id}]",
        )
        filtered_records.append(new_rec)
        seen_taxids.add(taxid)

    if not filtered_records:
        raise ValueError("No alignment sequences remain after filtering!")

    filtered_alignment = MultipleSeqAlignment(filtered_records)
    logger.info("Filtered alignment to %d sequences", len(filtered_records))

    tree_terminal_count = len(pruned_tree.get_terminals())
    aln_seq_count = len(filtered_alignment)

    if tree_terminal_count != aln_seq_count:
        logger.warning(
            "Count mismatch: tree has %d terminals, alignment has %d sequences",
            tree_terminal_count,
            aln_seq_count,
        )

    logger.info("✓ Tree and alignment matched successfully via tax_ids")

    return (
        pruned_tree,
        filtered_alignment,
        tree_taxid_to_sp,
        aln_taxid_to_sp,
        synthetic_taxids,
    )
