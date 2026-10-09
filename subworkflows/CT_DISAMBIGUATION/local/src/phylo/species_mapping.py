# species_mapping.py — Read the species-to-tax_id table and match a tree and an alignment through it.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/phylo/

"""
Species name to taxon ID mapping utilities.

Reads the species name -> tax_id table and uses it to match the species of a tree and of an alignment: both are pruned
or filtered to their common tax_ids and relabeled with them.

Imported by: src/asr/asr_single.py
Inputs: a tab-separated tax_id table (columns tax_id and species), a Bio.Phylo tree, a Bio.Align alignment
Outputs: the matched tree and alignment and the tax_id -> original name maps (in memory)
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

# Per-process deduplication: each unmatched species is warned about once per worker, to avoid
# flooding stderr when thousands of genes share the same species.
_WARNED_TREE_UNMATCHED_SPECIES: Set[str] = set()
_WARNED_ALIGNMENT_UNMATCHED_SPECIES: Set[str] = set()


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
    # Accept a plain string path
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

    # Species name -> tax_id
    mapping = {}
    for _, row in df.iterrows():
        species_name = str(row["species"]).strip()
        taxon_id = str(row["tax_id"]).strip()

        if species_name and taxon_id:
            mapping[species_name] = taxon_id

    logger.info(f"Loaded {len(mapping)} species → taxon ID mappings")
    logger.debug(f"Sample mappings: {dict(list(mapping.items())[:5])}")

    return mapping


def match_tree_alignment_by_taxid(
    tree: Tree,
    alignment: MultipleSeqAlignment,
    tax_mapping: Dict[str, str],
) -> Tuple[
    Tree,
    MultipleSeqAlignment,
    Dict[str, str],
    Dict[str, str],
]:
    """Match tree and alignment species using tax_id mapping.

    Species are matched by exact name through `tax_mapping`; tree and alignment names that are not in it are
    dropped. Each species must have its own tax_id (the map of NAME_CURATION guarantees it), so each sequence has a
    distinct label; two species sharing a tax_id raise a ValueError.

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

    Raises:
        ValueError: If no species match between tree and alignment, or if two species share a tax_id
    """

    logger.info("Matching tree and alignment species using tax_id mapping...")

    # Sorted so every mapping below is independent of set iteration order, i.e. of PYTHONHASHSEED.
    tree_species = sorted({tip.name for tip in tree.get_terminals()})
    aln_species = sorted({rec.id for rec in alignment})

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

    if tree_unmatched:
        # Report each missing species once, not for every gene.
        unseen_tree_unmatched = sorted(
            [sp for sp in tree_unmatched if sp not in _WARNED_TREE_UNMATCHED_SPECIES]
        )
        _WARNED_TREE_UNMATCHED_SPECIES.update(tree_unmatched)

        if unseen_tree_unmatched:
            preview = ", ".join(unseen_tree_unmatched[:10])
            suffix = "" if len(unseen_tree_unmatched) <= 10 else ", ..."
            logger.warning(
                "Pruning %d tree species missing in tax_id mapping (showing %d new): %s%s",
                len(tree_unmatched),
                len(unseen_tree_unmatched),
                preview,
                suffix,
            )
        else:
            logger.debug(
                "Pruning %d tree species missing in tax_id mapping (all previously reported)",
                len(tree_unmatched),
            )

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

    if aln_unmatched:
        unseen_aln_unmatched = sorted(
            [
                sp
                for sp in aln_unmatched
                if sp not in _WARNED_ALIGNMENT_UNMATCHED_SPECIES
            ]
        )
        _WARNED_ALIGNMENT_UNMATCHED_SPECIES.update(aln_unmatched)

        if unseen_aln_unmatched:
            preview = ", ".join(unseen_aln_unmatched[:10])
            suffix = "" if len(unseen_aln_unmatched) <= 10 else ", ..."
            logger.warning(
                "Ignoring %d alignment species missing in tax_id mapping (showing %d new): %s%s",
                len(aln_unmatched),
                len(unseen_aln_unmatched),
                preview,
                suffix,
            )
        else:
            logger.debug(
                "Ignoring %d alignment species missing in tax_id mapping (all previously reported)",
                len(aln_unmatched),
            )

    # Every species has its own tax_id: NAME_CURATION writes a map with one unique tax_id per
    # species. A tax_id shared by several tree or alignment species means the map is not that one.
    shared = {tid: sorted(sps) for tid, sps in aln_taxid_duplicates.items() if len(sps) > 1}
    tree_by_taxid: Dict[str, List[str]] = defaultdict(list)
    for sp, tid in tree_sp_to_taxid.items():
        tree_by_taxid[tid].append(sp)
    shared.update({tid: sorted(sps) for tid, sps in tree_by_taxid.items() if len(sps) > 1})
    if shared:
        raise ValueError(
            "tax_id shared by several species: "
            + "; ".join(f"{tid}: {', '.join(sps)}" for tid, sps in sorted(shared.items()))
            + ". Use the tax_id map written by NAME_CURATION (name_curation/species_taxid_map.tsv), "
            "which gives every species a unique tax_id."
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

    # Drop the tree tips without a tax_id mapping before the common-species prune,
    # so that the drop is explicit and logged.
    if tree_unmatched:
        tree = prune_tree(tree, sorted(tree_sp_to_taxid.keys()))

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
    )
