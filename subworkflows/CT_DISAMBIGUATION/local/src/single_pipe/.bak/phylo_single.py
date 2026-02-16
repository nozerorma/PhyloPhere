"""
Phylogenetic utilities for single-gene CAAS analysis.

This module provides functions to handle tree and alignment preprocessing,
including species-to-tax_id mapping and ensuring consistency between
tree and alignment data.
"""

import logging
from pathlib import Path
from typing import Dict, Tuple, List
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from src.phylo.species_mapping import (
    read_taxid_mapping as _read_taxid_mapping,
    match_tree_alignment_by_taxid,
    validate_taxids_in_tree
)
from src.phylo.tree_utils import load_tree
from src.utils.io_utils import read_alignment
from src.asr.tree_parser import (
    build_node_mapping as _build_node_mapping,
    get_tip_labels as _get_tip_labels
)

logger = logging.getLogger(__name__)

def preprocess_alignment_with_taxid(
    alignment: MultipleSeqAlignment,
    taxid_mapping: Dict[str, str]
) -> MultipleSeqAlignment:
    """
    Replace species names in the alignment with their corresponding tax_id values.

    Args:
        alignment: MultipleSeqAlignment object.
        taxid_mapping: Dictionary mapping species names to tax_ids.

    Returns:
        Updated MultipleSeqAlignment with species names replaced by tax_ids.
    """
    updated_records = []
    for record in alignment:
        species_name = record.id
        if species_name in taxid_mapping:
            record.id = taxid_mapping[species_name]
            updated_records.append(record)
        else:
            logger.warning(f"Species {species_name} not found in taxid mapping. Skipping.")

    if not updated_records:
        raise ValueError("No species in the alignment matched the taxid mapping.")

    return MultipleSeqAlignment(updated_records)

def load_and_preprocess_alignment(
    alignment_path: Path,
    taxid_mapping_path: Path
) -> Tuple[MultipleSeqAlignment, Dict[str, str]]:
    """
    Load alignment and preprocess it by replacing species names with tax_ids.

    Args:
        alignment_path: Path to the alignment file.
        taxid_mapping_path: Path to the taxid mapping file.

    Returns:
        Tuple containing the preprocessed alignment and the taxid mapping.
    """
    # Load taxid mapping
    taxid_mapping = _read_taxid_mapping(taxid_mapping_path)

    # Load alignment
    alignment = read_alignment(alignment_path, format='phylip-relaxed')

    # Preprocess alignment
    preprocessed_alignment = preprocess_alignment_with_taxid(alignment, taxid_mapping)

    return preprocessed_alignment, taxid_mapping

def prepare_tree_and_alignment(
    tree_path: Path,
    alignment: MultipleSeqAlignment,
    taxid_mapping: Dict[str, str]
) -> Tuple[Tree, MultipleSeqAlignment, Dict[str, str], Dict[str, str]]:
    """
    Load and preprocess the tree and alignment to ensure consistency.

    Args:
        tree_path: Path to the phylogenetic tree file.
        alignment: MultipleSeqAlignment object with species names replaced by tax_ids.
        taxid_mapping: Dictionary mapping species names to tax_ids.

    Returns:
        Tuple containing the matched tree, matched alignment, taxid_to_species mapping,
        and species_to_taxid mapping.
    """
    # Load the tree
    tree = load_tree(tree_path)

    # Match tree and alignment
    matched_tree, matched_alignment, taxid_to_species, aln_taxid_to_sp, _ = match_tree_alignment_by_taxid(
        tree,
        alignment,
        taxid_mapping
    )

    species_to_taxid = {species: taxid for taxid, species in aln_taxid_to_sp.items()}

    # Validate species mapping
    unmatched_species, missing_species = validate_taxids_in_tree(
        list(taxid_to_species.values()),
        [tip.name for tip in matched_tree.get_terminals()],
        group_name="alignment"
    )

    if unmatched_species:
        logger.warning(f"Unmatched species in alignment: {unmatched_species}")
    if missing_species:
        logger.warning(f"Species in tree but missing from alignment: {missing_species}")

    return matched_tree, matched_alignment, taxid_to_species, species_to_taxid

def load_and_preprocess_tree(
    tree_path: Path,
    alignment: MultipleSeqAlignment,
    taxid_mapping: Dict[str, str]
) -> Tuple[Tree, MultipleSeqAlignment, Dict[str, str], Dict[str, str]]:
    """
    Load and preprocess the tree to ensure consistency with the alignment.

    Args:
        tree_path: Path to the tree file.
        alignment: Preprocessed alignment with tax_ids.
        taxid_mapping: Dictionary mapping species names to tax_ids.

    Returns:
        Tuple containing the matched tree, matched alignment, taxid_to_species mapping,
        and species_to_taxid mapping.
    """
    return prepare_tree_and_alignment(tree_path, alignment, taxid_mapping)


def load_taxid_mapping(taxid_file: Path) -> Dict[str, str]:
    """Wrapper around core read_taxid_mapping for consistent imports."""
    return _read_taxid_mapping(taxid_file)


def build_tree_node_mapping(tree_file: Path, rst_file: Path = None):
    """
    Expose ASR tree parser build_node_mapping via phylo module.
    
    Args:
        tree_file: Path to Newick tree file
        rst_file: Optional path to RST file (preferred source for PAML node IDs)
    
    Returns:
        Tuple of (ordered_nodes, id_mapping)
    """
    return _build_node_mapping(tree_file=tree_file, rst_file=rst_file)


def extract_tip_labels(root_node) -> List[str]:
    """Expose get_tip_labels via phylo module for consistent access."""
    return _get_tip_labels(root_node)
