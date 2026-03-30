#!/usr/bin/env python3
"""
ASR (Ancestral State Reconstruction) Module for Single Gene Analysis

Handles all ancestral state reconstruction operations including:
- Alignment loading and parsing
- Tree loading and species matching
- PAML ASR execution
- Posterior probability parsing
- Phylogenetic context building

Author: Refactored from test_nutm2a_real_caas.py
Date: 2025-11-24
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Mapping
from dataclasses import dataclass
import logging
import json
from io import StringIO
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo

# Add src to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

# Core imports
from src.utils.io_utils import read_alignment
from src.asr.reconstruct import ASRReconstructor, ASRConfig
from src.asr.posterior import parse_paml_rst, parse_paml_rst_node_level
from src.asr.tree_parser import (
    get_tip_labels,
    parse_newick,
    assign_postorder_ids,
    get_node_order,
)
from src.phylo.species_mapping import match_tree_alignment_by_taxid, read_taxid_mapping
from src.phylo.tree_utils import load_tree

logger = logging.getLogger(__name__)


@dataclass
class SingleGeneASRConfig:
    """Configuration for single gene ASR pipeline execution."""

    alignment_path: Path
    tree_path: Path
    taxid_path: Optional[Path] = None
    output_dir: Path = Path("./output")
    model: str = "lg"
    posterior_threshold: float = 0.7
    threads: int = 1
    run_diagnostics: bool = False


@dataclass
class AlignmentData:
    """Container for alignment and mapping data."""

    alignment: MultipleSeqAlignment  # BioPython alignment object
    # Mapping orientations kept consistent:
    #   species_to_taxid: species name -> taxid (for trait/metadata files)
    #   taxid_to_species: taxid -> species name (for alignment/tree lookups)
    species_to_taxid: Dict[str, str]
    taxid_to_species: Dict[str, str]
    seq_by_id: Dict[str, str]
    seq_by_species: Dict[str, str]
    gene_name: str


@dataclass
class TreeData:
    """Container for tree and phylogenetic data."""

    tree: Phylo.BaseTree.Tree  # BioPython tree object
    tip_set: set
    nodes: list
    root: object
    taxid_mapping: Dict[str, str]
    node_mapping: Optional[Mapping[int, object]] = None


@dataclass
class ASRResults:
    """Container for ASR execution results."""

    gene: str
    posteriors_site: Dict[int, Dict[str, float]]
    posteriors_node: Optional[Dict[int, Dict[int, Dict[str, float]]]] = None
    rst_file: Optional[Path] = None
    tree_file: Optional[Path] = None
    taxid_to_species: Optional[Dict[str, str]] = None
    node_id_map: Optional[Dict[int, int]] = None


def _write_node_id_map(
    node_id_map: Optional[Dict[int, int]], tree_file: Optional[Path]
) -> None:
    """Write node-id remapping JSON next to the provided tree."""
    if not node_id_map or tree_file is None:
        return
    map_path = tree_file.with_suffix(tree_file.suffix + ".node_map.json")
    try:
        with map_path.open("w", encoding="utf-8") as handle:
            json.dump(node_id_map, handle, indent=2)
    except Exception as exc:
        logger.warning(f"Could not write node-id map to {map_path}: {exc}")


def load_alignment_and_mappings(
    alignment_path: Path, taxid_path: Optional[Path] = None, gene_name: str = "unknown"
) -> AlignmentData:
    """
    Load alignment and create lookup mappings.

    Args:
        alignment_path: Path to the input alignment file
        taxid_path: Optional path to taxid-species mapping file
        gene_name: Name of the gene for logging

    Returns:
        AlignmentData object with loaded alignment and mappings
    """
    logger.debug(f"Loading alignment for {gene_name}: {alignment_path}")

    if not alignment_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {alignment_path}")

    # Load alignment
    alignment: MultipleSeqAlignment = read_alignment(alignment_path, format="auto")
    logger.debug(
        f"Alignment loaded: {len(alignment)} sequences, {alignment.get_alignment_length()} positions"
    )

    # Build sequence lookups
    seq_by_id = {}
    seq_by_species = {}
    for rec in alignment:
        rec_id = str(rec.id)
        seq = str(rec.seq)
        seq_by_id[rec_id] = seq
        seq_by_species[rec_id] = seq  # Fallback mapping

    species_to_taxid: Dict[str, str] = {}
    taxid_to_species: Dict[str, str] = {}
    if taxid_path and taxid_path.exists():
        logger.debug("Loading TaxID to species mappings")
        species_to_taxid = read_taxid_mapping(taxid_path)  # species -> taxid
        taxid_to_species = {
            taxid: species for species, taxid in species_to_taxid.items()
        }

        # Update species lookups using mapping
        for species, taxid in species_to_taxid.items():
            if species in seq_by_id:
                seq_by_species[species] = seq_by_id[species]
            if taxid in seq_by_id:
                seq_by_species[species] = seq_by_id[taxid]

        logger.debug(f"Loaded {len(species_to_taxid)} TaxID mappings")

    return AlignmentData(
        alignment=alignment,
        species_to_taxid=species_to_taxid,
        taxid_to_species=taxid_to_species,
        seq_by_id=seq_by_id,
        seq_by_species=seq_by_species,
        gene_name=gene_name,
    )


def load_and_match_tree(
    tree_path: Path, alignment_data: AlignmentData, taxid_path: Optional[Path] = None
) -> TreeData:
    """
    Load phylogenetic tree and match to alignment species.

    Args:
        tree_path: Path to newick/nexus tree file
        alignment_data: Alignment data with mappings
        taxid_path: Optional TaxID mapping path

    Returns:
        TreeData object with matched tree and phylogenetic context
    """
    logger.debug(f"Loading tree: {tree_path}")

    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    tree = load_tree(tree_path)
    logger.debug(f"Tree loaded: {len(tree.get_terminals())} tips")

    # Build taxid mapping
    taxid_mapping = {}
    if taxid_path and taxid_path.exists():
        taxid_mapping = read_taxid_mapping(taxid_path)
        logger.debug("Attempting tree-alignment matching via TaxID...")

        # Match tree to alignment
        (
            matched_tree,
            matched_alignment,
            tree_taxid_to_sp,
            aln_taxid_to_sp,
            synthetic_taxids,
        ) = match_tree_alignment_by_taxid(tree, alignment_data.alignment, taxid_mapping)

        tree = matched_tree

        # Update alignment and mappings after relabeling to taxids
        alignment_data.alignment = matched_alignment
        alignment_data.taxid_to_species = aln_taxid_to_sp  # taxid -> species
        alignment_data.species_to_taxid = {
            species: taxid for taxid, species in aln_taxid_to_sp.items()
        }

        # Rebuild sequence lookups with new IDs (taxids)
        seq_by_id = {}
        seq_by_species = {}
        for rec in alignment_data.alignment:
            rec_id = str(rec.id)  # taxid after matching
            seq = str(rec.seq)
            seq_by_id[rec_id] = seq

            species_name = alignment_data.taxid_to_species.get(rec_id, rec_id)
            seq_by_species[species_name] = seq

        alignment_data.seq_by_id = seq_by_id
        alignment_data.seq_by_species = seq_by_species

        taxid_mapping = aln_taxid_to_sp

        logger.debug(
            f"After matching: {len(tree.get_terminals())} tips, "
            f"{len(synthetic_taxids)} synthetic taxids created"
        )

    # Build tree structure data (use matched tree tips as valid taxids)
    taxid_mapping_dict = taxid_mapping

    # Convert matched tree to Newick string and build node mapping so MRCA lookups work
    newick_buffer = StringIO()
    from Bio import Phylo

    # Ensure tree is unrooted to avoid PAML/rooting issues
    try:
        logger.debug(f"Tree rooted flag before ASR: {getattr(tree, 'rooted', None)}")
        if getattr(tree, "rooted", False):
            logger.debug("Input tree is rooted; unrooting before ASR")
            tree.rooted = False
            if hasattr(tree, "root"):
                tree.root = None
    except Exception:
        logger.debug("Could not determine rooting; proceeding as-is")

    Phylo.write(tree, newick_buffer, "newick")
    tree_newick = newick_buffer.getvalue().strip()

    parsed_root = parse_newick(tree_newick)
    assign_postorder_ids(parsed_root, start_id=0)
    ordered_nodes = get_node_order(parsed_root)
    id_mapping = {n.node_id: n for n in ordered_nodes if n.node_id is not None}
    tree_tip_set = set(get_tip_labels(parsed_root))

    return TreeData(
        tree=tree,
        tip_set=tree_tip_set,
        nodes=ordered_nodes,
        root=parsed_root,
        taxid_mapping=taxid_mapping_dict,
        node_mapping=id_mapping,
    )


def run_asr_pipeline(
    gene: str,
    config: SingleGeneASRConfig,
    skip_if_exists: bool = True,
    alignment_data=None,
    tree_data=None,
) -> ASRResults:
    """
    Execute complete ASR pipeline for a gene.

    Args:
        gene: Gene name
        config: ASR configuration
        skip_if_exists: Skip PAML execution if RST and RST1 files exist

    Returns:
        ASRResults object with all ASR data
    """
    logger.debug(f"Starting ASR pipeline for {gene}")

    config.output_dir.mkdir(parents=True, exist_ok=True)
    # PAML creates output in asr_{gene} subdirectory
    gene_asr_dir = config.output_dir / f"asr_{gene}"
    rst_file = gene_asr_dir / "rst"
    rst1_file = gene_asr_dir / "rst1"

    # Check for existing results (both rst and rst1 must be present for complete run)
    if skip_if_exists and rst_file.exists() and rst1_file.exists():
        logger.debug(f"Found existing RST/RST1 files, using cached results: {rst_file}")
    else:
        logger.debug(f"Running PAML ASR for {gene}...")

        # Load alignment and tree unless provided
        if alignment_data is None:
            alignment_data = load_alignment_and_mappings(
                config.alignment_path, config.taxid_path, gene
            )

        if tree_data is None:
            tree_data = load_and_match_tree(
                config.tree_path, alignment_data, config.taxid_path
            )

        # Configure and run PAML (exact same params as working test)
        asr_config = ASRConfig(
            model=config.model,
            compute_asr=True,
            output_dir=config.output_dir,
            fix_blength=2,
            threads=config.threads,
        )

        reconstructor = ASRReconstructor(asr_config)

        rst_file = reconstructor.reconstruct_gene(
            gene=gene,
            alignment=alignment_data.alignment,
            tree=tree_data.tree,
            output_dir=config.output_dir,
        )

        if not rst_file.exists():
            raise RuntimeError(f"PAML failed: RST file not created at {rst_file}")

        logger.debug("PAML ASR completed")

    # Parse posteriors
    logger.debug("Parsing ASR posteriors...")

    # Determine if node-level data is available
    # Prefer standard tree_paml.nwk (what reconstruct writes); allow gene-prefixed legacy name
    tree_paml_file = rst_file.parent / "tree_paml.nwk"
    if not tree_paml_file.exists():
        legacy_tree = rst_file.parent / f"{gene.lower()}_tree_paml.nwk"
        if legacy_tree.exists():
            tree_paml_file = legacy_tree
            logger.debug(f"Using legacy gene-prefixed PAML tree: {tree_paml_file}")
        else:
            logger.debug(
                "No tree_paml.nwk found; posterior parsing may lack node remap"
            )

    node_id_map = None
    if tree_paml_file.exists():
        parsed = parse_paml_rst_node_level(
            rst_file,
            tree_paml_file,
            threshold=config.posterior_threshold,
            return_mapping=True,
        )
        if isinstance(parsed, tuple):
            posteriors_node, node_id_map = parsed
        else:
            posteriors_node = parsed
            node_id_map = None
    else:
        posteriors_node = None
        node_id_map = None

    # Always parse site-level for backward compatibility
    posteriors_site = parse_paml_rst(rst_file)

    logger.debug(
        "Posteriors parsed "
        f"(node-level: {len(posteriors_node) if isinstance(posteriors_node, dict) else 0} nodes, "
        f"site-level: {len(posteriors_site)} sites)"
    )

    return ASRResults(
        gene=gene,
        posteriors_site=posteriors_site,
        posteriors_node=posteriors_node if isinstance(posteriors_node, dict) else None,
        rst_file=rst_file,
        tree_file=tree_paml_file if tree_paml_file.exists() else None,
        taxid_to_species=None,  # Will be set by caller if needed
        node_id_map=node_id_map if isinstance(node_id_map, dict) else None,
    )


def load_precomputed_asr(
    gene: str, config: SingleGeneASRConfig, alignment_data: AlignmentData
) -> ASRResults:
    """
    Load pre-computed ASR results without running PAML.

    Args:
        gene: Gene name
        config: ASR configuration (with RST file paths)
        alignment_data: Alignment data with mappings

    Returns:
        ASRResults object with loaded data
    """
    logger.debug(f"Loading pre-computed ASR results for {gene}")

    # Try canonical location first: config.output_dir / asr_{GENE} / rst
    rst_file = config.output_dir / f"asr_{gene}" / "rst"
    tree_file = config.output_dir / f"asr_{gene}" / "tree_paml.nwk"

    # Fall back to legacy location: config.output_dir / {GENE} / asr_{GENE} / rst
    if not rst_file.exists():
        legacy_rst_file = config.output_dir / gene / f"asr_{gene}" / "rst"
        if legacy_rst_file.exists():
            rst_file = legacy_rst_file
            tree_file = config.output_dir / gene / f"asr_{gene}" / "tree_paml.nwk"

    if not rst_file.exists():
        raise FileNotFoundError(
            f"Pre-computed RST file not found at:\n"
            f"  Canonical: {config.output_dir / f'asr_{gene}' / 'rst'}\n"
            f"  Legacy: {config.output_dir / gene / f'asr_{gene}' / 'rst'}"
        )

    logger.debug(f"Found RST file at: {rst_file}")

    # Parse posteriors
    posteriors_site = parse_paml_rst(rst_file)

    posteriors_node = None
    node_id_map = None
    if tree_file.exists():
        parsed = parse_paml_rst_node_level(
            rst_file,
            tree_file,
            threshold=config.posterior_threshold,
            return_mapping=True,
        )
        if isinstance(parsed, tuple):
            posteriors_node, node_id_map = parsed
        else:
            posteriors_node = parsed

    logger.debug("Pre-computed ASR results loaded")

    _write_node_id_map(
        node_id_map if isinstance(node_id_map, dict) else None,
        tree_file if tree_file and tree_file.exists() else None,
    )

    return ASRResults(
        gene=gene,
        posteriors_site=posteriors_site,
        posteriors_node=posteriors_node if isinstance(posteriors_node, dict) else None,
        rst_file=rst_file,
        tree_file=tree_file if tree_file.exists() else None,
        taxid_to_species=alignment_data.taxid_to_species,
        node_id_map=node_id_map if isinstance(node_id_map, dict) else None,
    )


def validate_asr_inputs(config: SingleGeneASRConfig) -> bool:
    """
    Validate that all required inputs for ASR are present.

    Args:
        config: ASR configuration

    Returns:
        True if all inputs are valid

    Raises:
        FileNotFoundError: If required files don't exist
        ValueError: If configuration is invalid
    """
    required_files = [
        ("alignment_path", config.alignment_path),
        ("tree_path", config.tree_path),
    ]

    if config.taxid_path:
        required_files.append(("taxid_path", config.taxid_path))

    for name, path in required_files:
        if not path.exists():
            raise FileNotFoundError(f"Required {name} not found: {path}")

    if config.posterior_threshold < 0 or config.posterior_threshold > 1:
        raise ValueError(
            f"Invalid posterior threshold: {config.posterior_threshold} (must be 0-1)"
        )

    return True
