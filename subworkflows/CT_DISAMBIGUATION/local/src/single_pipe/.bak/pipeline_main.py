#!/usr/bin/env python3
"""
Main Pipeline Coordinator for Single Gene CAAS Analysis

Orchestrates the modular CAAS analysis pipeline:
- ASR (Ancestral State Reconstruction) - optional
- Biochemical analysis and convergence detection
- Quality scoring and reporting

Supports both full pipeline execution and pre-computed ASR loading.

Author: Refactored from test_nutm2a_real_caas.py
Date: 2025-11-24
"""

import sys
import argparse
import json
import logging
import random
import csv
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple, cast, Set

# Add src to path for local imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))
from src.utils.logger import configure_logging

# Import pipeline modules
try:
    # Try relative imports first (when running as part of package)
    from .asr_single import (
        SingleGeneASRConfig,
        run_asr_pipeline,
        load_precomputed_asr,
        validate_asr_inputs,
        load_alignment_and_mappings,
        load_and_match_tree
    )

    from .biochem_single import (
        analyze_gene_biochemistry,
        summarize_biochemical_results
    )

    from .phylo_single import (
        load_taxid_mapping,
        build_tree_node_mapping,
        extract_tip_labels
    )

    # Import from report_single module
    from .report_single import (
        create_complete_gene_report
    )

except ImportError:
    # Fallback to absolute imports (when running script directly)
    from asr_single import (
        SingleGeneASRConfig,
        run_asr_pipeline,
        load_precomputed_asr,
        validate_asr_inputs,
        load_alignment_and_mappings,
        load_and_match_tree
    )

    from biochem_single import (
        analyze_gene_biochemistry,
        summarize_biochemical_results
    )

    from phylo_single import (
        load_taxid_mapping,
        build_tree_node_mapping,
        extract_tip_labels
    )

    # Import from report_single module
    from report_single import (
        create_complete_gene_report
    )


logger = logging.getLogger(__name__)


def load_ensembl_genes(ensembl_genes_file: Path) -> Set[str]:
    """Load Ensembl gene names from a TSV/CSV file (expects a 'gene' column)."""
    if not ensembl_genes_file.exists():
        raise FileNotFoundError(f"Ensembl genes file not found: {ensembl_genes_file}")

    with ensembl_genes_file.open(newline='') as handle:
        sample = handle.read(2048)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None or 'gene' not in reader.fieldnames:
            raise ValueError("Ensembl genes file must contain a 'gene' column")
        genes = {row['gene'].strip() for row in reader if row.get('gene')}
    if not genes:
        raise ValueError("No genes found in Ensembl genes file")
    return genes


def load_caas_positions_for_gene(caas_metadata_path: Path, gene: str) -> List[int]:
    """
    Load all CAAS positions available for a specific gene from metadata file.

    Args:
        caas_metadata_path: Path to CAAS metadata (.output) file
        gene: Gene name (e.g., "NUTM2A")

    Returns:
        List of zero-based MSA positions for the gene
    """
    from src.data.loaders import load_caas_metadata_df

    logger.info(f"Loading CAAS positions for {gene} from {caas_metadata_path}")
    df = load_caas_metadata_df(caas_metadata_path, gene)
    logger.debug(f"Metadata has {len(df)} entries for {gene}")

    positions = []

    for gene_pos in df['GenePos']:
        # Parse NUTM2A_85 -> pos=85
        token = str(gene_pos).strip()
        if "_" in token:
            gene_name, pos_str = token.rsplit("_", 1)
            if gene_name == gene:
                try:
                    pos = int(pos_str)
                    positions.append(pos)
                except ValueError:
                    logger.debug(f"Could not parse position from {token}")
                    continue

    logger.info(f"Found {len(positions)} positions for {gene}: {sorted(positions[:5])}...")
    return positions


def find_gene_alignment(
    alignment_dir: Optional[Path],
    gene: str,
    ensembl_genes: Optional[Set[str]] = None,
) -> Path:
    """
    Find alignment file for a specific gene.

    Args:
        alignment_dir: Directory containing alignment files
        gene: Gene name

    Returns:
        Path to alignment file

    Raises:
        FileNotFoundError: If no suitable alignment file is found
    """
    if ensembl_genes is not None and gene not in ensembl_genes:
        raise FileNotFoundError(
            f"Gene '{gene}' not present in provided Ensembl list; refusing to match alignments"
        )

    if alignment_dir and alignment_dir.exists():
        # Require exact prefix match before the first dot to avoid partials (e.g., MHS1 vs MHS12)
        candidates = alignment_dir.glob("**/*.phy")
        for path in candidates:
            prefix = path.name.split(".", 1)[0]
            if prefix == gene:
                return path

    raise FileNotFoundError(f"No alignment file found for gene {gene}")


def _select_debug_results(results: List[Any], count: int = 3) -> List[Any]:
    """Choose up to `count` results with valid positions for deep debugging."""
    valid = [r for r in results if getattr(r, 'position', None) is not None]
    if not valid:
        return []
    if len(valid) <= count:
        return valid
    rng = random.Random()
    return rng.sample(valid, count)


def _build_debug_record(result: Any, site_posteriors: Dict[str, float]) -> Dict[str, Any]:
    """
    Convert a BiochemResults entry into a JSON-friendly dict enriched with site posteriors.
    """
    record = asdict(result)
    record['site_posteriors'] = site_posteriors
    record.setdefault('tree_plot', None)
    return record


def analyze_single_gene(
    gene: str,
    alignment_path: Path,
    tree_path: Path,
    caas_positions: List[int],
    output_dir: Path,
    asr_cache_dir: Path,
    taxid_path: Optional[Path] = None,
    caas_metadata_path: Optional[Path] = None,
    trait_file_path: Optional[Path] = None,
    asr_mode: str = "precomputed",
    model: str = "lg",
    threads: int = 1,
    posterior_threshold: float = 0.7,
    run_diagnostics: bool = False,
    focus_positions: bool = False,
    debug_positions: bool = False,
    allow_low_confidence: bool = False,
    include_non_significant: bool = False,
    convergence_mode: str = 'focal_clade'
) -> Dict[str, Any]:
    """
    Run complete CAAS analysis pipeline for a single gene.

    Args:
        gene: Gene name
        alignment_path: Path to alignment file for this gene
        tree_path: Path to tree file
        caas_positions: List of CAAS positions to analyze
        output_dir: Output directory for results
        taxid_path: Optional TaxID mapping file
        caas_metadata_path: Optional CAAS metadata file
        trait_file_path: Optional trait file for contrast definitions
        asr_mode: 'compute' to run ASR, 'precomputed' to load existing
        model: Evolutionary model for ASR
        threads: Number of threads to use for codeml

    Returns:
        Dictionary with analysis results and file paths

    Updated to delegate tree and alignment preprocessing to phylo_single.
    """
    # Initialize generated_files to avoid unbound variable issues
    generated_files = []

    # Add safeguard for caas_positions None
    if caas_positions is None:
        caas_positions = []

    # If no positions specified, try to load all available positions for this gene from metadata
    if not caas_positions and caas_metadata_path and caas_metadata_path.exists():
        try:
            # Load all available positions for this gene from metadata
            # We'll need to import get_caas_info or create a utility function
            caas_positions = load_caas_positions_for_gene(caas_metadata_path, gene)
            logger.info(f"Loaded {len(caas_positions)} CAAS positions from metadata for {gene}")
        except Exception as e:
            logger.warning(f"Could not load CAAS positions from metadata: {e}")

    logger.info(f"Starting analysis for gene: {gene}")
    logger.info(f"Alignment: {alignment_path}")
    logger.info(f"Tree: {tree_path}")
    logger.info(f"CAAS positions: {len(caas_positions)}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"ASR mode: {asr_mode}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_output_dir = output_dir / f"{gene.lower()}_analysis"
    gene_output_dir.mkdir(exist_ok=True)
    
    # Add convergence mode suffix to report directory
    convergence_suffix = "_mrca" if convergence_mode == "mrca" else "_focal"
    report_output_dir = gene_output_dir / f"{gene.lower()}_reports{convergence_suffix}"
    report_output_dir.mkdir(exist_ok=True)

    # Setup per-gene log file with cleanup
    log_file = gene_output_dir / f"{gene.lower()}_analysis.log"
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)

    try:
        # Load alignment and tree data
        logger.info("Loading alignment and tree data...")

        alignment_data = load_alignment_and_mappings(
            alignment_path,
            taxid_path,
            gene
        )

        # Match tree and alignment
        logger.info("Matching tree and alignment...")
        tree_data = load_and_match_tree(
            tree_path,
            alignment_data,
            taxid_path
        )

        # Update references for compatibility
        alignment = alignment_data.alignment
        tree = tree_data.tree

        # Handle ASR step
        asr_results = None
        posterior_data = None

        if asr_mode == 'compute':
            logger.info("Running ASR pipeline...")

            asr_config = SingleGeneASRConfig(
                alignment_path=alignment_path,
                tree_path=tree_path,
                taxid_path=taxid_path,
                output_dir=asr_cache_dir,
                model=model,
                threads=threads,
            posterior_threshold=posterior_threshold,
            run_diagnostics=run_diagnostics
        )

            try:
                validate_asr_inputs(asr_config)
            except (FileNotFoundError, ValueError) as e:
                logger.error(f"ASR validation failed: {e}")
                raise

            asr_results = run_asr_pipeline(
                gene=gene,
                config=asr_config,
                skip_if_exists=True,
                alignment_data=alignment_data,
                tree_data=tree_data
            )

            if asr_results and asr_results.posteriors_node:
                posterior_data = asr_results.posteriors_node

            if asr_results and asr_results.tree_file and asr_results.tree_file.exists():
                original_tip_set = tree_data.tip_set
                try:
                    # Use RST file as source of truth for PAML node IDs
                    ordered_nodes, id_mapping = build_tree_node_mapping(
                        tree_file=asr_results.tree_file,
                        rst_file=asr_results.rst_file
                    )
                    root_node = ordered_nodes[-1]
                    tree_data.nodes = ordered_nodes
                    tree_data.root = root_node
                    tree_data.node_mapping = cast(Dict[int, object], id_mapping)
                    parsed_tips = set(extract_tip_labels(root_node))
                    tree_data.tip_set = parsed_tips if (parsed_tips & original_tip_set) else original_tip_set
                    logger.info("Updated tree_data from ASR (PAML) tree for node/posterior alignment")
                except Exception as e:
                    logger.warning(f"Could not rebuild tree_data from ASR tree: {e}")

        elif asr_mode == 'precomputed':
            logger.info("Loading pre-computed ASR results...")

            asr_config = SingleGeneASRConfig(
                alignment_path=alignment_path,
                tree_path=tree_path,
                taxid_path=taxid_path,
                output_dir=asr_cache_dir,
                model=model,
                threads=threads,
                posterior_threshold=posterior_threshold,
                run_diagnostics=run_diagnostics
            )

            asr_results = load_precomputed_asr(
                gene=gene,
                config=asr_config,
                alignment_data=alignment_data
            )
            posterior_data = asr_results.posteriors_node if asr_results else None

            if asr_results and asr_results.tree_file and asr_results.tree_file.exists():
                original_tip_set = tree_data.tip_set
                try:
                    # Use RST file as source of truth for PAML node IDs
                    ordered_nodes, id_mapping = build_tree_node_mapping(
                        tree_file=asr_results.tree_file,
                        rst_file=asr_results.rst_file
                    )
                    root_node = ordered_nodes[-1]
                    tree_data.nodes = ordered_nodes
                    tree_data.root = root_node
                    tree_data.node_mapping = cast(Dict[int, object], id_mapping)
                    parsed_tips = set(extract_tip_labels(root_node))
                    tree_data.tip_set = parsed_tips if (parsed_tips & original_tip_set) else original_tip_set
                    logger.info("Updated tree_data from precomputed ASR (PAML) tree for node/posterior alignment")
                except Exception as e:
                    logger.warning(f"Could not rebuild tree_data from precomputed ASR tree: {e}")

        else:
            raise ValueError(f"Unknown asr_mode: {asr_mode}")

        # Taxid mapping for biochemical analysis
        try:
            if taxid_path:
                taxid_mapping_param = load_taxid_mapping(taxid_path)
            else:
                taxid_mapping_param = alignment_data.species_to_taxid or {}
            logger.info(f"Using taxid mapping with {len(taxid_mapping_param)} species")
        except Exception as e:
            logger.warning(f"Could not load full taxid mapping: {e}")
            taxid_mapping_param = alignment_data.species_to_taxid or {}

        biochem_results, biochem_diagnostics = analyze_gene_biochemistry(
            gene=gene,
            alignment_data=alignment_data,
            tree_data=tree_data,
            caas_positions=caas_positions,
            caas_metadata_path=caas_metadata_path,
            trait_file_path=trait_file_path,
            taxid_mapping=taxid_mapping_param,
            posterior_data=posterior_data if posterior_data else None,
            posterior_threshold=posterior_threshold,
            allow_low_confidence=allow_low_confidence,
            diagnostics_dir=diag_root if run_diagnostics else None,
            convergence_mode=convergence_mode,
            asr_mode=asr_mode,
            include_non_significant=include_non_significant
        )

        # Use biochem_results directly without scoring
        scored_results = biochem_results

        logger.info(f"Biochemical analysis complete: {len(biochem_results)} positions analyzed")

        debug_records: List[Dict[str, Any]] = []
        debug_position_values: List[int] = []
        selected_results: List[Any] = []
        focus_requested = bool(focus_positions or debug_positions)
        if focus_requested and scored_results:
            selected_results = _select_debug_results(scored_results, count=3)
            debug_position_values = [
                r.position for r in selected_results if getattr(r, 'position', None) is not None
            ]
            if debug_position_values:
                label = "Debug" if debug_positions else "Focus"
                logger.info(f"{label} positions selected: {sorted(debug_position_values)}")
                focus_set = set(debug_position_values)
                for res in scored_results:
                    if getattr(res, 'position', None) in focus_set:
                        try:
                            res.is_focus = True
                        except Exception:
                            pass
        if debug_positions and selected_results:
            for res in selected_results:
                record = _build_debug_record(res, {})
                debug_records.append(record)

        if biochem_diagnostics.get('skipped_positions'):
            logger.info(
                f"Skipped positions: {biochem_diagnostics['skipped_positions']} "
                f"({dict(biochem_diagnostics.get('skip_reasons', {}))})"
            )
        if biochem_diagnostics.get('low_confidence_positions'):
            logger.info(f"Positions with low-confidence nodes: {biochem_diagnostics['low_confidence_positions']}")
        if biochem_diagnostics.get('tip_dump_file'):
            logger.info(f"Tip details file: {biochem_diagnostics['tip_dump_file']}")

        # Generate summary statistics
        stats = summarize_biochemical_results(biochem_results)

        # Filter results based on significance if requested
        if not include_non_significant:
            original_count = len(biochem_results)
            biochem_results = [r for r in biochem_results if getattr(r, 'is_significant', False)]
            scored_results = biochem_results  # Update scored_results to match
            filtered_count = len(biochem_results)
            logger.info(f"Filtered to significant positions only: {filtered_count}/{original_count} positions retained")
        else:
            logger.info("Including all positions (significant and non-significant)")

        # Try to get pass rate and mean score from stats if available
        pass_rate = stats.get('pass_rate', 0) if isinstance(stats, dict) else 0

        # Generate reports and visualizations through the report_single module
        logger.info("Generating reports and visualizations...")

        try:
            reporting_results = create_complete_gene_report(
                gene=gene,
                results=scored_results,
                output_dir=report_output_dir,
                include_plots=True,
                include_heatmaps=True,
                tree_file=asr_results.tree_file if asr_results else None,
                highlight_positions=debug_position_values if (focus_requested and debug_position_values) else None,
                include_non_significant=include_non_significant
            )
        except Exception as report_err:
            logger.warning(f"Reporting failed to execute: {report_err}")
            reporting_results = None

        generated_files = reporting_results.generated_files if reporting_results and reporting_results.generated_files else []
        public_results = reporting_results.public_results if reporting_results else []

        # Compute simple statistics from public results (no quality metrics)
        if public_results:
            significant_count = sum(1 for r in public_results if r.get('is_significant'))
            total_count = len(public_results)
            pass_rate = significant_count / total_count if total_count > 0 else 0.0
            logger.info(f"Analysis summary - Significant positions: {significant_count}/{total_count} ({pass_rate:.1%})")
        else:
            pass_rate = 0.0
            logger.info("No results available for statistics")

        if debug_records:
            tree_plot_map: Dict[int, str] = {}
            if reporting_results:
                for file_info in reporting_results.generated_files or []:
                    if file_info.get('type') == 'gene_tree_plot' and file_info.get('position') is not None:
                        path = file_info.get('path')
                        if path is not None:
                            tree_plot_map[int(file_info['position'])] = path
            for record in debug_records:
                pos = record.get('position')
                if pos is not None and pos in tree_plot_map:
                    record['tree_plot'] = tree_plot_map[pos]
            diag_dir = report_output_dir / "diagnostics"
            diag_dir.mkdir(parents=True, exist_ok=True)
            debug_file = diag_dir / f"{gene.lower()}_debug_positions.json"
            with open(debug_file, 'w', encoding='utf-8') as debug_handle:
                json.dump(debug_records, debug_handle, indent=2, ensure_ascii=False, default=str)
            generated_files.append({'type': 'debug_positions', 'path': str(debug_file)})
            logger.info(f"Debug details written to {debug_file}")

        results_path = None
        for file_info in generated_files:
            if file_info.get('type') == 'json_results' and file_info.get('path'):
                results_path = file_info['path']
                break

        if reporting_results and not reporting_results.success:
            logger.warning(f"Reporting completed with errors: {reporting_results.error_message}")

        # Summary
        success_count = len([r for r in biochem_results if r.pattern_type != 'unknown'])
        logger.info(f"Analysis complete for {gene}:")
        logger.info(f"  - Positions analyzed: {len(biochem_results)}")
        logger.info(f"  - Convergent events detected: {success_count}")
        logger.info(f"  - Report files generated: {len(generated_files)}")

        return {
            'gene': gene,
            'results': public_results,
            'results_file': str(results_path) if results_path else None,
            'statistics': stats,
            'generated_files': generated_files,
            'output_directory': str(gene_output_dir),
            'asr_mode': asr_mode,
            'success': True,
            'positions_analyzed': len(biochem_results),
            'convergent_events': success_count
        }

    except Exception as e:
        logger.error(f"Analysis failed for gene {gene}: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())

        return {
            'gene': gene,
            'success': False,
            'asr_mode': asr_mode,
            'error': str(e),
            'output_directory': str(gene_output_dir),
            'generated_files': generated_files
        }
    finally:
        # Clean up per-gene file handler to prevent accumulation
        root_logger.removeHandler(file_handler)
        file_handler.close()


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description="Single Gene CAAS Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline with ASR
  ./pipeline_main.py --genes BRCA1 --alignment_dir <alignment_directory> --tree data/species_tree.nex --asr-mode compute

  # Use pre-computed ASR results
  ./pipeline_main.py --genes BRCA1 --alignment_dir <alignment_directory> --tree data/species_tree.nex --asr-mode precomputed
        """
    )

    # Required arguments
    parser.add_argument(
        "--genes", "-g",
        nargs='+',
        required=True,
        help="Gene names to analyze (space-separated)"
    )

    parser.add_argument(
        "--alignment_dir", "-a",
        required=True,
        type=Path,
        help="Directory containing alignment files (.phy format) - will auto-discover gene-specific files"
    )

    parser.add_argument(
        "--ensembl-genes-file",
        type=Path,
        help="TSV/CSV with columns (gene, chr, start, end, strand, length); used to strictly whitelist genes"
    )

    parser.add_argument(
        "--tree", "-t",
        required=True,
        type=Path,
        help="Path to phylogenetic tree file (.nex or .nwk format)"
    )

    # Optional arguments
    parser.add_argument(
        "--taxid_mapping",
        type=Path,
        help="Path to TaxID to species name mapping file"
    )

    parser.add_argument(
        "--caas_metadata",
        type=Path,
        help="Path to CAAS metadata file (.output format)"
    )

    parser.add_argument(
        "--trait_file",
        type=Path,
        help="Path to trait file with species contrasts (.tab format)"
    )

    parser.add_argument(
        "--caas_positions",
        nargs='+',
        type=int,
        help="CAAS positions to analyze (overrides metadata if provided)"
    )

    parser.add_argument(
        "--output_dir", "-o",
        type=Path,
        default=Path("./pipeline_output"),
        help="Output directory for results (default: ./pipeline_output)"
    )

    parser.add_argument(
        "--asr-cache-dir",
        type=Path,
        help="Optional directory containing precomputed ASR outputs (rst/tree). Defaults to output_dir."
    )

    parser.add_argument(
        "--asr-mode",
        choices=['compute', 'precomputed'],
        default='precomputed',
        help="ASR handling: 'compute' run ASR, 'precomputed' load existing (default: precomputed)"
    )

    parser.add_argument(
        "--with-asr",
        action="store_true",
        default=False,
        help="Run ASR reconstruction (default: load pre-computed results)"
    )

    parser.add_argument(
        "--asr_model",
        default="lg",
        help="Evolutionary model for ASR (poisson, proportional, dayhoff, jtt, wag, lg; default: lg)"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for codeml/ASR (default: 1)"
    )
    parser.add_argument(
        "--posterior-threshold",
        type=float,
        default=0.7,
        help="Posterior probability threshold for ASR node states (default: 0.7; use 0.0 for debugging)"
    )

    parser.add_argument(
        "--run-diagnostics",
        action='store_true',
        help="Enable diagnostics output: node-level posteriors (TSV) and per-position tip details (JSON) to diagnostics/"
    )

    parser.add_argument(
        "--allow-low-confidence",
        action='store_true',
        help="Continue analysis even when node-level ASR posteriors are low or missing"
    )

    parser.add_argument(
        "--focus-positions",
        action='store_true',
        help="Highlight three sampled positions in plots and tables for visual inspection"
    )

    parser.add_argument(
        "--debug-positions",
        action='store_true',
        help="Enable extra debug output (first three analyzed positions) including annotated trees and TSV dumps"
    )

    parser.add_argument(
        "--include-non-significant",
        action='store_true',
        default=False,
        help="Include non-significant positions in plots and analysis outputs (default: only significant positions)"
    )

    parser.add_argument(
        "--convergence-mode",
        type=str,
        choices=['mrca', 'focal_clade'],
        default='focal_clade',
        help="Source of truth for convergence pattern assessment: 'mrca' uses mrca_contrast node, 'focal_clade' uses all focal_nodes (dynamic multi-pair; default: focal_clade)"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )

    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional path to write pipeline logs (overwrites if exists)"
    )

    parser.add_argument(
        "--config_file",
        type=Path,
        help="Path to JSON configuration file (overrides command line args)"
    )

    args = parser.parse_args()

    # Legacy alias: --with-asr forces compute mode
    if getattr(args, 'with_asr', False):
        args.asr_mode = 'compute'

    configure_logging(verbose=args.verbose, log_file=args.log_file)

    # Validate input directory
    if not args.alignment_dir.exists():
        logger.error(f"Alignment directory not found: {args.alignment_dir}")
        sys.exit(1)

    if not args.tree.exists():
        logger.error(f"Tree file not found: {args.tree}")
        sys.exit(1)

    if args.taxid_mapping and not args.taxid_mapping.exists():
        logger.error(f"TaxID mapping file not found: {args.taxid_mapping}")
        sys.exit(1)

    # Load Ensembl gene whitelist if provided
    ensembl_genes: Optional[Set[str]] = None
    if args.ensembl_genes_file:
        try:
            ensembl_genes = load_ensembl_genes(args.ensembl_genes_file)
            logger.info(f"Loaded {len(ensembl_genes)} genes from {args.ensembl_genes_file}")
        except Exception as exc:
            logger.error(f"Failed to load Ensembl genes file: {exc}")
            sys.exit(1)

    # Verify gene alignment files can be found
    missing_genes = []
    gene_alignments = {}
    for gene in args.genes:
        try:
            alignment_file = find_gene_alignment(args.alignment_dir, gene, ensembl_genes)
            gene_alignments[gene] = alignment_file
            logger.info(f"Found alignment for {gene}: {alignment_file}")
        except FileNotFoundError:
            missing_genes.append(gene)

    if missing_genes:
        logger.error(f"No alignment files found for genes: {', '.join(missing_genes)}")
        logger.error(f"Searched in: {args.alignment_dir}")
        logger.error("Expected patterns: {GENE}*.phy, {GENE}.phy, *{GENE}*.phy")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    asr_cache_dir = args.asr_cache_dir if args.asr_cache_dir else args.output_dir

    # Log start
    logger.info("=" * 80)
    logger.info("SINGLE GENE CAAS ANALYSIS PIPELINE")
    logger.info("=" * 80)
    logger.info(f"Genes to analyze: {', '.join(args.genes)}")
    logger.info(f"Alignment directory: {args.alignment_dir}")
    logger.info(f"Tree: {args.tree}")
    logger.info(f"ASR mode: {args.asr_mode}")
    logger.info(f"Output directory: {args.output_dir}")
    caas_pos_count = len(args.caas_positions) if args.caas_positions else 0
    logger.info(f"CAAS positions: {caas_pos_count} specified")
    logger.info("")

    # Run analysis for each gene
    all_results = []

    for gene in args.genes:
        logger.info("-" * 40)
        logger.info(f"ANALYZING GENE: {gene}")
        logger.info("-" * 40)

        gene_result = analyze_single_gene(
            gene=gene,
            alignment_path=gene_alignments[gene],
            tree_path=args.tree,
            caas_positions=args.caas_positions,
            output_dir=args.output_dir,
            asr_cache_dir=asr_cache_dir,
            taxid_path=args.taxid_mapping,
            caas_metadata_path=args.caas_metadata,
            trait_file_path=args.trait_file,
            asr_mode=args.asr_mode,
            model=args.asr_model,
            threads=args.threads,
            posterior_threshold=args.posterior_threshold,
            run_diagnostics=args.run_diagnostics,
            debug_positions=args.debug_positions,
            focus_positions=args.focus_positions,
            allow_low_confidence=args.allow_low_confidence,
            include_non_significant=args.include_non_significant,
            convergence_mode=args.convergence_mode
        )

        all_results.append(gene_result)

        if gene_result['success']:
            logger.info(f"✓ Analysis successful for {gene}")
        else:
            logger.error(f"✗ Analysis failed for {gene}: {gene_result.get('error', 'Unknown error')}")

        logger.info("")

    # Write summary
    summary_file = args.output_dir / "pipeline_summary.json"
    try:
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False, default=str)
        logger.info(f"✓ Pipeline summary written to: {summary_file}")
    except Exception as e:
        logger.error(f"Failed to write summary: {e}")

    # Final statistics
    successful_runs = sum(1 for r in all_results if r['success'])
    total_positions = sum(r.get('positions_analyzed', 0) for r in all_results if r['success'])
    total_convergent = sum(r.get('convergent_events', 0) for r in all_results if r['success'])

    logger.info("=" * 80)
    logger.info("PIPELINE EXECUTION SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total genes processed: {len(args.genes)}")
    logger.info(f"Successful runs: {successful_runs}/{len(args.genes)}")
    logger.info(f"Total CAAS positions analyzed: {total_positions}")
    logger.info(f"Convergent evolution events detected: {total_convergent}")
    logger.info(f"Results saved to: {args.output_dir}")
    logger.info("=" * 80)

    # Exit with appropriate code
    if successful_runs == len(args.genes):
        logger.info("🎉 All gene analyses completed successfully!")
        sys.exit(0)
    else:
        logger.warning(f"⚠️  {len(args.genes) - successful_runs} gene analyses failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
