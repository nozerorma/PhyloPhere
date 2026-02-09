#!/usr/bin/env python3
"""
CAAS Cluster Filtering Script

This script identifies and filters high-density clusters of CAAS (Context-Aware Amino Acid 
Substitutions) positions within genes. Positions are flagged as "Discarded" if they fall 
within regions where the density of CAAS exceeds a threshold over a minimum interval length.

Algorithm:
-----------
For each gene, the script evaluates all possible position intervals [start, end] where:
  - Interval span ≥ minlen (minimum length threshold)
  - Density = count_positions / span ≥ maxcaas (maximum density threshold)

Positions within any interval exceeding the density threshold are marked for removal,
preventing false CAAS discovery from tightly clustered substitutions ("trains").

Usage:
------
  python filter_caas_clusters-param.py -i input.caas -c 0.7 -l 3 [-v]

Input Format:
-------------
  Tab-separated file with columns: Gene, Position (integer), [other columns...]

Output:
-------
  - Filtered TSV file: input.filtered.minlenX.maxcaasY.tsv
    Columns: Gene, Position, ClusteringFlag ("Good" or "Discarded")
  - Log file: input.minlenX.maxcaasY.log
"""

import pandas as pd
import argparse
import logging
import sys
from pathlib import Path
import numpy as np

# ============================================================================
# Argument Parsing
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="CAAS Train Hack with dynamic pruning and logging")
    parser.add_argument(
        "--inputfile", "-i",
        type=str,
        required=True,
        help="Path to input CAAS discovery (.caas) file"
    )
    parser.add_argument(
        "--maxcaas", "-c",
        type=float,
        default=0.7,
        help="Maximum CAAS density threshold (0 <= maxcaas <= 1)"
    )
    parser.add_argument(
        "--minlen", "-l",
        type=int,
        default=3,
        help="Minimum interval length (>=1)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose debug output"
    )
    args = parser.parse_args()
    if not (0.0 <= args.maxcaas <= 1.0):
        parser.error("--maxcaas must be between 0 and 1.")
    if args.minlen < 1:
        parser.error("--minlen must be at least 1.")
    return args

# ============================================================================
# Logger Setup
# ============================================================================

def setup_logger(input_path: Path, maxcaas: float, minlen: int, verbose: bool):
    """
    Configure logging to both file and console.
    
    Args:
        input_path: Path to input file (used to name log file)
        maxcaas: Maximum density threshold (for log filename)
        minlen: Minimum interval length (for log filename)
        verbose: If True, set DEBUG level; otherwise INFO
        
    Returns:
        Configured logger instance
    """
    name = input_path.stem
    log_file = input_path.with_name(f"{name}.minlen{minlen}.maxcaas{int(maxcaas*100)}.log")
    level = logging.DEBUG if verbose else logging.INFO
    logger = logging.getLogger("CAAS_HACK")
    logger.setLevel(level)
    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logger.addHandler(ch)
    return logger

# ============================================================================
# CAAS Clustering Detection Algorithm
# ============================================================================

def ctrain(position_list, maxcaas=0.5, minlen=10, logger=None):
    """
    Identify CAAS positions within high-density clusters ("trains").
    
    This function evaluates all possible intervals of CAAS positions within a gene,
    identifying those where the density (positions per span) exceeds the threshold.
    Positions in any high-density interval are marked for removal.
    
    Algorithm:
    ----------
    1. Sort and extract unique positions
    2. For each pair of positions (l, r):
       - Calculate interval span: end - start + 1
       - Calculate density: position_count / span
       - If density ≥ maxcaas and span ≥ minlen, mark all positions in [l, r] as bad
    3. Return sorted list of bad positions
    
    Args:
        position_list: List of integer positions for a single gene
        maxcaas: Maximum allowed density threshold (0.0 to 1.0)
        minlen: Minimum interval length to consider
        logger: Optional logger for debug output
        
    Returns:
        Sorted list of positions to discard (empty list if none found)
        
    Example:
        Positions: [10, 11, 12, 50]
        With maxcaas=0.7, minlen=3:
        - Interval [10, 12]: span=3, count=3, density=1.0 > 0.7 → discard [10, 11, 12]
        - Interval [10, 50]: span=41, count=4, density=0.098 < 0.7 → keep
    """
    unique = np.sort(np.unique(position_list))
    n = len(unique)
    
    # Early exit if not enough positions to form a cluster
    if n < minlen:
        return []
    
    bad_positions = set()
    
    # Evaluate all position pairs (l, r) where l ≤ r
    for r in range(n):
        for l in range(r + 1):
            start, end = unique[l], unique[r]
            span = end - start + 1
            
            # Skip intervals shorter than minimum length
            if span < minlen:
                continue
            
            count = r - l + 1
            density = count / span
            
            if logger:
                logger.debug(
                    f"Interval [{start}, {end}]: "
                    f"count={count}, span={span}, density={density:.3f}"
                )
            
            # Flag positions in high-density intervals
            if density >= maxcaas:
                bad_positions.update(unique[l:r+1].tolist())
    
    return sorted(bad_positions)

# ============================================================================
# Main Filtering Function
# ============================================================================

def filterCAAS(infile, maxcaas, minlen, logger):
    """
    Filter CAAS positions by identifying and marking high-density clusters.
    
    Reads a tab-separated CAAS file, processes each gene independently to identify
    clustered positions, and outputs a filtered file with clustering flags.
    
    Args:
        infile: Path to input .caas file (tab-separated with Gene, Position columns)
        maxcaas: Maximum density threshold (0.0 to 1.0)
        minlen: Minimum interval length for clustering detection
        logger: Configured logger instance
        
    Returns:
        Path object to the output filtered file
        
    Raises:
        SystemExit: If input file is invalid or missing required columns
    """
    path = Path(infile)
    if not path.is_file():
        logger.error(f"Input file not found: {infile}")
        sys.exit(1)
    
    # Validate input file format and required columns
    try:
        df = pd.read_csv(path, sep="\t", header=0)
        
        # Check for required columns
        for col in ("Gene", "Position"):
            if col not in df.columns:
                logger.error(f"Missing required column: {col}")
                sys.exit(1)
        
        # Ensure Position column contains integers
        if not pd.api.types.is_integer_dtype(df["Position"]):
            try:
                df["Position"] = pd.to_numeric(
                    df["Position"], 
                    downcast="integer", 
                    errors="raise"
                )
            except ValueError:
                logger.error("Position column contains non-integer values")
                sys.exit(1)
                
    except pd.errors.ParserError:
        logger.error(f"Input file is not a valid tab-separated file: {infile}")
        sys.exit(1)

    discarded = []
    genes = df["Gene"].unique()
    total_genes = len(genes)
    
    logger.info(
        f"Processing {total_genes} genes with "
        f"maxcaas={maxcaas}, minlen={minlen}"
    )
    
    # Cache column access for performance
    position_col = df["Position"]
    gene_col = df["Gene"]
    
    # Process each gene independently
    for i, gene in enumerate(genes, 1):
        logger.info(f"Gene [{i}/{total_genes}]: {gene}")
        
        # Extract positions for current gene
        poses = position_col[gene_col == gene].tolist()
        
        # Identify clustered positions
        bad = ctrain(poses, maxcaas, minlen, logger)
        
        if bad:
            logger.debug(
                f"  Discarded {len(bad)} positions: "
                f"{', '.join(map(str, bad))}"
            )
        else:
            logger.debug("  No clustered positions found")
        
        # Track discarded positions
        discarded += [(gene, p) for p in bad]
    
    # Create clustering flag annotations
    if discarded:
        disc_df = pd.DataFrame(discarded, columns=["Gene", "Position"])
        disc_df["ClusteringFlag"] = "Discarded"
    else:
        disc_df = pd.DataFrame(columns=["Gene", "Position", "ClusteringFlag"])
    
    # Merge flags with original data and fill missing as "Good"
    out = (
        df.merge(disc_df, on=["Gene", "Position"], how="left")
        .assign(ClusteringFlag=lambda d: d["ClusteringFlag"].fillna("Good"))
    )
    
    # Output only essential columns
    out = out[["Gene", "Position", "ClusteringFlag"]]
    
    # Generate output filename
    out_file = path.with_suffix(
        f".filtered.minlen{minlen}.maxcaas{int(maxcaas*100)}.tsv"
    )
    out.to_csv(out_file, sep="\t", index=False)
    
    # Summary statistics
    logger.info(f"Output written to: {out_file}")
    logger.info(f"Total genes processed: {total_genes}")
    logger.info(f"Total positions discarded: {len(discarded)}")
    logger.info(
        f"Positions retained: {len(out[out['ClusteringFlag'] == 'Good'])}"
    )
    
    return out_file

# ============================================================================
# Entry Point
# ============================================================================

if __name__ == "__main__":
    args = parse_args()
    logger = setup_logger(
        Path(args.inputfile), 
        args.maxcaas, 
        args.minlen, 
        args.verbose
    )
    
    try:
        output_path = filterCAAS(
            args.inputfile, 
            args.maxcaas, 
            args.minlen, 
            logger
        )
        logger.info("✓ Clustering analysis completed successfully")
        print(f"\nFiltered output: {output_path}")
        
    except Exception as e:
        logger.error(f"✗ Processing failed: {e}")
        sys.exit(1)