#!/usr/bin/env python3

import pandas as pd
import argparse
import logging
import sys
from pathlib import Path
import numpy as np

# ---------------------- Argument Parsing ----------------------
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

# ---------------------- Logger Setup ----------------------
def setup_logger(input_path: Path, maxcaas: float, minlen: int, verbose: bool):
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

# ---------------------- Pruning Logic ----------------------
def ctrain(position_list, maxcaas=0.5, minlen=10, logger=None):
    unique = np.sort(np.unique(position_list))  # Use numpy for faster sorting and unique
    n = len(unique)
    bad_positions = set()
    intervals = []  # merged bad intervals
    if n < minlen:
        return []
    for r in range(n):
        for l in range(r+1):
            start, end = unique[l], unique[r]
            span = end - start + 1
            if span < minlen:
                continue
            count = r - l + 1
            density = count / span
            if logger:
                logger.debug(f"Pair: {start}-{end}, Count: {count}, Span: {span}, Density: {density:.3f}")
            if density >= maxcaas:
                bad_positions.update(unique[l:r+1].tolist())  # Convert to list for update
                # Merge new interval
                ns, ne = start, end
                merged = []
                for ds, de in intervals:
                    if de < ns or ds > ne:
                        merged.append((ds, de))
                    else:
                        ns, ne = min(ns, ds), max(ne, de)
                merged.append((ns, ne))
                intervals = merged
    return sorted(bad_positions)

# ---------------------- Main Filter ----------------------
def filterCAAS(infile, maxcaas, minlen, logger):
    path = Path(infile)
    if not path.is_file():
        logger.error(f"Input '{infile}' not found.")
        sys.exit(1)
    
    # Validate input file format
    try:
        df = pd.read_csv(path, sep="\t", header=0)
        for col in ("Gene", "Position"):
            if col not in df.columns:
                logger.error(f"Missing column: {col}")
                sys.exit(1)
        # Validate Position column contains integers
        if not pd.api.types.is_integer_dtype(df["Position"]):
            try:
                df["Position"] = pd.to_numeric(df["Position"], downcast="integer", errors="raise")
            except ValueError:
                logger.error("Position column contains non-integer values.")
                sys.exit(1)
    except pd.errors.ParserError:
        logger.error(f"Input file '{infile}' is not a valid tab-separated file.")
        sys.exit(1)

    discarded = []
    genes = df["Gene"].unique()
    total_genes = len(genes)
    logger.info(f"Processing {total_genes} genes with maxcaas={maxcaas}, minlen={minlen}")
    
    # Cache Position column access
    position_col = df["Position"]
    gene_col = df["Gene"]
    
    for i, gene in enumerate(genes, 1):
        logger.info(f"[{i}/{total_genes}] {gene}")
        poses = position_col[gene_col == gene].tolist()
        bad = ctrain(poses, maxcaas, minlen, logger)
        if bad:  # Log discarded positions if any
            logger.debug(f"Discarded positions for gene {gene}: {', '.join(map(str, bad))}")
        else:
            logger.debug(f"No positions discarded for gene {gene}")
        discarded += [(gene, p) for p in bad]
    
    # Label
    if discarded:
        disc_df = pd.DataFrame(discarded, columns=["Gene", "Position"])
        disc_df["ClusteringFlag"] = "Discarded"
    else:
        disc_df = pd.DataFrame(columns=["Gene", "Position", "ClusteringFlag"])
    
    out = (df.merge(disc_df, on=["Gene", "Position"], how="left")
        .assign(ClusteringFlag=lambda d: d["ClusteringFlag"].fillna("Good")))
    # Select only Gene, Position, and ClusteringFlag for output
    out = out[["Gene", "Position", "ClusteringFlag"]]
    out_file = path.with_suffix(f".filtered.minlen{minlen}.maxcaas{int(maxcaas*100)}.tsv")
    out.to_csv(out_file, sep="\t", index=False)
    
    logger.info(f"Filtered file written to: {out_file}")
    logger.info(f"Total discarded positions: {len(discarded)}")
    logger.info(f"Completed clustering of {total_genes} genes.")
    
    return out_file

# ---------------------- Entry Point ----------------------
if __name__ == "__main__":
    args = parse_args()
    logger = setup_logger(Path(args.inputfile), args.maxcaas, args.minlen, args.verbose)
    try:
        output_path = filterCAAS(args.inputfile, args.maxcaas, args.minlen, logger)
        logger.info("Completed successfully.")
        print(f"Filtered file written to: {output_path}")
    except Exception as e:
        logger.error(f"Failed: {e}")
        sys.exit(1)