#!/usr/bin/env python3
# filter_caas_clusters-param.py — Flag the CAAS positions that lie in dense clusters ("trains").
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Cluster filter: flags the positions of a gene that lie in a high-density run of CAAS.

The positions of a tight cluster ("train") of substitutions are flagged "Discarded", so that
clustered substitutions do not pass as independent CAAS.
For each gene (and each caap_group when that column exists) every interval [start, end]
of the sorted positions with span >= minlen and count / span >= maxcaas flags all the
positions it holds (core.postproc.ctrain, shared with the permulation null). With --map-dir
the span is measured in untrimmed alignment columns.

Called by:  CT_FILTER process (ctpp_clustfilter.nf), once per (minlen, maxcaas) pair
Inputs:     -i  TSV with at least Gene and Position (caap_group is used when present)
            -l, -c  minimum interval span and maximum density (maxcaas)
            --map-dir  optional directory of per-gene MAP tables (untrimmed coordinates)
Outputs:    <input stem>.filtered.minlen<L>.maxcaas<C*100>.tsv with Gene, Position, [caap_group],
            clustering_flag ("Good" or "Discarded"); <input stem>.minlen<L>.maxcaas<C*100>.log

Usage:
  python filter_caas_clusters-param.py -i input.tsv -c 0.7 -l 3 [-v]
"""

import pandas as pd
import argparse
import logging
import sys
from pathlib import Path

# core.postproc holds the one implementation of trains, shared by the observed chain and the null.
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION" / "local"))
from src.core.columns import gene_columns, index_files  # noqa: E402
from src.core.postproc import ctrain  # noqa: E402

# ── Argument Parsing ──────────────────────────────────────────────────────────

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
        "--map-dir",
        type=str,
        default=None,
        help="Directory of the trimmer's per-gene MAP tables. Trains then measure their span in untrimmed "
             "alignment columns; a gene without a MAP keeps trimmed coordinates. The permulation null must "
             "use the same directory (--train-map-dir)."
    )
    parser.add_argument(
        "--map-suffix",
        type=str,
        default=".map.tsv",
        help="File-name tail of the MAP tables"
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

# ── Logger Setup ──────────────────────────────────────────────────────────────

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

# ── Main Filtering Function ───────────────────────────────────────────────────

def filterCAAS(infile, maxcaas, minlen, logger, map_dir=None, map_suffix=".map.tsv"):
    """
    Filter CAAS positions by identifying and marking high-density clusters.
    
    Reads a tab-separated CAAS file, processes each gene independently to identify
    clustered positions, and outputs a filtered file with clustering flags.
    
    Args:
        infile: Path to input .caas file (tab-separated with Gene, Position columns)
        maxcaas: Maximum density threshold (0.0 to 1.0)
        minlen: Minimum interval length for clustering detection
        logger: Configured logger instance
        map_dir: optional directory of MAP tables; with it, spans are measured in untrimmed columns
        map_suffix: file-name tail of the MAP tables
        
    Returns:
        Path object to the output filtered file
        
    Raises:
        SystemExit: If input file is invalid or missing required columns
    """
    path = Path(infile)
    if not path.is_file():
        logger.error(f"Input file not found: {infile}")
        sys.exit(1)
    
    # Validate the file format and the required columns
    try:
        # keep_default_na=False + na_values=[""]: the CAAS table has categorical amino-acid
        # columns (caas, amino_encoded, derived_residues) whose values can be NA-sentinel
        # strings ("N/A" is Asn on the changed side against Ala), and the default NA parsing
        # of pandas would blank them. Only an empty cell is missing here.
        df = pd.read_csv(path, sep="\t", header=0,
                         keep_default_na=False, na_values=["", "nan", "NaN"])
        
        for col in ("Gene", "Position"):
            if col not in df.columns:
                logger.error(f"Missing required column: {col}")
                sys.exit(1)
        
        # Position must be integer
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
    map_index = index_files(map_dir, map_suffix) if map_dir else None
    genes_without_map = []
    if map_index is not None:
        logger.info(f"Untrimmed-coordinate trains: {len(map_index)} MAP files in {map_dir}")
    genes = df["Gene"].unique()
    total_genes = len(genes)
    
    # With a caap_group column every group of a gene is processed on its own
    has_caap_group = "caap_group" in df.columns
    
    if has_caap_group:
        caap_groups = df["caap_group"].unique()
        logger.info(
            f"CAAP mode detected: Processing {total_genes} genes across "
            f"{len(caap_groups)} CAAP groups ({', '.join(caap_groups)}) with "
            f"maxcaas={maxcaas}, minlen={minlen}"
        )
    else:
        logger.info(
            f"Processing {total_genes} genes with "
            f"maxcaas={maxcaas}, minlen={minlen}"
        )
    
    position_col = df["Position"]
    gene_col = df["Gene"]
    caap_group_col = df["caap_group"] if has_caap_group else None
    
    # Each gene is independent
    for i, gene in enumerate(genes, 1):
        gene_df = df[gene_col == gene]
        columns = gene_columns(map_index, gene, map_suffix) if map_index is not None else None
        if map_index is not None and columns is None:
            genes_without_map.append(gene)
        
        if has_caap_group:
            groups_in_gene = gene_df["caap_group"].unique()
            logger.info(f"Gene [{i}/{total_genes}]: {gene} (Groups: {', '.join(groups_in_gene)})")
            
            for group in groups_in_gene:
                group_df = gene_df[gene_df["caap_group"] == group]
                positions = sorted(group_df["Position"].unique())
                
                if len(positions) < minlen:
                    continue
                
                logger.debug(
                    f"  Group {group}: {len(positions)} positions: "
                    f"{positions[:5]}{'...' if len(positions) > 5 else ''}"
                )
                
                group_discarded = ctrain(positions, maxcaas, minlen, columns)
                
                if group_discarded:
                    logger.info(
                        f"  Group {group}: {len(group_discarded)} positions flagged "
                        f"(density threshold exceeded)"
                    )
                    for pos in group_discarded:
                        discarded.append((gene, pos, group))
        else:
            # No caap_group column: all positions of the gene form one unit
            logger.info(f"Gene [{i}/{total_genes}]: {gene}")
            positions = sorted(gene_df["Position"].unique())
            
            if len(positions) < minlen:
                continue
            
            logger.debug(
                f"  {len(positions)} positions: "
                f"{positions[:5]}{'...' if len(positions) > 5 else ''}"
            )
            
            gene_discarded = ctrain(positions, maxcaas, minlen, columns)
            
            if gene_discarded:
                logger.info(
                    f"  {len(gene_discarded)} positions flagged "
                    f"(density threshold exceeded)"
                )
                for pos in gene_discarded:
                    discarded.append((gene, pos, None))
    
    if genes_without_map:
        logger.warning(
            f"{len(genes_without_map)} of {total_genes} genes have no MAP file and keep trimmed coordinates, "
            f"e.g. {sorted(genes_without_map)[:5]}"
        )

    out = df.copy()
    
    # A row is Discarded when its (Gene, Position[, caap_group]) was flagged
    if has_caap_group:
        out["clustering_flag"] = out.apply(
            lambda row: "Discarded" if (row["Gene"], row["Position"], row["caap_group"]) in discarded else "Good",
            axis=1
        )
    else:
        discarded_set = {(gene, pos) for gene, pos, _ in discarded}
        out["clustering_flag"] = out.apply(
            lambda row: "Discarded" if (row["Gene"], row["Position"]) in discarded_set else "Good",
            axis=1
        )
    
    # Only the keys and the flag are written (caap_group included when present)
    output_cols = ["Gene", "Position", "clustering_flag"]
    if "caap_group" in df.columns:
        output_cols = ["Gene", "Position", "caap_group", "clustering_flag"]
        out = out[output_cols]
    else:
        out = out[output_cols]
    
    out_file = path.with_suffix(
        f".filtered.minlen{minlen}.maxcaas{int(maxcaas*100)}.tsv"
    )
    out.to_csv(out_file, sep="\t", index=False)
    
    logger.info(f"Output written to: {out_file}")
    logger.info(f"Total genes processed: {total_genes}")
    logger.info(f"Total positions discarded: {len(discarded)}")
    logger.info(
        f"Positions retained: {len(out[out['clustering_flag'] == 'Good'])}"
    )
    
    return out_file

# ── Entry Point ───────────────────────────────────────────────────────────────

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
            logger,
            args.map_dir,
            args.map_suffix
        )
        logger.info("✓ Clustering analysis completed successfully")
        print(f"\nFiltered output: {output_path}")
        
    except Exception as e:
        logger.error(f"✗ Processing failed: {e}")
        sys.exit(1)