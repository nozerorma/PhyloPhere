#!/usr/bin/env python3
# filter_caas_clusters-param.py — Flag the CAAS positions that lie in dense clusters ("trains").
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Cluster filter: flags the positions of a gene that lie in a high-density run of CAAS.

The positions of a tight cluster ("train") of substitutions are flagged "Discarded", so that
clustered substitutions do not pass as independent CAAS.
For each gene every interval [start, end]
of the sorted positions with span >= minlen and count / span >= maxcaas flags all the
positions it holds (core.postproc.ctrain, shared with the permulation null). With --map-dir
the span is measured in untrimmed alignment columns.

Called by:  CT_FILTER process (ctpp_clustfilter.nf), once per (minlen, maxcaas) pair
Inputs:     -i  TSV with at least Gene and Position
            -l, -c  minimum interval span and maximum density (maxcaas)
            --map-dir  optional directory of per-gene MAP tables (untrimmed coordinates)
Outputs:    <input stem>.filtered.minlen<L>.maxcaas<C*100>.tsv with Gene, Position,
            clustering_flag ("Good" or "Discarded"); <input stem>.minlen<L>.maxcaas<C*100>.log

Usage:
  python filter_caas_clusters-param.py -i input.tsv -c 0.7 -l 3 [-v]
"""

from postproc_io import read_table, positions as validated_positions, parameter_token
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
    parser.add_argument("--min-positions", type=int, default=None, help="Minimum candidate count; default follows --minlen")
    args = parser.parse_args()
    if not (0.0 <= args.maxcaas <= 1.0):
        parser.error("--maxcaas must be between 0 and 1.")
    if args.min_positions is not None and args.min_positions < 1:
        parser.error("--min-positions must be positive")
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
    log_file = input_path.with_name(f"{name}.minlen{minlen}.maxcaas{parameter_token(maxcaas)}.log")
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

def filterCAAS(infile, maxcaas, minlen, logger, map_dir=None, map_suffix=".map.tsv", min_positions=None):
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
    
    df = read_table(path)
    if not {'Gene', 'Position'} <= set(df):
        raise ValueError('Missing Gene/Position')
    df['Position'] = validated_positions(df.Position)

    discarded = []
    map_index = index_files(map_dir, map_suffix, df.Gene.unique()) if map_dir else None
    genes_without_map = []
    if map_index is not None:
        logger.info(f"Untrimmed-coordinate trains: {len(map_index)} MAP files in {map_dir}")
    genes = df["Gene"].unique()
    total_genes = len(genes)
    
    logger.info(
        f"Processing {total_genes} genes with "
        f"maxcaas={maxcaas}, minlen={minlen}"
    )

    # All candidates, including side=none; score-independent quality mask.
    for i, (gene, gene_df) in enumerate(df.groupby("Gene", sort=False), 1):
        columns = gene_columns(map_index, gene, map_suffix) if map_index is not None else None
        if map_index is not None and columns is None:
            genes_without_map.append(gene)

        logger.info(f"Gene [{i}/{total_genes}]: {gene}")
        positions = sorted(gene_df["Position"].unique())

        if len(positions) < (minlen if min_positions is None else min_positions):
            continue

        logger.debug(
            f"  {len(positions)} positions: "
            f"{positions[:5]}{'...' if len(positions) > 5 else ''}"
        )

        gene_discarded = ctrain(positions, maxcaas, minlen, columns, min_positions)

        if gene_discarded:
            logger.info(
                f"  {len(gene_discarded)} positions flagged "
                f"(density threshold exceeded)"
            )
            for pos in gene_discarded:
                discarded.append((gene, pos))

    if genes_without_map:
        logger.warning(
            f"{len(genes_without_map)} of {total_genes} genes have no MAP file and keep trimmed coordinates, "
            f"e.g. {sorted(genes_without_map)[:5]}"
        )

    # Only the keys and the flag are written; a row is Discarded when its (Gene, Position) was flagged
    out = df[["Gene", "Position"]].drop_duplicates().copy()
    is_discarded = pd.MultiIndex.from_frame(out).isin(set(discarded))
    out["clustering_flag"] = ["Discarded" if d else "Good" for d in is_discarded]

    out_file = path.with_suffix(
        f".filtered.minlen{minlen}.maxcaas{parameter_token(maxcaas)}.tsv"
    )
    out.to_csv(out_file, sep="\t", index=False)
    import json
    from src.core.scoring_provenance import digest
    Path(str(out_file)+'.scoring.json').write_text(json.dumps(dict(
        enabled=True,minlen=minlen,maxcaas=maxcaas,min_positions=min_positions,
        maps={str(g):digest(map_index.get(g)) if map_index else None for g in genes})))
    
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
            args.map_suffix,
            args.min_positions
        )
        logger.info("✓ Clustering analysis completed successfully")
        print(f"\nFiltered output: {output_path}")
        
    except Exception as e:
        logger.error(f"✗ Processing failed: {e}")
        sys.exit(1)