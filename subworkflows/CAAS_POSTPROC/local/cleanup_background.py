#!/usr/bin/env python3
"""
Background Gene List Cleanup Module

Removes problematic genes (extreme/dubious) from CT background files to ensure
consistency between discovery and background datasets for downstream enrichment analyses.

Background files contain all genes tested in CAAS discovery, regardless of significance.
This cleanup ensures that genes removed from discovery (due to extreme density or clustering)
are also removed from background, preventing bias in enrichment statistics.

Author: PhyloPhere Pipeline
License: GPL-3.0
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def load_removed_genes(summary_file, trait_col='Trait'):
    """
    Load list of genes to remove from background.
    
    Args:
        summary_file: Path to removed_genes_summary.tsv
        trait_col: Name of trait column
    
    Returns:
        DataFrame with columns [Trait, Gene, Category]
    """
    if not Path(summary_file).exists():
        print(f"Error: Removed genes summary file not found: {summary_file}", file=sys.stderr)
        sys.exit(1)
    
    removed_df = pd.read_csv(summary_file, sep='\t')
    
    # Validate required columns
    required_cols = [trait_col, 'Gene']
    missing_cols = [col for col in required_cols if col not in removed_df.columns]
    if missing_cols:
        print(f"Error: Removed genes summary missing columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(removed_df)} removed gene entries across {removed_df[trait_col].nunique()} traits", file=sys.stderr)
    
    return removed_df


def load_background_file(bg_file, trait_name, trait_col='Trait'):
    """
    Load background gene list from CT discovery output.
    
    Args:
        bg_file: Path to background file (TSV or single-column TXT)
        trait_name: Trait name to annotate
        trait_col: Name of trait column
    
    Returns:
        DataFrame with columns [Trait, Gene]
    """
    if not Path(bg_file).exists():
        print(f"Warning: Background file not found: {bg_file}", file=sys.stderr)
        return pd.DataFrame(columns=[trait_col, 'Gene'])
    
    # Try reading as TSV first, then as single-column text
    try:
        # Attempt TSV with header
        bg_df = pd.read_csv(bg_file, sep='\t')
        if 'Gene' not in bg_df.columns:
            # No 'Gene' column, assume first column is gene names
            bg_df = pd.read_csv(bg_file, sep='\t', header=None, names=['Gene'])
    except Exception:
        # Fallback: single-column text file
        bg_df = pd.read_csv(bg_file, header=None, names=['Gene'])
    
    bg_df[trait_col] = trait_name
    bg_df = bg_df[[trait_col, 'Gene']]  # Reorder columns
    
    print(f"Loaded {len(bg_df)} background genes for trait: {trait_name}", file=sys.stderr)
    
    return bg_df


def clean_background(bg_df, removed_genes_df, trait_col='Trait'):
    """
    Remove flagged genes from background dataset.
    
    Uses anti-join to exclude genes present in removed_genes_df.
    
    Args:
        bg_df: Background gene DataFrame [Trait, Gene]
        removed_genes_df: Removed genes DataFrame [Trait, Gene, ...]
        trait_col: Name of trait column
    
    Returns:
        Cleaned background DataFrame
    """
    # Get unique genes to remove for this trait
    trait_name = bg_df[trait_col].iloc[0]
    genes_to_remove = removed_genes_df[removed_genes_df[trait_col] == trait_name]['Gene'].unique()
    
    # Anti-join: keep only genes NOT in removed list
    cleaned_df = bg_df[~bg_df['Gene'].isin(genes_to_remove)].copy()
    
    n_removed = len(bg_df) - len(cleaned_df)
    print(f"  {trait_name}: Removed {n_removed} genes ({len(genes_to_remove)} unique flagged genes)", file=sys.stderr)
    print(f"  Retained: {len(cleaned_df)}/{len(bg_df)} genes ({100*len(cleaned_df)/len(bg_df):.1f}%)", file=sys.stderr)
    
    return cleaned_df


def write_cleaned_background(cleaned_df, output_file, trait_name):
    """
    Write cleaned background gene list.
    
    Args:
        cleaned_df: Cleaned background DataFrame
        output_file: Output file path
        trait_name: Trait name for logging
    """
    # Write gene list only (single column, no header)
    cleaned_df['Gene'].to_csv(output_file, index=False, header=False, sep='\n')
    print(f"Wrote cleaned background for {trait_name}: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Clean background gene lists by removing problematic genes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Clean single background file
  python cleanup_background.py -s removed_genes_summary.tsv \
    -b malignant.background.txt -t malignant_prevalence \
    -o cleaned_malignant.txt
  
  # Clean multiple backgrounds via directory
  python cleanup_background.py -s removed_genes_summary.tsv \
    -d backgrounds/ -o cleaned_backgrounds/
"""
    )
    
    # Input files
    parser.add_argument('-s', '--summary', required=True,
                        help='Removed genes summary file from gene filtering')
    parser.add_argument('-b', '--background-file', default=None,
                        help='Single background file to clean')
    parser.add_argument('-d', '--background-dir', default=None,
                        help='Directory containing background files (alternative to -b)')
    parser.add_argument('-p', '--pattern', default='*.background.txt',
                        help='File pattern for background files (default: *.background.txt)')
    
    # Trait mapping
    parser.add_argument('-t', '--trait', default=None,
                        help='Trait name for single background file (required with -b)')
    parser.add_argument('--trait-col', default='Trait',
                        help='Name of trait column (default: Trait)')
    parser.add_argument('--extract-trait-from-filename', action='store_true',
                        help='Extract trait name from filename (for -d mode)')
    
    # Output
    parser.add_argument('-o', '--output', required=True,
                        help='Output file (for -b mode) or directory (for -d mode)')
    
    args = parser.parse_args()
    
    # Validate input mode
    if args.background_file and args.background_dir:
        print("Error: Specify either --background-file or --background-dir, not both", file=sys.stderr)
        sys.exit(1)
    
    if not args.background_file and not args.background_dir:
        print("Error: Must specify either --background-file or --background-dir", file=sys.stderr)
        sys.exit(1)
    
    if args.background_file and not args.trait:
        print("Error: --trait required when using --background-file", file=sys.stderr)
        sys.exit(1)
    
    # Load removed genes
    print(f"Loading removed genes summary: {args.summary}", file=sys.stderr)
    removed_genes_df = load_removed_genes(args.summary, trait_col=args.trait_col)
    
    # Process backgrounds
    if args.background_file:
        # Single file mode
        print(f"\nProcessing single background file: {args.background_file}", file=sys.stderr)
        bg_df = load_background_file(args.background_file, args.trait, trait_col=args.trait_col)
        
        if len(bg_df) == 0:
            print(f"Warning: No genes loaded from {args.background_file}", file=sys.stderr)
            return
        
        cleaned_df = clean_background(bg_df, removed_genes_df, trait_col=args.trait_col)
        write_cleaned_background(cleaned_df, args.output, args.trait)
        
    else:
        # Directory mode
        print(f"\nProcessing background directory: {args.background_dir}", file=sys.stderr)
        bg_dir = Path(args.background_dir)
        
        if not bg_dir.exists():
            print(f"Error: Background directory not found: {bg_dir}", file=sys.stderr)
            sys.exit(1)
        
        # Find background files
        bg_files = list(bg_dir.glob(args.pattern))
        
        if len(bg_files) == 0:
            print(f"Warning: No background files found matching pattern: {args.pattern}", file=sys.stderr)
            # Create a dummy output file to satisfy Nextflow expectations
            output_dir = Path(args.output)
            output_dir.mkdir(parents=True, exist_ok=True)
            dummy_file = output_dir / "cleaned_background_none.txt"
            dummy_file.write_text("# No background files found to process\n")
            print(f"Created placeholder file: {dummy_file}", file=sys.stderr)
            return
        
        print(f"Found {len(bg_files)} background files", file=sys.stderr)
        
        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each background file
        for bg_file in bg_files:
            # Extract trait name
            if args.extract_trait_from_filename:
                # Handle different naming formats:
                # - background_genes.output -> "genes"
                # - {trait}.background.txt -> {trait}
                # - background.output -> "background"
                stem = bg_file.stem  # filename without extension
                if 'background_genes' in stem:
                    trait_name = "genes"
                elif '.background' in stem:
                    trait_name = stem.replace('.background', '')
                else:
                    trait_name = stem
            else:
                # Try to match with trait names in removed_genes_df
                trait_names = removed_genes_df[args.trait_col].unique()
                trait_name = None
                for tn in trait_names:
                    if tn in str(bg_file):
                        trait_name = tn
                        break
                
                if trait_name is None:
                    print(f"Warning: Could not determine trait for {bg_file}, skipping", file=sys.stderr)
                    continue
            
            print(f"\nProcessing {bg_file.name} (trait: {trait_name})", file=sys.stderr)
            
            # Load, clean, and write
            bg_df = load_background_file(bg_file, trait_name, trait_col=args.trait_col)
            
            if len(bg_df) == 0:
                print(f"  Skipping (no genes loaded)", file=sys.stderr)
                continue
            
            cleaned_df = clean_background(bg_df, removed_genes_df, trait_col=args.trait_col)
            
            # Output filename
            output_file = output_dir / f"cleaned_background_{trait_name}.txt"
            write_cleaned_background(cleaned_df, output_file, trait_name)
    
    print("\n=== Background Cleanup Complete ===", file=sys.stderr)


if __name__ == '__main__':
    main()
