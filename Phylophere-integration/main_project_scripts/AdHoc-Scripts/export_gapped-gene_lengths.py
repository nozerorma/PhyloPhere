#!/usr/bin/env python3

import argparse
from Bio import AlignIO
import os
from io import StringIO

def main():
    parser = argparse.ArgumentParser(description='Process MSAs to compute adjusted lengths based on gap positions.')
    parser.add_argument('--species-file', required=True, help='Path to the species list file.')
    parser.add_argument('--alignment-dir', required=True, help='Directory containing the alignment files.')
    parser.add_argument('--alignment-format', default='phylip-relaxed', help='Format of the alignment files (default: phylip-relaxed).')
    parser.add_argument('--output-metadata', required=True, help='Path to the metadata output file.')
    parser.add_argument('--output-summary', required=True, help='Path to the summary output file.')
    parser.add_argument('--output-background', required=True, help='Path to the gene background output file.')
    args = parser.parse_args()

    # Read and deduplicate species list
    species_set = set()
    with open(args.species_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                species_set.add(parts[0])
    species_list = sorted(species_set)
    analyzed_species_str = ','.join(species_list)

    # Open output files
    with open(args.output_metadata, 'w') as metadata_out, \
        open(args.output_summary, 'w') as summary_out, \
        open(args.output_background, 'w') as gene_background_out:

        metadata_out.write("Analyzed species\tGene\tOriginal length\tGap positions count\tAdjusted length\n")
        summary_out.write("Gene\tAdjusted length\n")

        # Process each alignment file in the directory
        for filename in os.listdir(args.alignment_dir):
            filepath = os.path.join(args.alignment_dir, filename)
            if not os.path.isfile(filepath):
                continue

            try:
                # Read entire content to extract header and parse alignment
                with open(filepath, 'r') as f:
                    content = f.read()

                # Extract reference length from the first line header
                first_line = content.split('\n', 1)[0].strip()
                header_parts = first_line.split()
                if len(header_parts) < 2:
                    print(f"Skipping {filename}: malformed header line")
                    continue
                try:
                    ref_length = int(header_parts[1])
                except ValueError:
                    print(f"Skipping {filename}: invalid length in header")
                    continue

                # Parse alignment
                alignment = AlignIO.read(StringIO(content), args.alignment_format)

                # Check if all species are present
                species_in_alignment = {record.id for record in alignment}
                missing_species = [s for s in species_list if s not in species_in_alignment]
                if missing_species:
                    print(f"Skipping {filename}: missing species {', '.join(missing_species)}")
                    continue

                # Collect all gap positions from the species list
                gap_positions = set()
                for record in alignment:
                    if record.id in species_list:
                        for idx, char in enumerate(record.seq):
                            if char == '-':
                                gap_positions.add(idx)

                adjusted_length = ref_length - len(gap_positions)

                # Write results to files
                metadata_out.write(f"{analyzed_species_str}\t{filename}\t{ref_length}\t{len(gap_positions)}\t{adjusted_length}\n")
                summary_out.write(f"{filename}\t{adjusted_length}\n")

                # Write another file for gene background
                gene_name = filename.split('.')[0]
                gene_background_out.write(f"{gene_name}\n")

            except Exception as e:
                print(f"Skipping {filename}: error processing file ({str(e)})")
                continue

if __name__ == "__main__":
    main()
