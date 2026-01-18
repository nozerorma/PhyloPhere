#!/usr/bin/env python3

import sys
from Bio import AlignIO

def extract_segments(msa_file, target_species, target_position, context_range):
    # Read the MSA file
    alignment = AlignIO.read(msa_file, "phylip")

    # Iterate over each sequence in the alignment
    for record in alignment:
        species_name = record.id
        if species_name in target_species:
            sequence = str(record.seq)
            # Calculate the start and end of the segment
            start = max(0, target_position - context_range - 1)  # Python is 0-based
            end = min(len(sequence), target_position + context_range)
            segment = sequence[start:end]
            # Highlight the target position with brackets
            segment = (segment[:context_range] + 
                       "[" + segment[context_range] + "]" + 
                       segment[context_range+1:])
            print(f">{species_name}\n{segment}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python extract_segments.py <msa_file> <target_species> <target_position> <context_range>")
        print("Example: python extract_segments.py ZFHX3.phy 'Genus_species1,Genus_species2' 1363 10")
        sys.exit(1)

    msa_file = sys.argv[1]
    target_species = sys.argv[2].split(',')
    target_position = int(sys.argv[3])
    context_range = int(sys.argv[4])

    extract_segments(msa_file, target_species, target_position, context_range)
