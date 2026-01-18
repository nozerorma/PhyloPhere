#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <msa_file> <target_species> <target_position> <context_range>"
    echo "Example: $0 ZFHX3.phy 'Genus_species1,Genus_species2' 1363 10"
    exit 1
fi

msa_file=$1
target_species=$(echo $2 | tr ',' ' ')
target_position=$3
context_range=$4

# Function to extract and highlight the segment
extract_and_highlight() {
    species_name=$1
    sequence=$2
    target_position=$3
    context_range=$4

    start=$((target_position - context_range - 1))
    end=$((target_position + context_range))
    segment=${sequence:$start:$((end - start))}

    # Highlight the target position
    pos_in_segment=$((target_position - start - 1))
    highlighted_segment="${segment:0:$pos_in_segment}[${segment:$pos_in_segment:1}]${segment:$((pos_in_segment + 1))}"

    echo ">$species_name"
    echo "$highlighted_segment"
}

# Read the MSA file line by line
species=""
sequence=""
while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        if [ -n "$species" ]; then
            for target in $target_species; do
                if [ "$species" == "$target" ]; then
                    extract_and_highlight "$species" "$sequence" "$target_position" "$context_range"
                fi
            done
        fi
        species=$(echo $line | sed 's/^>//')
        sequence=""
    else
        sequence+=$line
    fi
done < "$msa_file"

# Check the last sequence in the file
for target in $target_species; do
    if [ "$species" == "$target" ]; then
        extract_and_highlight "$species" "$sequence" "$target_position" "$context_range"
    fi
done
