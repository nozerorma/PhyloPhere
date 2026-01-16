#!/bin/sh

mkdir -p Gene-background

for i in $(find . -type f -name "df*"); do
    pair=$(basename ${i} .tab | cut -d "_" -f 3-)
    python export_gapped-gene_lengths.py \
        --species-file $i \
        --alignment-dir ~/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Primate_alignments \
        --output-metadata Gene-background/metadata-${pair}.txt \
        --output-summary Gene-background/summary-${pair}.txt \
        --output-background Gene-background/background-${pair}.txt

done
