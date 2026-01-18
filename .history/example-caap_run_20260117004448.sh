#!/bin/bash

# Example run script for CAAP mode in PhyloPhere
# This demonstrates how to enable property-based convergence detection

nextflow run main.nf \
    --ct_tool "discovery,resample,bootstrap" \
    --alignment "Data/protein_alignments/**/*" \
    --traitfile "Data/CAAS_traitfiles/traitfile.tab" \
    --ali_format "phylip-relaxed" \
    --patterns "1,2,3" \
    --maxbggaps "0" \
    --maxfggaps "0" \
    --maxgaps "NO" \
    --maxgapsperposition "0.5" \
    --maxbgmiss "1" \
    --maxfgmiss "1" \
    --maxmiss "NO" \
    --paired_mode true \
    --miss_pair true \
    --max_conserved "1" \
    --caap_mode true \
    --tree "Data/Phylogeny/tree.nwk" \
    --traitvalues "Data/CAAS_traitfiles/traitvalues.tab" \
    --strategy "BM" \
    --perm_strategy "random" \
    --cycles "1000" \
    --chunk_size "500" \
    --outdir "Out/Phylophere_CAAP" \
    -profile apptainer \
    -resume

# CAAP Mode Features:
# - Tests 5 grouping schemes (GS0-GS4) for each position
# - GS0 = Identity (classical CAAS, validation)
# - GS1 = Physicochemical properties (6 groups)
# - GS2 = Charge and size (7 groups)
# - GS3 = Hydrophobicity (6 groups)
# - GS4 = Fine-grained biochemical (12 groups)
#
# Output format:
# - Discovery: Gene | Mode | CAAP_Group | Trait | Position | Substitution | ...
# - Bootstrap: Position | CAAP_Group | Count | Cycles | EmpiricalPval
#
# Each position generates one row per matching scheme
