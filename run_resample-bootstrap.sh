#!/bin/bash
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: example-sbatch_discovery-array.sh
#
# DO NOT MODIFY THE QUEUE UNLESS YOU ARE COMPLETELY SURE ABOUT WHAT YOU ARE DOING




# Define the directory where trait files are located
DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too
TRAIT_DIR="/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/Data/6.CAAS_traitfiles"
TREE_FILE="/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/Data/5.Phylogeny/science.abn7829_data_s4.nex.pruned.tree"
TRAIT_VALUES="/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/Data/10.Bootstrap/traitfile.tab"

# Define workdir
WORK_DIR="PhyloPhere-work"
mkdir -p $WORK_DIR
RESAMPLE_FILE="/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/CAAS-RESULTS-PHYLOPHERE/Ranked/df8_sp.tab/20240708-183505/CT/resample/Ranked/df8_sp.tab/science.abn7829_data_s4.nex.pruned.tree.resampled.output"



# Assuming your trait files are named in a way that they can be sorted and indexed
## Modify according to the needs and folder structure, this assumes traitdir/subdir/traits, for traitdir/traits remove **
TRAIT_FILE=$TRAIT_DIR/Ranked/df4_sp.tab

# Run nextflow pipeline
RESULTS_DIR="CAAS-RESULTS-PHYLOPHERE/Ranked/df8_sp.tab/phylogeny/$DATE_WITH_TIME"
nextflow run main.nf -with-tower -profile normal_local \
    -w $WORK_DIR \
    --ct_tool "resample,bootstrap" \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --maxbgmiss "0" \
    --maxfgmiss "0" \
    --maxmiss "0" \
    --patterns "1,2,3" \
    --tree $TREE_FILE \
    --strategy "BM" \
    --perm_strategy "phylogeny" \
    --template $TRAIT_FILE \
    --traitvalues $TRAIT_VALUES \
    --cycles "1000"
