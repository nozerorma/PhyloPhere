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
# File: discovery-local.sh
#

BASEDIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/caas_new_algorithm"

TRAIT="malignant_prevalence"
WORK_DIR="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/CAAS_2.0/Work"
mkdir -p $WORK_DIR


## CAASTOOLS DISCOVERY
ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Primate_alignments" # Alignment dir in the specified format
TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Discovery/1.Traitfiles/$TRAIT/traitfile.tab" # Directory where your trait files are located

## CAASTOOLS RESAMPLE
TREE_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree" # Path to the tree file
TRAIT_VALUES="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Discovery/1.5.Bootstrap_traitfiles/$TRAIT/boot_traitfile.tab" # Path to the trait values file
CHUNK_SIZE="1000" # Cycles per file (creates directory with multiple resample_*.tab files)

## CAASTOOLS BOOTSTRAP
#RESAMPLED_DIR="/path/goes/here/user/Resampled_files"                # If resample already done, path to the resampled directory
#DISCOVERY_DIR="/path/to/discovery/output"                           # Optional: for position filtering optimization (100-1000x speedup)
#PROGRESS_LOG="bootstrap_progress.log"                                # Optional: progress tracking with timestamps and ETA


## If running locally (not using Singularity)
CT_PATH="./ct"                                      # Path to the CAASTOOLS binary

# Pipeline execution

### If using Singularity, add the -with-singularity flag and the -profile singularity flag
### If using Nextflow Tower, add the -with-tower flag

# Apptainer run
## nextflow run main.nf -with-tower -with-apptainer -profile local,apptainer

# Local run
## nextflow run main.nf -with-tower -profile local

# Slurm localized run (using haswell, with apptainer)
## nextflow run main.nf -with-tower -with-apptainer -profile local,apptainer

# Slurm delocalized run (using haswell, using slurm, with apptainer). Careful with slurm queueing, not very efficient due to overhead.
# Queue should be set in the sbatch directly.
## nextflow run main.nf -with-tower -with-apptainer -profile slurm,apptainer

# Singularity is legacy. Apptainer should be used unless necessary.
# For some reason, Nextflow does not mount HOME directly no more. We need to flag it.
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true
# CAASTOOLS DISCOVERY
echo "Running CAASTOOLS DISCOVERY for trait: $TRAIT"
RESULTS_DIR="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/CAAS_2.0/1.Discovery/$TRAIT"   # Directory where results will be stored

mkdir -p $RESULTS_DIR

nextflow run main.nf -with-tower -profile local \
    -w $WORK_DIR \
    --ct_tool "discovery,resample,bootstrap" \
    --paired_mode \
    --alignment $ALI_DIR \
    --ali_format "phylip-relaxed" \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --maxfgmiss "1" \
    --maxbgmiss "1" \
    --maxgaps "0" \
    --miss_pair \
    --max_conserved "1" \
    --caap_mode \
    --patterns "1,2,3" \
    --tree $TREE_FILE \
    --strategy "BM" \
    --perm_strategy "random" \
    --traitvalues $TRAIT_VALUES \
    --cycles "1000000" \
    --chunk_size "$CHUNK_SIZE"

nextflow clean -f # This flag should be disabled if debugging
