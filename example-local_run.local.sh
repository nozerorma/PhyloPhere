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

BASEDIR="/media/miguel/TMB_XPrize4/caas"

TRAIT="malignant_prevalence"
WORK_DIR="/tmp/Results/work"
mkdir -p $WORK_DIR


## CAASTOOLS DISCOVERY
ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Ali_toy" # Alignment dir in the specified format
TRAIT_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Discovery/1.Traitfiles/$TRAIT/traitfiles" # Directory where your trait files are located

## CAASTOOLS RESAMPLE
TREE_FILE="$BASEDIR/Data/5.Phylogeny" # Path to the tree file
TRAIT_VALUES="$BASEDIR/Out/2.CAAS/1.Discovery/1.5.Bootstrap_traitfiles/$TRAIT/traitfile.tab" # Path to the trait values file

## CAASTOOLS BOOTSTRAP
#RESAMPLED_DIR="/path/goes/here/user/Resampled_files"                # If resample already done, path to the resampled files


## If running locally (not using Singularity)
CT_PATH="./ct"                                      # Path to the CAASTOOLS binary

# Pipeline execution

### If using Singularity, add the -with-singularity flag and the -profile singularity flag
### If using Nextflow Tower, add the -with-tower flag
### If running in Marvin normal partition, add the -profile normal flag
### If running in Marvin haswell partition, add the -profile haswell flag
##### IF WE WANT TO USE SLURM AS EXECUTOR, ADD THE -profile normal_slurm OR -profile haswell_slurm FLAG AND MODIFY SBATCH PARAMETERS ACCORDINGLY

#nextflow run main.nf -with-singularity -with-tower -profile singularity,haswell \

# Singularity is legacy. Use apptainer for recent Nextflow version.
# For some reason, Nextflow does not mount HOME directly no more. We need to flag it.

export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true


for TRAIT_FILE in $(find $TRAIT_DIR -type f -name "*1.tab"); do

    TRAIT_FILENAME=($(basename $TRAIT_FILE))
    echo "Working with traitfile $TRAIT_FILENAME"

    RESULTS_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Discovery/toy/$TRAIT/$TRAIT_FILENAME/"   # Directory where results will be stored

    mkdir -p $RESULTS_DIR

    nextflow run main.nf -with-tower -with-apptainer -profile apptainer,local \
        -w $WORK_DIR \
        --ct_tool "discovery" \
        --alignment $ALI_DIR \
        --ali_format "phylip-relaxed" \
        --traitfile $TRAIT_FILE \
        --outdir $RESULTS_DIR \
        --maxbgmiss "0" \
        --maxfgmiss "0" \
        --maxmiss "0" \
        --patterns "1,2,3" \
        --tree $TREE_FILE \
        --strategy "BM" \
        --perm_strategy "random" \
        --traitvalues $TRAIT_VALUES \
        --cycles "1000"

    nextflow clean -f

done
