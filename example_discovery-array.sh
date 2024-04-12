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
# File: SBATCH_discovery-array.sh
#

# change according to necessities and cluster characteristics

#SBATCH --job-name=nfct-discovery-array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p haswell
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=user@mail.com
#SBATCH --array=1-10%5

# Define the directory where trait files are located
DATE_WITH_TIME=`date "+%Y%m%d-%H%M%S"`

TRAIT_DIR="/my/trait/dir"
RESULTS_DIR="/my/result/dir/${DATE_WITH_TIME}"
TOWER_TOKEN="mytowertoken"

#WORK_DIR="/if/you/need/to/specify/different/work/dir"

# Calculate file index
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Assuming your trait files are named in a way that they can be sorted and indexed, in a two level structure
TRAIT_FILES=($(ls $TRAIT_DIR/**/*.tab))
TRAIT_FILE=${TRAIT_FILES[$FILE_INDEX]}

module load Nextflow
nextflow run main.nf -with-singularity -with-tower -profile singularity \
    --accessToken $TOWER_TOKEN \
    --ct_tool discovery \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --maxbgmiss "2" \
    --maxfgmiss "2"  \
    --maxmiss "2" \
    --patterns "1,2,3,4"

