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
# File: example-sbatch_rerconverge.sh
#
# DO NOT MODIFY THE QUEUE UNLESS YOU ARE COMPLETELY SURE ABOUT WHAT YOU ARE DOING

# Ensure the Slurm log directory exists
mkdir -p Slurm

#SBATCH --job-name=name                                             # job name
#SBATCH -p haswell                                                  # partition (In Marvin: normal, haswell)
#SBATCH --nodes=1                                                   # If run in local mode, set to 1. Else, to whatever needed
#SBATCH --ntasks=1                                                  # If run in local mode, set to 1. Else, to whatever needed
#SBATCH --cpus-per-task=8                                           # If run in local mode, set to whatever needed. Else, to 1
#SBATCH --mem=16G                                                   # Resource hungry, 16GB should be fine
#SBATCH -t 24:00:00                                                 # time (D-HH:MM)
#SBATCH --exclude=mr-00-[01-25],mr-01-[01-02]                       # exclude nodes (opt)
#SBATCH -e Slurm/slurm-%A_%a.err                                    # error file
#SBATCH -o Slurm/slurm-%A_%a.out                                    # output file
#SBATCH --mail-type=START,END,FAIL                                  # notifications for job to be sent
#SBATCH --mail-user=user@provider.edu                               # email to send notifications to

DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too

# RERConverge specifics
MY_TRAITS="/path/goes/here/user/trait_data.csv"
RESULTS_DIR="/path/goes/here/user/outdir/RER/${DATE_WITH_TIME}"
TRAIT_COLUMN="trait_column_name"
SPECIES_COLUMN="species_column_name"

## Nextflow specifics
WORK_DIR="/path/goes/here/user/PhyloPhere-work"                     # Working directory (tmp)
mkdir -p $WORK_DIR

## If running in Marvin Normal partition                            # Load the modules needed for run in normal partition (older arch)
#module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2
#module load DendroPy/4.4.0-foss-2018b-Python-3.6.6
#module load Biopython/1.75-foss-2018b-Python-3.6.6
#module load R/3.6.0-foss-2018b
#module load Java/11.0.2
#module load Singularity/3.7.1-foss-2016b
#module load Nextflow/22.10.5

## Else if running in Marvin Haswell partition                      # Load the modules needed for run in haswell partition (newer arch)
module load Nextflow

# Pipeline execution

### If using Singularity, add the -with-singularity flag and the -profile singularity flag
### If using Nextflow Tower, add the -with-tower flag
### If running in Marvin normal partition, add the -profile normal flag
### If running in Marvin haswell partition, add the -profile haswell flag
##### IF WE WANT TO USE SLURM AS EXECUTOR, ADD THE -profile normal_slurm OR -profile haswell_slurm FLAG AND MODIFY SBATCH PARAMETERS ACCORDINGLY

nextflow run main.nf -with-singularity -with-tower -profile singularity,haswell \
    -w $WORK_DIR \
    --outdir $RESULTS_DIR \
    --rer_tool "build_trait,build_trees,build_matrix,continuous" \
    --my_traits $MY_TRAITS \
    --traitname $TRAIT_NAME \
    --sp_name $SPECIES_COLUMN \
    --rer_minsp 10 \
    --winsorizeRER 3 \
    --winsorizeTrait 3

# If the tool has already been run, other parameters can be used to rerun the tool
## ie. --trait_out, --trees_out

nextflow clean -f

