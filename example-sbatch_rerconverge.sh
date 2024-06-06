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

# Function to submit the array job
submit_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere-rerc-mal
#SBATCH -p haswell
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --time=24:00:00

# send mail if needed
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=mail@institution.edu

# Define the directory where trait files are located

DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too
TRAIT_DIR="/path/goes/here/scratch/user/traitdir"
RESULTS_DIR="/path/goes/here/user/outdir/RER/${DATE_WITH_TIME}"
TRAIT_COLUMN="trait_column_name"
SPECIES_COLUMN="species_column_name"

# Define workdir
WORK_DIR="/path/goes/here/user/PhyloPhere-work"
mkdir -p $WORK_DIR

# DIRS
MY_TRAITS="/path/goes/here/user/trait_data.csv"

# Load nextflow
module load Nextflow

# Run rerconverge instance of Phylophere
nextflow run main.nf -with-singularity -with-tower -profile singularity \
    -w $WORK_DIR \
    --outdir $RESULTS_DIR \
    --rer_tool build_trait,build_trees,build_matrix,continuous \
    --my_traits $MY_TRAITS \
    --traitname $TRAIT_NAME \
    --sp_name $SPECIES_COLUMN \
    --rer_minsp 10 \
    --winsorizeRER 3 \
    --winsorizeTrait 3
EOF
}

# Function to submit the cleanup job
submit_cleanup_job() {
    local job_id=$1
    sbatch --dependency=afterok:$job_id <<'EOF'
#!/bin/bash
#SBATCH --job-name=cleanup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p haswell
#SBATCH --mem=2G
#SBATCH -e Slurm/cleanup-%j.err
#SBATCH -o Slurm/cleanup-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mail@institution.edu

rm -rf ${WORK_DIR}/*

EOF
}

# Main script logic
job_id=$(submit_job)
submit_cleanup_job $job_id


