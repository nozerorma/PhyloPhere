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

# Ensure the Slurm log directory exists
mkdir -p Slurm

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere-discovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH -p haswell
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=mail@institution.edu
#SBATCH --array=1-8%8

# Define the directory where trait files are located
DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too
TRAIT_DIR="/path/goes/here/scratch/user/traitdir"
RESULTS_DIR="/path/goes/here/user/outdir/CT/${DATE_WITH_TIME}"

# Define workdir
WORK_DIR="/path/goes/here/user/PhyloPhere-work"
mkdir -p $WORK_DIR

# Calculate file index
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Assuming your trait files are named in a way that they can be sorted and indexed
## Modify according to the needs and folder structure, this assumes traitdir/subdir/traits, for traitdir/traits remove **
TRAIT_FILES=($(ls $TRAIT_DIR/**/*.tab))
TRAIT_FILE=${TRAIT_FILES[$FILE_INDEX]}

# Load nextflow module
module load Nextflow

# Run nextflow pipeline
nextflow run main.nf -with-singularity -with-tower -profile singularity \
    -w $WORK_DIR \
    --ct_tool discovery \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --maxbgmiss "0" \
    --maxfgmiss "0" \
    --maxmiss "0" \
    --patterns "1,2,3"
EOF
}

# Function to submit the cleanup job
submit_cleanup_job() {
    local array_job_id=$1
    sbatch --dependency=afterok:$array_job_id <<'EOF'
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
array_job_id=$(submit_array_job)
submit_cleanup_job $array_job_id
