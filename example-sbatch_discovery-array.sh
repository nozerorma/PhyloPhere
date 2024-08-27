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
# File: SBATCH_discovery.sh
#

# Ensure the Slurm log directory exists
mkdir -p Slurm

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash                 
#SBATCH --job-name=name                                             # job name
#SBATCH -p haswell                                                  # partition (In Marvin: normal, haswell)
#SBATCH --nodes=1                                                   # If run in local mode, set to 1. Else, to whatever needed
#SBATCH --ntasks=1                                                  # If run in local mode, set to 1. Else, to whatever needed
#SBATCH --cpus-per-task=12                                          # If run in local mode, set to whatever needed. Else, to 1
#SBATCH --mem-per-cpu=1G                                            # 1GB enough for CAASTOOLS run, might need to increase for other tools
#SBATCH -t 24:00:00                                                 # time (D-HH:MM)
#SBATCH --exclude=mr-00-[01-25],mr-01-[01-02]                       # exclude nodes (opt)
#SBATCH -e Slurm/slurm-%A_%a.err                                    # error file
#SBATCH -o Slurm/slurm-%A_%a.out                                    # output file
#SBATCH --mail-type=START,END,FAIL                                  # notifications for job to be sent
#SBATCH --mail-user=user@provider.edu                               # email to send notifications to
#SBATCH --array=1-2%2                                               # array job (for parallelization)

DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S")                             # add %3N as we want milliseconds too

# Calculate file index
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

## CAASTOOLS DISCOVERY
TRAIT_DIR="/path/goes/here/user/CAAS_traitfiles"                    # Directory where your trait files are located
RESULTS_DIR="/path/goes/here/user/Results/$DATE_WITH_TIME"          # Directory where results will be stored

### Assuming your trait files are named in a way that they can be sorted and indexed
TRAIT_FILES=($(ls $TRAIT_DIR/Scenario/*.tab))                       # Get all trait files for specific scenario
TRAIT_FILE=${TRAIT_FILES[$FILE_INDEX]}                              # Get the trait file for the current array job
TRAIT_NAME=($(basename $TRAIT_FILE))                                # Get the trait name

## CAASTOOLS RESAMPLE
TREE_FILE="/path/goes/here/user/tree.nwk"                           # Path to the tree file
TRAIT_VALUES="/path/goes/here/user/traitfile.tab"                   # Path to the trait values file

## CAASTOOLS BOOTSTRAP
RESAMPLED_DIR="/path/goes/here/user/Resampled_files"                # If resample already done, path to the resampled files

## Nextflow specifics
WORK_DIR="/path/goes/here/user/PhyloPhere-work"                     # Working directory (tmp)
mkdir -p $WORK_DIR

## If running locally (not using Singularity)
#CT_PATH="/binary/location/ct"                                      # Path to the CAASTOOLS binary

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
    --ct_tool "discovery,resample,bootstrap" \
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
EOF
}

# Function to submit the cleanup job
## Needed if running parallel scenarios, as nextflow clean will mess up everything

submit_cleanup_job() {
    local array_job_id=$1
    sbatch --dependency=afterok:$array_job_id <<'EOF'
#!/bin/bash
#SBATCH --job-name=cleanup
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH -e Slurm/cleanup-%j.err
#SBATCH -o Slurm/cleanup-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu

rm -rf /path/goes/here/user/PhyloPhere-work/*

EOF
}

# Main script logic
array_job_id=$(submit_array_job)
submit_cleanup_job $array_job_id
