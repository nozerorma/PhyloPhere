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
#SBATCH --job-name=boot-mp
#SBATCH -p haswell
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # node count
#SBATCH --cpus-per-task=8               # only 1 cpu cores is needed to run the nextflow head process
#SBATCH --mem-per-cpu=1G                # memory per cpu (4G is default for most partitions)
#SBATCH -t 12:00:00
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-3%3

# Define the directory where trait files are located
DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too
TRAIT_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/6.CAAS_traitfiles/malignant_prevalence/"
TREE_FILE="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/5.Phylogeny/science.abn7829_data_s4.nex.pruned.tree"
TRAIT_VALUES="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/10.Bootstrap/mp.tab"
ALIGNMENT_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/2.Alignments/*.phy"
CT_PATH="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/0.Phylophere/ct"

WORK_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/PhyloPhere-work2"

# Calculate file index
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Assuming your trait files are named in a way that they can be sorted and indexed
TRAIT_FILES=($(ls $TRAIT_DIR/*.tab))
TRAIT_FILE=${TRAIT_FILES[$FILE_INDEX]}
TRAIT_NAME=($(basename $TRAIT_FILE))

RESULTS_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/Bootstrap/malignant_prevalence/$TRAIT_NAME/$DATE_WITH_TIME"

# Copy ct binary to run folder
# CASE LOCAL
# module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2
# module load DendroPy/4.4.0-foss-2018b-Python-3.6.6
# module load Biopython/1.75-foss-2018b-Python-3.6.6
# module load R/3.6.0-foss-2018b
# module load Java/11.0.2

# Define ct as alias
# module load Singularity/3.7.1-foss-2016b
# module load Nextflow/22.10.5

module load Nextflow

nextflow run main.nf -with-singularity -with-tower -profile singularity,haswell \
    -w $WORK_DIR \
    --ct_tool "resample,bootstrap" \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --alignment "$ALIGNMENT_DIR" \
    --ali_format "phylip-relaxed" \
    --maxbgmiss "0" \
    --maxfgmiss "0" \
    --maxmiss "0" \
    --maxbggaps "0" \
    --maxfggaps "0" \
    --maxgapps "0" \
    --patterns "1,2,3" \
    --tree $TREE_FILE \
    --strategy "BM" \
    --perm_strategy "random" \
    --traitvalues $TRAIT_VALUES \
    --cycles "1000"
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
#SBATCH --mail-user=miguel.ramon@upf.edu

rm -rf /gpfs42/robbyfs/scratch/lab_anavarro/mramon/PhyloPhere-work2/*

EOF
}

# Main script logic
array_job_id=$(submit_array_job)
submit_cleanup_job $array_job_id
