#!/bin/bash
#SBATCH --job-name=bmge_translate
#SBATCH -p haswell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G               # memory per cpu (4G is default for most partitions)
#SBATCH -t 24:00:00
#SBATCH -e /gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Slurm/slurm-%A_%a.err
#SBATCH -o /gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-5%5

# Define directories
CODON_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Mammals_alignments/codon_alignments_phase_3"
PROT_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Mammals_alignments/protein_alignments_phase_3"
TRACK_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Mammals_alignments/protein_trackfiles_phase_3"
TRIM_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/NEOPLASIA_DATA/Mammals_alignments/codon_alignments_trimmed_phase_3"
CONDA_DIR="/gpfs42/robbyfs/homes/users/mramon/.local/bin/"

# Create directories for this portion
mkdir -p "${PROT_DIR}/portion_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TRACK_DIR}/portion_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TRIM_DIR}/portion_${SLURM_ARRAY_TASK_ID}"

# Source Micromamba from the specified directory
baseshell=$(basename $SHELL)
export MAMBA_ROOT_PREFIX=/gpfs42/robbyfs/homes/users/mramon/micromamba
eval "$(/gpfs42/robbyfs/homes/users/mramon/.local/bin/micromamba shell hook --shell $baseshell)"
micromamba activate bmge

# Get the portion index from the job array ID
i=$SLURM_ARRAY_TASK_ID

# Define the function for processing files
do_bmge() {
    file=$1
    PROT=$2
    TRACK=$3
    TRIM=$4
    index=$5
    echo "Processing file: $file, out_prot=$PROT, out_track=$TRACK, out_trim=$TRIM, index=$index"
    gene_name=$(basename $file .fa)
    bmge -i $file -t CODON -g 0.1 -o ${TRIM}/portion_${index}/${gene_name}.phy \
         -oh ${TRACK}/portion_${index}/${gene_name}.html \
         -oaa ${PROT}/portion_${index}/${gene_name}.phy
}
export -f do_bmge

# Find all files in the current portion and process them in parallel
find "${CODON_DIR}/portion_${i}" -name "*.fa" | parallel -j 12 do_bmge {} "$PROT_DIR" "$TRACK_DIR" "$TRIM_DIR" "$i"
