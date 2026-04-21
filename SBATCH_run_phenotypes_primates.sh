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
# File: SBATCH_run_phenotypes.sh
#

# Ensure the Slurm log directory exists
mkdir -p Slurm

# ============================================================
# CONFIGURATION  (edit before submitting)
# ============================================================
REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"
WORK_BASE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs"

# Phenotype catalogue (array 1-10):
#   1.  neoplasia_prevalence  [CLASS 1]
#   2.  malignant_prevalence  [CLASS 1]
#   3.  frug_idx              [CLASS 2]
#   4.  fol_idx               [CLASS 2]
#   5.  ins_idx               [CLASS 2]
#   6.  omn_idx               [CLASS 2]
#   7.  omn_spec_idx          [CLASS 2]
#   8.  herb_idx              [CLASS 2]
#   9.  folfrug_idx           [CLASS 2]
#  10.  Ethanol               [CLASS 2]

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere-cons
#SBATCH --partition=haswell
#SBATCH -t 144:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-2%2

module load modulepath/haswell
module load Nextflow
module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate phylophere

# ── Toy mode ─────────────────────────────────────────────────────────────────
#export IS_TOY=true    # small alignments + 100 cycles; unset or false for full run

# ── Source run sub-directory ──────────────────────────────────────────────────
# Pre-computed CT outputs (discovery, resample, bootstrap, gene sets) are read
# from:  CAAS_RESULTS/<TRAIT>/<SOURCE_RUN_SUBDIR>/filter/
# Set to the timestamped dir name of a completed run, or leave as "runtime" if
# you maintain a stable symlink/alias there.
export SOURCE_RUN_SUBDIR="20260318_201152"

# ── Selection analysis toggles ────────────────────────────────────────────────
# Gene-set paths are auto-derived from SOURCE_RUN_SUBDIR; enable the tools here.
export RUN_FADE=true
export FADE_MODE="gene_set"

# ── RERConverge ───────────────────────────────────────────────────────────────
export RUN_RER=true # Note that we need to fix the gene_trees. They have species which dont match the species in the trait file, and RERConverge fails when that happens. We can either fix the gene trees or modify the RERConverge code to ignore those species. For now, we will set this to false to avoid errors.
export RER_TOOL="build_trait,build_tree,build_matrix,continuous"
export GENE_TREES="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"
export RER_PERM_BATCHES=100
export RER_PERMS_PER_BATCH=100
export RER_GMT_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere/subworkflows/RERCONVERGE/dat/c2.cp.pid.v2026.1.Hs.symbols.gmt"

# ── VEP ───────────────────────────────────────────────────────────────────────
export RUN_VEP=true
export CDS_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons/alignments/Primates_BMGE/TRIM"
export TRACK_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons/alignments/Primates_BMGE/TRACK"
export VEP_REFVERSION="hg38"

# ── CT accumulation knobs ────────────────────────────────────────────────────
export ACCUM_RANDOMIZATION_TYPE="cons_decile"
export ACCUM_FDR=0.1
export ACCUM_LOG_LEVEL="INFO"

# ── Scoring ───────────────────────────────────────────────────────────────────
export RUN_SCORING=true
export RUN_SCORING_ORA=true
export RUN_SCORING_STRESS=true
export SCORING_WINDOW_SIZE_BP=1000000

REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"
SINGLE_RUNNER="${REPO_DIR}/run_phenotype_single_primates.sh"

# ── Phenotype catalogue ───────────────────────────────────────────────────────
# Args passed to run_phenotype_single.sh:
#   CLASS  TRAIT  SECONDARY  CTRAIT  PRUNE  PRUNE_SEC  DISCRETE_METHOD
case $SLURM_ARRAY_TASK_ID in
     1)  CLASS=1; TRAIT="neoplasia_prevalence"; SECONDARY="malignant_prevalence"; CTRAIT="neoplasia_necropsy"; PRUNE="neoplasia_exclude.txt"; PRUNE_SEC="malignant_exclude.txt"; DISCRETE="quintile"   ;;
     2)  CLASS=1; TRAIT="malignant_prevalence"; SECONDARY="neoplasia_prevalence"; CTRAIT="malignant_count";    PRUNE="malignant_exclude.txt"; PRUNE_SEC="neoplasia_exclude.txt"; DISCRETE="quintile"   ;;
#      3)  CLASS=2; TRAIT="frug_idx";     SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      4)  CLASS=2; TRAIT="fol_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      5)  CLASS=2; TRAIT="ins_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      6)  CLASS=2; TRAIT="omn_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      7)  CLASS=2; TRAIT="omn_spec_idx"; SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      8)  CLASS=2; TRAIT="herb_idx";     SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#      9)  CLASS=2; TRAIT="folfrug_idx";  SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     10)  CLASS=2; TRAIT="Ethanol";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="quintile" ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo " NODE          : $(hostname)"
echo "======================================================"

bash "$SINGLE_RUNNER" "$CLASS" "$TRAIT" "$SECONDARY" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"
EOF
}

# # Function to submit the cleanup job
# submit_cleanup_job() {
#     local array_job_id=$1
#     sbatch --dependency=afterok:$array_job_id <<EOF
# #!/bin/bash
# #SBATCH --job-name=phylophere-cleanup
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=2G
# #SBATCH -t 01:00:00
# #SBATCH -e Slurm/cleanup-%j.err
# #SBATCH -o Slurm/cleanup-%j.out
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=miguel.ramon@upf.edu
#
# echo "Cleaning up Nextflow work directories under: ${WORK_BASE}"
# rm -rf "${WORK_BASE:?}/"*
# echo "Cleanup done."
# EOF
# }

# Main script logic
array_job_id=$(submit_array_job)
echo "Submitted array job  : ${array_job_id}  (10 tasks)"
#
# cleanup_job_id=$(submit_cleanup_job "$array_job_id")
# echo "Submitted cleanup job: ${cleanup_job_id}  (depends on afterok:${array_job_id})"

echo ""
echo "Monitor with: squeue -u \$USER"
echo "Logs in     : Slurm/"
