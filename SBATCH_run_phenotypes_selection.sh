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
# PHYLOPHERE: Selection-Only SLURM Array Submitter
#
# PURPOSE
# -------
# Submits a SLURM array job that re-runs ONLY the selection tools
# (FADE, MoleRate, and optionally RERConverge) across one or more phenotypes,
# using gene sets produced by a prior completed full run.
#
# All CT, signification, disambiguation, postproc, ORA, string, and
# accumulation steps are skipped.
#
# USAGE
# -----
# 1. Set SOURCE_RUN_SUBDIR to the timestamped directory of the completed run
#    from which gene sets should be read.
# 2. Enable/disable tools with RUN_FADE, RUN_MOLERATE, RUN_RER.
# 3. Uncomment the desired phenotype entries in the case statement.
# 4. Adjust --array= range and %1 concurrency limit as needed.
# 5. Run:  bash SBATCH_run_phenotypes_selection.sh
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: SBATCH_run_phenotypes_selection.sh
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
#SBATCH --job-name=phylophere-selection
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH -t 48:00:00
#SBATCH -e Slurm/selection-%A_%a.err
#SBATCH -o Slurm/selection-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-1%1

module load Nextflow
module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate phylophere

# ── Toy mode ─────────────────────────────────────────────────────────────────
# export IS_TOY=true    # small alignments; unset or false for full run

# ── Source run sub-directory ──────────────────────────────────────────────────
# Gene sets are read from:  CAAS_RESULTS/<TRAIT>/<SOURCE_RUN_SUBDIR>/filter/
# Set to the timestamped dir name of a completed full run.
export SOURCE_RUN_SUBDIR="20260226_204433"

# ── Selection tool toggles ────────────────────────────────────────────────────
export RUN_FADE=true
export RUN_MOLERATE=true
export FADE_MODE="gene_set"       # "gene_set" | "all"
export MOLERATE_MODE="gene_set"   # "gene_set" | "all"

# ── RERConverge ───────────────────────────────────────────────────────────────
# Note: gene trees must have species names matching the trait file.
# Keep false until species-name consistency has been verified.
export RUN_RER=false
export RER_TOOL="build_trait,build_tree,build_matrix,continuous"
export RER_GENE_SET_MODE="gene_set"   # "gene_set" | "all"
export GENE_TREES="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"

REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"
SINGLE_RUNNER="${REPO_DIR}/run_phenotype_single_selection.sh"

# ── Phenotype catalogue ───────────────────────────────────────────────────────
# Args passed to run_phenotype_single_selection.sh:
#   CLASS  TRAIT  SECONDARY  CTRAIT  PRUNE  PRUNE_SEC  DISCRETE_METHOD
case $SLURM_ARRAY_TASK_ID in
    1)  CLASS=1; TRAIT="neoplasia_prevalence"; SECONDARY="malignant_prevalence"; CTRAIT="neoplasia_necropsy"; PRUNE="neoplasia_exclude.txt"; PRUNE_SEC="malignant_exclude.txt"; DISCRETE="decile"   ;;
#     2)  CLASS=1; TRAIT="malignant_prevalence"; SECONDARY="neoplasia_prevalence"; CTRAIT="malignant_count";    PRUNE="malignant_exclude.txt"; PRUNE_SEC="neoplasia_exclude.txt"; DISCRETE="decile"   ;;
#     3)  CLASS=2; TRAIT="frug_idx";     SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     4)  CLASS=2; TRAIT="fol_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     5)  CLASS=2; TRAIT="ins_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     6)  CLASS=2; TRAIT="omn_idx";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     7)  CLASS=2; TRAIT="omn_spec_idx"; SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     8)  CLASS=2; TRAIT="herb_idx";     SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#     9)  CLASS=2; TRAIT="folfrug_idx";  SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="decile"   ;;
#    10)  CLASS=2; TRAIT="Ethanol";      SECONDARY=""; CTRAIT=""; PRUNE=""; PRUNE_SEC=""; DISCRETE="quintile" ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo " NODE          : $(hostname)"
echo " MODE          : selection only"
echo "======================================================"

bash "$SINGLE_RUNNER" "$CLASS" "$TRAIT" "$SECONDARY" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"
EOF
}

# Main script logic
array_job_id=$(submit_array_job)
echo "Submitted selection-only array job: ${array_job_id}"

echo ""
echo "Monitor with: squeue -u \$USER"
echo "Logs in     : Slurm/"
