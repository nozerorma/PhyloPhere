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
# ─── OVERVIEW ─────────────────────────────────────────────────────────────────
# Submits one SLURM array job per phenotype (one phenotype per node).
# Each array task calls run_phenotype_single.sh with the appropriate
# phenotype parameters.
#
# Phenotype catalogue (10 total):
#   CLASS 1 — Cancer / Pruned-Secondary  (tasks 1-2)
#     1. neoplasia_prevalence
#     2. malignant_prevalence
#   CLASS 2 — Diet / Simple  (tasks 3-10)
#     3. frug_idx        6. omn_idx         9. folfrug_idx
#     4. fol_idx         7. omn_spec_idx   10. Ethanol
#     5. ins_idx         8. herb_idx
#
# Usage:
#   bash SBATCH_run_phenotypes.sh [--toy] [--no-cleanup]
#
#   --toy        : enable toy mode (fewer cycles / randomisations)
#   --no-cleanup : skip the post-run work-directory cleanup job
#
# ─── TOGGLES (edit before submitting) ────────────────────────────────────────
#   RUN_CLASS1   : submit CLASS 1 phenotypes (cancer)
#   RUN_CLASS2   : submit CLASS 2 phenotypes (diet)
#   RUN_FADE     : enable FADE analysis in each task
#   RUN_MOLERATE : enable MoleRate analysis in each task
# ─────────────────────────────────────────────────────────────────────────────

set -Eeuo pipefail

# Ensure the Slurm log directory exists
mkdir -p Slurm

# ============================================================
# COMMAND-LINE FLAGS
# ============================================================
IS_TOY=false
RUN_CLEANUP=true

for arg in "$@"; do
    case "$arg" in
        --toy)        IS_TOY=true     ;;
        --no-cleanup) RUN_CLEANUP=false ;;
        *)
            echo "Unknown argument: $arg"
            echo "Usage: $0 [--toy] [--no-cleanup]"
            exit 1
            ;;
    esac
done

# ============================================================
# USER TOGGLES
# ============================================================
RUN_CLASS1=true    # neoplasia_prevalence | malignant_prevalence
RUN_CLASS2=true    # diet / dietary-index phenotypes (María Sánchez Bermúdez)
RUN_FADE=false     # HyPhy FADE (directional amino-acid selection)
RUN_MOLERATE=false # HyPhy MoleRate (relative evolutionary rate)

# ============================================================
# CLUSTER PATH CONFIGURATION
# These must match the paths set in run_phenotype_single.sh
# ============================================================
REPO_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/PhyloPhere"
WORK_BASE="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/PhyloPhere-work"
SINGLE_RUNNER="${REPO_DIR}/run_phenotype_single.sh"

if [ ! -f "$SINGLE_RUNNER" ]; then
    echo "Error: run_phenotype_single.sh not found at: $SINGLE_RUNNER"
    echo "  Make sure REPO_DIR is set correctly above."
    exit 1
fi

# ============================================================
# PHENOTYPE CATALOGUE
# ─────────────────────────────────────────────────────────────
# Parallel arrays — one slot per phenotype.
# Index 0-based; SLURM_ARRAY_TASK_ID is 1-based.
#
# Columns:
#   PHENO_CLASS           : 1 or 2
#   PHENO_TRAIT           : trait name (column in CSV)
#   PHENO_SECONDARY       : secondary trait  (CLASS 1 only; "" for CLASS 2)
#   PHENO_CTRAIT          : c_trait          (CLASS 1 only; "" for CLASS 2)
#   PHENO_PRUNE           : prune file basename (CLASS 1 only; "")
#   PHENO_PRUNE_SEC       : secondary prune file basename (CLASS 1 only; "")
#   PHENO_DISCRETE_METHOD : discretisation method (CLASS 2 only; "decile")
# ============================================================

PHENO_CLASS=()
PHENO_TRAIT=()
PHENO_SECONDARY=()
PHENO_CTRAIT=()
PHENO_PRUNE=()
PHENO_PRUNE_SEC=()
PHENO_DISCRETE_METHOD=()

# Helper: add a CLASS 1 (pruned-secondary) phenotype
add_class1() {
    local trait="$1" sec="$2" ctrait="$3" prune="$4" prune_sec="$5"
    PHENO_CLASS+=("1")
    PHENO_TRAIT+=("$trait")
    PHENO_SECONDARY+=("$sec")
    PHENO_CTRAIT+=("$ctrait")
    PHENO_PRUNE+=("$prune")
    PHENO_PRUNE_SEC+=("$prune_sec")
    PHENO_DISCRETE_METHOD+=("decile")
}

# Helper: add a CLASS 2 (simple) phenotype
add_class2() {
    local trait="$1" discrete_method="${2:-decile}"
    PHENO_CLASS+=("2")
    PHENO_TRAIT+=("$trait")
    PHENO_SECONDARY+=("")
    PHENO_CTRAIT+=("")
    PHENO_PRUNE+=("")
    PHENO_PRUNE_SEC+=("")
    PHENO_DISCRETE_METHOD+=("$discrete_method")
}

# ── CLASS 1 — Cancer phenotypes (Comparative Oncology Alliance · IBE-UPF) ────
if [ "$RUN_CLASS1" = true ]; then
    add_class1 "neoplasia_prevalence" "malignant_prevalence"  "neoplasia_necropsy" "neoplasia_exclude.txt" "malignant_exclude.txt"
    add_class1 "malignant_prevalence" "neoplasia_prevalence"  "malignant_count"    "malignant_exclude.txt" "neoplasia_exclude.txt"
fi

# ── CLASS 2 — Dietary phenotypes (María Sánchez Bermúdez · UOC / IBE-UPF) ───
if [ "$RUN_CLASS2" = true ]; then
    # Trophic specialisation indices
    add_class2 "frug_idx"     "decile"
    add_class2 "fol_idx"      "decile"
    add_class2 "ins_idx"      "decile"
    add_class2 "omn_idx"      "decile"
    add_class2 "omn_spec_idx" "decile"
    add_class2 "herb_idx"     "decile"
    add_class2 "folfrug_idx"  "decile"
    # Ethanol: fewer categories → quintile to avoid empty bins
    add_class2 "Ethanol"      "quintile"
fi

N_PHENOTYPES="${#PHENO_TRAIT[@]}"

if [ "$N_PHENOTYPES" -eq 0 ]; then
    echo "Error: no phenotypes selected (both RUN_CLASS1 and RUN_CLASS2 are false?)."
    exit 1
fi

# ============================================================
# SUMMARY
# ============================================================
echo "=========================================="
echo " PHYLOPHERE SBATCH MULTI-PHENOTYPE RUNNER"
echo "=========================================="
echo " IS_TOY       : $IS_TOY"
echo " RUN_CLASS1   : $RUN_CLASS1"
echo " RUN_CLASS2   : $RUN_CLASS2"
echo " RUN_FADE     : $RUN_FADE"
echo " RUN_MOLERATE : $RUN_MOLERATE"
echo " N_PHENOTYPES : $N_PHENOTYPES"
echo " Array range  : 1-${N_PHENOTYPES}%${N_PHENOTYPES}"
echo "------------------------------------------"
printf " %-3s %-2s  %-35s %s\n" "IDX" "CL" "TRAIT" "DISCRETE/SECONDARY"
echo " ---  --  -----------------------------------  ------------------"
for i in $(seq 0 $((N_PHENOTYPES - 1))); do
    extra="${PHENO_SECONDARY[$i]:-${PHENO_DISCRETE_METHOD[$i]}}"
    printf " %-3d %-2s  %-35s %s\n" "$((i+1))" "${PHENO_CLASS[$i]}" "${PHENO_TRAIT[$i]}" "$extra"
done
echo "=========================================="
echo ""

# ============================================================
# ARRAY JOB SUBMISSION
# ============================================================
submit_array_job() {
    sbatch --parsable \
        --array="1-${N_PHENOTYPES}%${N_PHENOTYPES}" \
        --export=ALL,IS_TOY="$IS_TOY",RUN_FADE="$RUN_FADE",RUN_MOLERATE="$RUN_MOLERATE",\
PHENO_CLASSES="$(IFS=','; echo "${PHENO_CLASS[*]}")",\
PHENO_TRAITS="$(IFS=','; echo "${PHENO_TRAIT[*]}")",\
PHENO_SECONDARIES="$(IFS=','; echo "${PHENO_SECONDARY[*]}")",\
PHENO_CTRAITS="$(IFS=','; echo "${PHENO_CTRAIT[*]}")",\
PHENO_PRUNES="$(IFS=','; echo "${PHENO_PRUNE[*]}")",\
PHENO_PRUNE_SECS="$(IFS=','; echo "${PHENO_PRUNE_SEC[*]}")",\
PHENO_DISCRETE_METHODS="$(IFS=','; echo "${PHENO_DISCRETE_METHOD[*]}")",\
SINGLE_RUNNER="$SINGLE_RUNNER" \
        <<'HEREDOC'
#!/bin/bash
#SBATCH --job-name=phylophere-pheno
#SBATCH -p haswell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t 48:00:00
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu

module load Nextflow

# ── Reconstruct arrays from comma-separated env vars ─────────────────────────
IFS=',' read -ra _CLASSES          <<< "$PHENO_CLASSES"
IFS=',' read -ra _TRAITS           <<< "$PHENO_TRAITS"
IFS=',' read -ra _SECONDARIES      <<< "$PHENO_SECONDARIES"
IFS=',' read -ra _CTRAITS          <<< "$PHENO_CTRAITS"
IFS=',' read -ra _PRUNES           <<< "$PHENO_PRUNES"
IFS=',' read -ra _PRUNE_SECS       <<< "$PHENO_PRUNE_SECS"
IFS=',' read -ra _DISCRETE_METHODS <<< "$PHENO_DISCRETE_METHODS"

# ── Select this task's phenotype (SLURM arrays are 1-based) ──────────────────
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))

CLASS="${_CLASSES[$IDX]}"
TRAIT="${_TRAITS[$IDX]}"
SECONDARY="${_SECONDARIES[$IDX]:-}"
CTRAIT="${_CTRAITS[$IDX]:-}"
PRUNE="${_PRUNES[$IDX]:-}"
PRUNE_SEC="${_PRUNE_SECS[$IDX]:-}"
DISCRETE="${_DISCRETE_METHODS[$IDX]:-decile}"

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID  (0-index: $IDX)"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo " NODE          : $(hostname)"
echo "======================================================"

# ── Export per-run env flags to single-runner ─────────────────────────────────
export IS_TOY RUN_FADE RUN_MOLERATE

bash "$SINGLE_RUNNER" \
    "$CLASS" "$TRAIT" "$SECONDARY" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"

HEREDOC
}

# ============================================================
# CLEANUP JOB (runs after all array tasks succeed)
# ============================================================
submit_cleanup_job() {
    local array_job_id="$1"
    sbatch --dependency=afterok:"$array_job_id" <<HEREDOC
#!/bin/bash
#SBATCH --job-name=phylophere-cleanup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p haswell
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH -e Slurm/cleanup-%j.err
#SBATCH -o Slurm/cleanup-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu

echo "Cleaning up Nextflow work directories under: ${WORK_BASE}"
rm -rf "${WORK_BASE:?}/"*
echo "Cleanup done."
HEREDOC
}

# ============================================================
# SUBMIT
# ============================================================
array_job_id=$(submit_array_job)
echo "Submitted array job: ${array_job_id}  (${N_PHENOTYPES} tasks)"

if [ "$RUN_CLEANUP" = true ]; then
    cleanup_job_id=$(submit_cleanup_job "$array_job_id")
    echo "Submitted cleanup job: ${cleanup_job_id}  (depends on afterok:${array_job_id})"
else
    echo "Cleanup skipped (--no-cleanup)."
fi

echo ""
echo "Monitor with: squeue -u \$USER"
echo "Logs in     : Slurm/"
