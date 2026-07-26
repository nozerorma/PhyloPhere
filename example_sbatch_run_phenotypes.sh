#!/bin/bash
#SBATCH --job-name=phylophere
#SBATCH --partition=standard
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --array=1-2%2

mkdir -p Slurm

source /etc/profile
module purge 2>/dev/null || true

# Load Nextflow and Conda modules if using cluster modules
module load Nextflow 2>/dev/null || true
module load Miniconda3 2>/dev/null || true

source ~/.bashrc
conda deactivate 2>/dev/null || true
conda activate phylophere

# --- GLOBAL RUN CONFIG ---
export IS_TOY=false
export ALI_WORSKPACE="/path/to/alignments_workspace"
export ALI_DIR="${ALI_WORSKPACE}/PROT/bmge"
export RESUME="${RESUME:-1}"

if [ "$IS_TOY" = true ]; then
    export TAG="_toy"
    export TOY_MODE=true
    export TOY_N="1000"
    export PERMS_CYCLES="1000"
    export CAAS_FULL_PERMS="100"
    export N_RANDOMIZATIONS="1000"
    export RER_PERM_BATCHES="10"
    export RER_PERMS_PER_BATCH="10"
else
    export TAG=""
    export TOY_MODE=false
    export TOY_N=""
    export PERMS_CYCLES="1000000"
    export CAAS_FULL_PERMS="${CAAS_FULL_PERMS:-1000}"
    export N_RANDOMIZATIONS="1000000"
    export RER_PERM_BATCHES="10"
    export RER_PERMS_PER_BATCH="100"
fi

# --- INDIVIDUAL TOOL TOGGLES ---
export RUN_REPORTING=true
export RUN_PRUNE=true
export RUN_CAAS=true
export RUN_DISAMBIGUATION=true
export RUN_ACCUMULATION=true
export RUN_VEP=true
export RUN_SCORING=true
export RUN_SCORING_STRESS=true
export RUN_ENRICHMENT=true
export RUN_POSENRICH=true
export RUN_RER=false
export RUN_FADE=false

# --- GENERAL AND INPUT DIRS ---
export DATADIR="/path/to/data"
export CAAS_OUTBASE="/path/to/results"
export WORK_BASE="/path/to/work_dirs"
export ASR_CACHE_DIR="${DATADIR}/ASR_cache"

export MAP_DIR="${ALI_WORSKPACE}/MAP"
export ENTROPY_DIR="${ALI_WORSKPACE}/PROT_VAR_RAW"

# --- CAAS & TRAIT INPUTS ---
export TRAIT_FILE="${DATADIR}/traits.csv"
export TREE_FILE="${DATADIR}/phylogeny.tree"
export PRUNE_DIR="${DATADIR}/prune_lists"
export SIMPLE_TRAIT_FILE="${DATADIR}/simple_traits.csv"
export GENE_ENSEMBL_FILE="${DATADIR}/ensembl_genes.tsv"
export ALI_SP_NAMES="${DATADIR}/ali_sp_names.txt"
export INPUT_TAX_ID="${DATADIR}/taxid_species.tsv"

# --- EXTERNAL DATABASES ---
export GMT_DIR="${DATADIR}/GMTs"
export STRING_DB_DIR="${DATADIR}/STRING"
export EGG_MEMBERS_FILE="${DATADIR}/eggNOG/members.tsv"
export EGG_ANNOTATIONS_FILE="${DATADIR}/eggNOG/annotations.tsv"
export DOMAIN_VARIABILITY_FILE="${ALI_WORSKPACE}/reports/variability/domain_variability.tsv"
export UCR_POSITIONS_FILE="${ALI_WORSKPACE}/reports/variability/ucr_positions.tsv"
export FUBAR_SITES_FILE="${ALI_WORSKPACE}/reports/positive_selection/fubar_sites.tsv"
export VEP_PRIMATEAI_DB="${DATADIR}/PAI3D/PrimateAI-3D.hg38.txt.gz"
export COSMIC_DB="${DATADIR}/COSMIC/Cosmic_MutantCensus.tsv.gz"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINGLE_RUNNER="${REPO_DIR}/example_run_phenotype_single.sh"

# --- PHENOTYPE CATALOGUE ---
export SLURM_ARRAY_TASK_ID
case $SLURM_ARRAY_TASK_ID in
     1)  CLASS=1; TRAIT="phenotype_1"; SECONDARY="phenotype_2"; NTRAIT="n_count"; CTRAIT="c_count"; PRUNE="prune_1.txt"; PRUNE_SEC="prune_2.txt"; DISCRETE="quintile"
         ;;
     2)  CLASS=1; TRAIT="phenotype_2"; SECONDARY="phenotype_1"; NTRAIT="n_count"; CTRAIT="c_count"; PRUNE="prune_2.txt"; PRUNE_SEC="prune_1.txt"; DISCRETE="quintile"
         ;;
     *)
         echo "Error: Invalid task ID $SLURM_ARRAY_TASK_ID" >&2
         exit 1
         ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo "======================================================"

bash "$SINGLE_RUNNER" "$CLASS" "$TRAIT" "$SECONDARY" "$NTRAIT" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"
