#!/bin/bash
#SBATCH --job-name=phylophere
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

mkdir -p Slurm

source /etc/profile
module purge

# Try to load Rocky 9 modulepaths and Nextflow version first
if { module load modulepath/UPF/apps 2>/dev/null || module load modulepath/haswell-rocky9 2>/dev/null; } && module load Nextflow/25.10.2 2>/dev/null; then
    echo "[sbatch] Loaded Nextflow/25.10.2 via Rocky 9 stack"
else
    # Fallback to Haswell stack
    module load modulepath/haswell 2>/dev/null
    module load Nextflow/25.04.6 2>/dev/null
    echo "[sbatch] Loaded Nextflow/25.04.6 via Haswell stack"
fi

# Load Conda
module load Miniconda3 2>/dev/null || module load Miniconda3/25.7.0-2 2>/dev/null || true

source ~/.bashrc
conda deactivate
conda activate phylophere

# --- GLOBAL RUN CONFIG ---
export IS_TOY=false
export ALI_WORSKPACE="/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027"
export ALI_DIR="${ALI_WORSKPACE}/PROT/bmge"
export RESUME="${RESUME:-0}"


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

# --- GENERAL AND INPUT DIRS ---
## --- INDIVIDUAL TOOLS ---
export RUN_REPORTING=true
export RUN_PRUNE=true
export RUN_CAAS=false
export RUN_DISAMBIGUATION=false
export RUN_ACCUMULATION=false
export RUN_VEP=false
export RUN_SCORING=false
export RUN_SCORING_STRESS=false
export RUN_ENRICHMENT=false
export RUN_POSENRICH=false
export RUN_RER=false
export RUN_FADE=false

## --- COMMON THINGS MAY BE HARDCODED IN THE CONFIG FILES WHEN NEEDED
export DATADIR="/homes/users/mramon/scratch/2.Primates/1.Primates_data"
export CAAS_OUTBASE="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_traitonly"
export WORK_BASE="/homes/users/mramon/scratch/3.Work_dirs3"
export ASR_CACHE_DIR="${DATADIR}/7.ASR_primates"

export MAP_DIR="${ALI_WORSKPACE}/MAP"
export ENTROPY_DIR="${ALI_WORSKPACE}/PROT_VAR_RAW"

# --- CAAS ---
export TRAIT_FILE="${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"
export TREE_FILE="${DATADIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
export PRUNE_DIR="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP"
export SIMPLE_TRAIT_FILE="${DATADIR}/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv"
export GENE_ENSEMBL_FILE="${DATADIR}/2.Alignments/ensembl_genes.output"
export ALI_SP_NAMES="${DATADIR}/2.Alignments/ali_sp_names.txt"
export INPUT_TAX_ID="${DATADIR}/5.Phylogeny/taxid_species_family_primates.tsv"

# --- OTHER TOOLS ---
## --- ENRICHMENT AND SCORING ---
export GMT_DIR="${DATADIR}/4.External_DBs/GMTs"
export STRING_DB_DIR="${DATADIR}/4.External_DBs/STRING"
export EGG_MEMBERS_FILE="${DATADIR}/4.External_DBs/eggNOG/9443_members.tsv"
export EGG_ANNOTATIONS_FILE="${DATADIR}/4.External_DBs/eggNOG/9443_annotations.tsv"
export DOMAIN_VARIABILITY_FILE="${ALI_WORSKPACE}/reports/variability/domain_variability.tsv"
export UCR_POSITIONS_FILE="${ALI_WORSKPACE}/reports/variability/ucr_positions.tsv"
export FUBAR_SITES_FILE="${ALI_WORSKPACE}/reports/positive_selection/fubar_sites.tsv"
export VEP_PRIMATEAI_DB="${DATADIR}/4.External_DBs/PAI3D/PrimateAI-3D.hg38.txt.gz"
export COSMIC_DB="${DATADIR}/4.External_DBs/COSMIC/Cosmic_MutantCensus_v104_GRCh38.tsv.gz"

## --- RER / FADE ---
export FADE_MODE="all"
export RER_TOOL="build_trait,build_tree,build_matrix,continuous"
export GENE_TREES="${DATADIR}/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"

REPO_DIR="/homes/users/mramon/scratch/0.Phylophere"
SINGLE_RUNNER="${REPO_DIR}/run_phenotype_single_primates_traitonly.sh"

# --- PHENOTYPE CATALOGUE ---
export SLURM_ARRAY_TASK_ID
case $SLURM_ARRAY_TASK_ID in
     1)  CLASS=1; TRAIT="neoplasia_prevalence"; SECONDARY="malignant_prevalence"; NTRAIT="adult_necropsy_count"; CTRAIT="neoplasia_necropsy"; PRUNE="neoplasia_exclude.txt"; PRUNE_SEC="malignant_exclude.txt"; DISCRETE="quintile"
         export SCORING_RER_INPUT="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/neoplasia_prevalence/runtime/filter/rerconverge/rer_results/rerconverge_summary_neoplasia_prevalence.tsv"
         export SCORING_FADE_SUMMARY_TOP="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/neoplasia_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
         export SCORING_FADE_SUMMARY_BOTTOM="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/neoplasia_prevalence/runtime/filter/selection/fade/bottom/fade_summary_bottom.tsv"
         ;;
     2)  CLASS=1; TRAIT="malignant_prevalence"; SECONDARY="neoplasia_prevalence"; NTRAIT="adult_necropsy_count"; CTRAIT="malignant_count"; PRUNE="malignant_exclude.txt"; PRUNE_SEC="neoplasia_exclude.txt"; DISCRETE="quintile"
         export SCORING_RER_INPUT="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/malignant_prevalence/runtime/filter/rerconverge/rer_results/rerconverge_summary_malignant_prevalence.tsv"
         export SCORING_FADE_SUMMARY_TOP="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/malignant_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
         export SCORING_FADE_SUMMARY_BOTTOM="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer/malignant_prevalence/runtime/filter/selection/fade/bottom/fade_summary_bottom.tsv"
         ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo "======================================================"

bash "$SINGLE_RUNNER" "$CLASS" "$TRAIT" "$SECONDARY" "$NTRAIT" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"

