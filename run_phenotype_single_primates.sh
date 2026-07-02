#!/bin/bash
set -Eeuo pipefail

# --- ARGUMENTS ---
CLASS="${1:?Usage: $0 <CLASS> <TRAIT> [SECONDARY_TRAIT] [C_TRAIT] [PRUNE_FILE] [PRUNE_SECONDARY_FILE] [DISCRETE_METHOD]}"
TRAIT="${2:?TRAIT must be provided}"
SECONDARY_TRAIT="${3:-}"
C_TRAIT="${4:-}"
PRUNE_FILE="${5:-}"
PRUNE_SECONDARY_FILE="${6:-}"
DISCRETE_METHOD="${7:-}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp=$(date +%Y%m%d_%H%M%S)

# --- CLUSTER / ENVIRONMENT CONFIGURATION ---
DATADIR="/homes/users/mramon/scratch/2.Primates/1.Primates_data"
CAAS_OUTBASE="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final"
WORK_BASE="/homes/users/mramon/scratch/3.Work_dirs_new"
ASR_CACHE_DIR="${ASR_CACHE_DIR:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/7.ASR_primates}"

TRAIT_FILE="${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"
TREE_FILE="${DATADIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
PRUNE_DIR="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP"
SIMPLE_TRAIT_FILE="${DATADIR}/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv"
ALI_SP_NAMES="${ALI_SP_NAMES:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/2.Alignments/ali_sp_names.txt}"
INPUT_TAX_ID="${INPUT_TAX_ID:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/5.Phylogeny/taxid_species_family_primates.tsv}"


# --- ENV ACTIVATION ---
source ~/.bashrc
conda deactivate || true
conda activate phylophere

# Setup unique Nextflow home directory to avoid file lock collisions on history database
export NXF_HOME="$WORK_BASE/.nextflow/.nextflow_${TRAIT}"
mkdir -p "$NXF_HOME"
ln -sfn /homes/users/mramon/.nextflow/plugins "$NXF_HOME/plugins"

# Create and cd into a unique run directory to prevent Nextflow session/cache database lock collisions
RUN_DIR="$WORK_BASE/${TRAIT}_nxf_run"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

# --- TOY / FULL MODE CONFIG ---
IS_TOY="${IS_TOY:-false}"
if [ "$IS_TOY" = true ]; then
    TAG="_toy"
    ALI_DIR="${DATADIR}/2.Alignments/Ali_toy"
    CYCLES="1000"
    N_RANDOMIZATIONS="1000"
    CAAS_FULL_PERMS="100"
else
    TAG=""
    ALI_DIR="/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/PROT/bmge"
    CYCLES="1000000"
    N_RANDOMIZATIONS="1000000"
    CAAS_FULL_PERMS="${CAAS_FULL_PERMS:-1000}"
fi

RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${timestamp}"
WORK_DIR="${WORK_BASE}/${TRAIT}${TAG}/${timestamp}"
mkdir -p "${RESULTS_BASE}" "${WORK_DIR}"

# --- PRECOMPUTED INPUT DETECTION ---
# Only pass precomputed directories/files to Nextflow if they actually exist on disk.
# This prevents Nextflow from initializing with empty channels and crashing downstream.
PRECOMPUTED_FLAGS=()
if [ -n "${SOURCE_RUN_SUBDIR:-}" ]; then
    SOURCE_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${SOURCE_RUN_SUBDIR}"
    if [ -f "${SOURCE_BASE}/caastools/discovery.tab" ]; then
        PRECOMPUTED_FLAGS+=(--discovery_out "${SOURCE_BASE}/caastools/discovery.tab")
    fi
    if [ -d "${SOURCE_BASE}/resample/nw_tree.resampled.output" ] || [ -f "${SOURCE_BASE}/resample/nw_tree.resampled.output" ]; then
        PRECOMPUTED_FLAGS+=(--resample_out "${SOURCE_BASE}/resample/nw_tree.resampled.output")
    fi
    if [ -f "${SOURCE_BASE}/caastools/bootstrap.tab" ]; then
        PRECOMPUTED_FLAGS+=(--bootstrap_input "${SOURCE_BASE}/caastools/bootstrap.tab")
    fi
    if [ -f "${SOURCE_BASE}/caastools/background_genes.output" ]; then
        PRECOMPUTED_FLAGS+=(--background_input "${SOURCE_BASE}/caastools/background_genes.output")
    fi
    if [ -f "${SOURCE_BASE}/data_exploration/2.CT/1.Traitfiles/traitfile.tab" ]; then
        PRECOMPUTED_FLAGS+=(--caas_config "${SOURCE_BASE}/data_exploration/2.CT/1.Traitfiles/traitfile.tab")
    fi
fi

# --- NEXTFLOW FLAGS ---
NF_FLAGS=(
    -profile slurm
    -with-tower
    -name "${TRAIT}_$(cat /proc/sys/kernel/random/uuid 2>/dev/null | tr -d '-' | cut -c1-8 || printf '%08x' $RANDOM)"
    --alignment "$ALI_DIR"
    --ali_format "fasta"
    --tree "$TREE_FILE"
    --cycles "$CYCLES"
    --accumulation_n_randomizations "$N_RANDOMIZATIONS"
    --accumulation_randomization_type "${ACCUM_RANDOMIZATION_TYPE:-cons_decile}"
    --accumulation_fdr "${ACCUM_FDR:-0.05}"
    --accumulation_log_level "${ACCUM_LOG_LEVEL:-INFO}"
    --ct_disambig_asr_mode "precomputed"
    --ct_disambig_asr_cache_dir "$ASR_CACHE_DIR"
    --ct_tool "discovery,resample,bootstrap"
    --ct_signification
    --ct_disambiguation
    --ct_postproc
    --ct_accumulation
    --caas_postproc_mode "filter"
    --outdir "$RESULTS_BASE"
    -resume
    "${PRECOMPUTED_FLAGS[@]+"${PRECOMPUTED_FLAGS[@]}"}"
)

# RERConverge Flags
if [ "${RUN_RER:-true}" = true ]; then
    NF_FLAGS+=(
        --rer_tool "${RER_TOOL:-build_trait,build_tree,build_matrix,continuous}"
        --gene_trees "${GENE_TREES:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt}"
        --rer_perm_batches "${RER_PERM_BATCHES:-10}"
        --rer_perms_per_batch "${RER_PERMS_PER_BATCH:-100}"
    )
fi

# FADE Flags
if [ "${RUN_FADE:-true}" = true ]; then
    NF_FLAGS+=(
        --fade
        --fade_mode "${FADE_MODE:-all}"
    )
fi

# VEP Flags
if [ "${RUN_VEP:-true}" = true ]; then
    NF_FLAGS+=(
        --vep
        --vep_map_dir "${MAP_DIR:-/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/MAP}"
    )
fi

# POSENRICH (position-wise enrichment) flags
if [ "${RUN_POSENRICH:-true}" = true ]; then
    NF_FLAGS+=(
        --posenrich
        --egg_members_file "${EGG_MEMBERS_FILE:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/eggNOG/9443_members.tsv}"
        --egg_annotations_file "${EGG_ANNOTATIONS_FILE:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/eggNOG/9443_annotations.tsv}"
        --cosmic_db "${COSMIC_DB:-${REPO_DIR}/to_integrate/Cosmic_MutantCensus_v104_GRCh38.tsv.gz}"
    )
fi

# Enrichment Flags
if [ "${RUN_ENRICHMENT:-true}" = true ]; then
    NF_FLAGS+=(
        --enrichment
        --gmt_dir "${GMT_DIR:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/GMTs}"
    )
fi

# Scoring Flags
if [ "${RUN_SCORING:-true}" = true ]; then
    NF_FLAGS+=(
        --scoring
        --scoring_window_size_bp "${SCORING_WINDOW_SIZE_BP:-1000000}"
        --gene_ensembl_file "${GENE_ENSEMBL_FILE:-/homes/users/mramon/scratch/2.Primates/1.Primates_data/2.Alignments/ensembl_genes.output}"
    )
    if [ "${RUN_SCORING_STRESS:-true}" = true ]; then
        NF_FLAGS+=(--scoring_stress)
    fi
fi

# CAAS Permulation Flags
if [ "${RUN_CAAS_PERMULATION:-true}" = true ]; then
    NF_FLAGS+=(
        --caas_permulation_enrichment true
        --caas_full_perms "$CAAS_FULL_PERMS"
    )
fi

# Species and Tax ID Flags
[ -n "${ALI_SP_NAMES:-}" ] && NF_FLAGS+=(--ali_sp_names "$ALI_SP_NAMES")
[ -n "${INPUT_TAX_ID:-}" ] && NF_FLAGS+=(--tax_id "$INPUT_TAX_ID")

# Precomputed FADE/RER summaries for scoring
[ -n "${SCORING_RER_INPUT:-}" ] && NF_FLAGS+=(--scoring_rer_input "$SCORING_RER_INPUT")
[ -n "${SCORING_FADE_SUMMARY_TOP:-}" ] && NF_FLAGS+=(--scoring_fade_summary_top "$SCORING_FADE_SUMMARY_TOP")
[ -n "${SCORING_FADE_SUMMARY_BOTTOM:-}" ] && NF_FLAGS+=(--scoring_fade_summary_bottom "$SCORING_FADE_SUMMARY_BOTTOM")

# --- CLASS-SPECIFIC FLAGS ---
if [ "$CLASS" = "1" ]; then
    PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
    PRUNE_SECONDARY_LIST="${PRUNE_DIR}/${PRUNE_SECONDARY_FILE}"
    NF_FLAGS+=(
        --my_traits "$TRAIT_FILE"
        --traitname "$TRAIT"
        --secondary_trait "$SECONDARY_TRAIT"
        --n_trait "adult_necropsy_count"
        --c_trait "$C_TRAIT"
        --branch_trait "LQ"
        --prune_data
        --prune_list "$PRUNE_LIST"
        --prune_list_secondary "$PRUNE_SECONDARY_LIST"
        --reporting
        --contrast_selection
    )
elif [ "$CLASS" = "2" ]; then
    NF_FLAGS+=(
        --my_traits "$SIMPLE_TRAIT_FILE"
        --traitname "$TRAIT"
        --branch_trait "LQ"
        --discrete_method "$DISCRETE_METHOD"
        --n_trait ""
        --c_trait ""
        --secondary_trait ""
        --reporting
        --contrast_selection
    )
else
    echo "Error: CLASS must be 1 or 2, got: $CLASS"
    exit 1
fi

echo "=========================================="
echo " Running Nextflow Single Phenotype: $TRAIT"
echo " Profile  : local"
echo " Output   : $RESULTS_BASE"
echo "=========================================="

nextflow run "${REPO_DIR}/main.nf" "${NF_FLAGS[@]}" -w "$WORK_DIR"

if [ -f "${RESULTS_BASE}/low_contrasts.skip" ]; then
    echo "Warning: Skipped '${TRAIT}': too few foreground contrasts."
    exit 0
fi

echo "Success: Completed $TRAIT [CLASS $CLASS]"
