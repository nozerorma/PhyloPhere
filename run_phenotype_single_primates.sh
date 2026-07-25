#!/bin/bash
set -Eeuo pipefail

# --- ARGUMENTS ---
CLASS="${1:?Usage: $0 <CLASS> <TRAIT> [SECONDARY_TRAIT] [N_TRAIT] [C_TRAIT] [PRUNE_FILE] [PRUNE_SECONDARY_FILE] [DISCRETE_METHOD]}"
TRAIT="${2:?TRAIT must be provided}"
SECONDARY_TRAIT="${3:-}"
N_TRAIT="${4:-}"
C_TRAIT="${5:-}"
PRUNE_FILE="${6:-}"
PRUNE_SECONDARY_FILE="${7:-}"
DISCRETE_METHOD="${8:-}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp=$(date +%Y%m%d_%H%M%S)

# --- ENV ACTIVATION ---
source ~/.bashrc
conda deactivate
conda activate phylophere

# Setup unique Nextflow home directory to avoid file lock collisions on history database
export NXF_HOME="$WORK_BASE/.nextflow/.nextflow_${TRAIT}"
mkdir -p "$NXF_HOME"
ln -sfn /homes/users/mramon/.nextflow/plugins "$NXF_HOME/plugins"

# Create and cd into a unique run directory to prevent Nextflow session/cache database lock collisions
RUN_DIR="$WORK_BASE/${TRAIT}_nxf_run"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}"
WORK_DIR="${WORK_BASE}/${TRAIT}${TAG}"
mkdir -p "${RESULTS_BASE}" "${WORK_DIR}"

# --- NEXTFLOW FLAGS ---
NF_FLAGS=(
    -profile slurm
    -with-tower
    -name "${TRAIT}_$(cat /proc/sys/kernel/random/uuid 2>/dev/null | tr -d '-' | cut -c1-8 || printf '%08x' $RANDOM)"
    --alignment "$ALI_DIR"
    --tree "$TREE_FILE"
    --outdir "$RESULTS_BASE"
)

# General flags
if [ "$RESUME" = 1 ]; then
    NF_FLAGS+=( -resume )
fi

if [ "${TOY_MODE:-true}" = true ]; then
    NF_FLAGS+=(
        --toy_mode
        --toy_n "$TOY_N"
    )
fi

if [ "${RUN_REPORTING:-true}" = true ]; then
    NF_FLAGS+=(
        --reporting
    )
fi

if [ "${RUN_PRUNE:-true}" = true ]; then
    NF_FLAGS+=(
        --prune_data
    )
fi

if [ "${RUN_CAAS:-true}" = true ]; then
    NF_FLAGS+=(
        --contrast_selection
        --ct_tool "discovery,resample,bootstrap"
        --caas_config "$TRAIT_FILE"
        --perms_cycles "$PERMS_CYCLES"
        --caas_permulation_enrichment true
        --caas_full_perms "$CAAS_FULL_PERMS"
    )
fi

if [ "${RUN_DISAMBIGUATION:-true}" = true ]; then
    NF_FLAGS+=(
        --ct_disambiguation
        --ct_disambig_asr_mode "precomputed"
        --ct_disambig_asr_cache_dir "$ASR_CACHE_DIR"
        --ct_postproc
    )
fi

if [ "${RUN_ACCUMULATION:-true}" = true ]; then
    NF_FLAGS+=(
        --ct_accumulation
        --accumulation_n_randomizations "$N_RANDOMIZATIONS"
        --ct_disambig_asr_mode "precomputed"
        --ct_disambig_asr_cache_dir "$ASR_CACHE_DIR"
        --accumulation_entropy_dir "$ENTROPY_DIR"
    )
fi

# RERConverge Flags
if [ "${RUN_RER:-true}" = true ]; then
    NF_FLAGS+=(
        --rer_tool "${RER_TOOL}"
        --gene_trees "${GENE_TREES}"
        --rer_perm_batches "${RER_PERM_BATCHES}"
        --rer_perms_per_batch "${RER_PERMS_PER_BATCH}"
    )
fi

# FADE Flags
if [ "${RUN_FADE:-true}" = true ]; then
    NF_FLAGS+=(
        --fade
        --fade_mode "${FADE_MODE}"
    )
fi


# Scoring Flags
if [ "${RUN_SCORING:-true}" = true ]; then
    NF_FLAGS+=(
        --scoring
        --scoring_window_size_bp "${SCORING_WINDOW_SIZE_BP:-1000000}"
        --gene_ensembl_file "${GENE_ENSEMBL_FILE}"
    )
    if [ "${RUN_SCORING_STRESS:-true}" = true ]; then
        NF_FLAGS+=(--scoring_stress)
    fi
fi

# VEP Flags
if [ "${RUN_VEP:-true}" = true ]; then
    NF_FLAGS+=(
        --vep
        --vep_primateai_db "${VEP_PRIMATEAI_DB}"
        --vep_map_dir "${MAP_DIR}"
        --scoring_vep_cosmic "${SCORING_VEP_COSMIC:-}"
    )
fi

# Enrichment Flags
if [ "${RUN_ENRICHMENT:-true}" = true ]; then
    NF_FLAGS+=(
        --enrichment
        --gmt_dir "${GMT_DIR}"
        --string_db_dir "${STRING_DB_DIR:-}"
    )
fi

# POSENRICH (position-wise enrichment) flags
if [ "${RUN_POSENRICH:-true}" = true ]; then
    NF_FLAGS+=(
        --posenrich
        --egg_members_file "${EGG_MEMBERS_FILE}"
        --egg_annotations_file "${EGG_ANNOTATIONS_FILE}"
        --cosmic_db "${COSMIC_DB}"
        --domain_variability_file "${DOMAIN_VARIABILITY_FILE}"
        --ucr_positions_file "${UCR_POSITIONS_FILE}"
        --fubar_sites_file "${FUBAR_SITES_FILE}"
    )
fi

# Species and Tax ID Flags
[ -n "${ALI_SP_NAMES:-}" ] && NF_FLAGS+=(--ali_sp_names "$ALI_SP_NAMES")
[ -n "${INPUT_TAX_ID:-}" ] && NF_FLAGS+=(--tax_id "$INPUT_TAX_ID")

# Precomputed FADE/RER summaries for scoring
[ -n "${SCORING_RER_INPUT:-}" ] && NF_FLAGS+=(--scoring_rer_input "$SCORING_RER_INPUT")
[ -n "${SCORING_RER_PERMS_INPUT:-}" ] && NF_FLAGS+=(--rer_perms_file "$SCORING_RER_PERMS_INPUT")
[ -n "${RER_UNIVERSE_FILE:-}" ] && NF_FLAGS+=(--rer_universe_file "$RER_UNIVERSE_FILE")
[ -n "${SCORING_FADE_SUMMARY_TOP:-}" ] && NF_FLAGS+=(--scoring_fade_summary_top "$SCORING_FADE_SUMMARY_TOP")
[ -n "${SCORING_FADE_SUMMARY_BOTTOM:-}" ] && NF_FLAGS+=(--scoring_fade_summary_bottom "$SCORING_FADE_SUMMARY_BOTTOM")
[ -n "${FADE_UNIVERSE_FILE:-}" ] && NF_FLAGS+=(--fade_universe_file "$FADE_UNIVERSE_FILE")

# --- CLASS-SPECIFIC FLAGS ---
if [ "$CLASS" = "1" ]; then
    PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
    PRUNE_SECONDARY_LIST="${PRUNE_DIR}/${PRUNE_SECONDARY_FILE}"
    NF_FLAGS+=(
        --my_traits "$TRAIT_FILE"
        --traitname "$TRAIT"
        --secondary_trait "$SECONDARY_TRAIT"
        --n_trait "$N_TRAIT"
        --c_trait "$C_TRAIT"
        --branch_trait "LQ"
        --prune_list "$PRUNE_LIST"
        --prune_list_secondary "$PRUNE_SECONDARY_LIST"
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
    )
else
    echo "Error: CLASS must be 1 or 2, got: $CLASS"
    exit 1
fi

echo "======================================================"
echo "Running Nextflow Single Phenotype: $TRAIT"
echo "SLURM TASK ID     : $SLURM_ARRAY_TASK_ID"

echo "Profile           : slurm"
echo "Trait             : $TRAIT  [CLASS $CLASS]"
echo "Alignments        : $ALI_DIR"
echo "Output            : $RESULTS_BASE"
echo "Work dir          : $WORK_DIR"

echo "Resuming?         : $([ "$RESUME" = 1 ] && echo yes || echo 'no (clean run)')"
echo "======================================================"

nextflow run "${REPO_DIR}/main.nf" "${NF_FLAGS[@]}" -w "$WORK_DIR"

if [ -f "${RESULTS_BASE}/low_contrasts.skip" ]; then
    echo "Warning: Skipped '${TRAIT}': too few foreground contrasts."
    exit 0
fi

echo "Success: Completed $TRAIT [CLASS $CLASS]"
