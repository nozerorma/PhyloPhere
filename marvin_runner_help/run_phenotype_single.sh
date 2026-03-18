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
# PHYLOPHERE: Single-Phenotype Runner
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Runs one complete PhyloPhere pipeline for a single phenotype.
# Called by SBATCH_run_phenotypes.sh (one SLURM array task per phenotype).
#
# Usage (direct):
#   run_phenotype_single.sh <CLASS> <TRAIT> [SECONDARY_TRAIT] [C_TRAIT] \
#                           [PRUNE_FILE] [PRUNE_SECONDARY_FILE] [DISCRETE_METHOD]
#
#   CLASS               : 1 (cancer / pruned-secondary) | 2 (simple / diet)
#   TRAIT               : trait column name in the trait CSV
#   SECONDARY_TRAIT     : (CLASS 1) secondary trait column name; empty for CLASS 2
#   C_TRAIT             : (CLASS 1) case-count column; empty for CLASS 2
#   PRUNE_FILE          : (CLASS 1) basename of primary prune list in PRUNE_DIR
#   PRUNE_SECONDARY_FILE: (CLASS 1) basename of secondary prune list in PRUNE_DIR
#   DISCRETE_METHOD     : (CLASS 2) discretisation method; default "decile"
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: run_phenotype_single.sh
#

set -Eeuo pipefail

# ─── ARGUMENTS ───────────────────────────────────────────────────────────────
CLASS="${1:?Usage: $0 <CLASS> <TRAIT> [SECONDARY_TRAIT] [C_TRAIT] [PRUNE_FILE] [PRUNE_SECONDARY_FILE] [DISCRETE_METHOD]}"
TRAIT="${2:?TRAIT must be provided as \$2}"
SECONDARY_TRAIT="${3:-}"
C_TRAIT="${4:-}"
PRUNE_FILE="${5:-}"
PRUNE_SECONDARY_FILE="${6:-}"
DISCRETE_METHOD="${7:-decile}"

# ─── REPO ROOT ───────────────────────────────────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# CLUSTER / ENVIRONMENT CONFIGURATION
# ─────────────────────────────────────────────────────────────
# All cluster-specific paths are controlled here.
# Adjust before submitting on a new cluster.
# ============================================================

DATADIR="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data"
CAAS_OUTBASE="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet"
WORK_BASE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs/diet"
ASR_CACHE_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/asr"

# ── Source run sub-directory ──────────────────────────────────────────────────
# Under CAAS_OUTBASE/<TRAIT>/  there must be a sub-directory containing a
# completed prior run's filter/ outputs (caastools/, postproc/, accumulation/).
# Set SOURCE_RUN_SUBDIR to "runtime" (default) or the timestamped dir name.
# The full resolved path becomes: CAAS_OUTBASE/<TRAIT>/<SOURCE_RUN_SUBDIR>/filter
SOURCE_RUN_SUBDIR="${SOURCE_RUN_SUBDIR:-runtime}"

TRAIT_FILE="${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"
TREE_FILE="${DATADIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
PRUNE_DIR="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP"
SIMPLE_TRAIT_FILE="${DATADIR}/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv"

N_TRAIT_PRUNED="adult_necropsy_count"
BRANCH_TRAIT="LQ"

# ── Toy / Full mode ──────────────────────────────────────────────────────────
IS_TOY="${IS_TOY:-false}"

if [ "$IS_TOY" = true ]; then
    TAG="_toy"
    ALI_DIR="${DATADIR}/2.Alignments/Ali_toy"
    CYCLES="100"
    N_RANDOMIZATIONS="1000"
else
    TAG=""
    ALI_DIR="${DATADIR}/2.Alignments/Primate_alignments"
    CYCLES="1000000"
    N_RANDOMIZATIONS="1000000"
fi

# ── Selection analysis toggles ───────────────────────────────────────────────
RUN_FADE="${RUN_FADE:-false}"
RUN_MOLERATE="${RUN_MOLERATE:-false}"
FADE_MODE="${FADE_MODE:-gene_set}"
MOLERATE_MODE="${MOLERATE_MODE:-gene_set}"

# ── RERConverge toggles ──────────────────────────────────────────────────────
RUN_RER="${RUN_RER:-false}"
RER_TOOL="${RER_TOOL:-build_trait,build_tree,build_matrix,continuous}"
RER_GENE_SET_MODE="${RER_GENE_SET_MODE:-gene_set}"
GENE_TREES="${GENE_TREES:-${DATADIR}/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt}"

# ============================================================
# ENVIRONMENT SETUP
# ============================================================
# module purge
# module load modulepath/haswell
# module load Nextflow
# module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate phylophere

timestamp=$(date +%Y%m%d_%H%M%S)

# ── Source-run derived paths (resolved here so TRAIT / TAG are known) ─────────
# All CT outputs (discovery, resample, bootstrap, background) and gene-set
# inputs (postproc top/bottom, accumulation CSV) are auto-derived from the
# prior completed run that lives under SOURCE_RUN_SUBDIR/filter/.
SOURCE_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${SOURCE_RUN_SUBDIR}"

# CT raw outputs
SOURCE_DISCOVERY="${SOURCE_BASE}/caastools/discovery.tab"
SOURCE_RESAMPLE_DIR="${SOURCE_BASE}/resample/nw_tree.resampled.output"
SOURCE_BOOTSTRAP="${SOURCE_BASE}/caastools/bootstrap.tab"
SOURCE_BACKGROUND="${SOURCE_BASE}/caastools/background_genes.output"
SOURCE_CONFIG="${SOURCE_BASE}/data_exploration/2.CT/1.Traitfiles/traitfile.tab"

# # Gene-set inputs for FADE / MoleRate / RER (used in gene_set mode)
# SOURCE_POSTPROC_TOP="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
# SOURCE_POSTPROC_BOTTOM="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"
# SOURCE_ACCUMULATION_CSV="${SOURCE_BASE}/accumulation/aggregation/accumulation_global.csv"

# Short random hex ID used to name the Nextflow run (8 chars)
NXF_RUN_ID="$(cat /proc/sys/kernel/random/uuid 2>/dev/null | tr -d '-' | cut -c1-8 \
              || printf '%08x' $RANDOM)"

# ============================================================
# INPUT VALIDATION
# ============================================================
for f in "$TREE_FILE"; do
    if [ ! -f "$f" ]; then
        echo "Error: required file not found: $f"
        exit 1
    fi
done

if [ "$CLASS" = "1" ] && [ ! -f "$TRAIT_FILE" ]; then
    echo "Error: TRAIT_FILE not found: $TRAIT_FILE"
    exit 1
fi

if [ "$CLASS" = "2" ] && [ ! -f "$SIMPLE_TRAIT_FILE" ]; then
    echo "Error: SIMPLE_TRAIT_FILE not found: $SIMPLE_TRAIT_FILE"
    exit 1
fi

if [ ! -d "$ALI_DIR" ]; then
    echo "Error: alignment directory not found: $ALI_DIR"
    exit 1
fi

# ============================================================
# FADE / MOLERATE / RER FLAGS
# Gene-set paths are auto-derived from the source run.
# Override SOURCE_POSTPROC_* / SOURCE_ACCUMULATION_CSV above if needed.
# ============================================================
FADE_NF_FLAGS=()
if [ "$RUN_FADE" = true ]; then
    FADE_NF_FLAGS=(
        --fade
        --fade_mode "$FADE_MODE"
    )
fi

MOLERATE_NF_FLAGS=()
if [ "$RUN_MOLERATE" = true ]; then
    MOLERATE_NF_FLAGS=(
        --molerate
        --molerate_mode "$MOLERATE_MODE"
    )
fi

# ── RERConverge flags ────────────────────────────────────────────────────────
RER_NF_FLAGS=()
if [ "$RUN_RER" = true ]; then
    RER_NF_FLAGS=(
        --rer_tool          "$RER_TOOL"
        --rer_gene_set_mode "$RER_GENE_SET_MODE"
        --gene_trees        "$GENE_TREES"
    )
fi

# ============================================================
# COMMON NEXTFLOW FLAGS
# ============================================================
COMMON_NF_FLAGS=(
    -with-tower
    -profile slurm
    -name "${TRAIT}_${NXF_RUN_ID}"
    --ct_tool "discovery,resample,bootstrap" #resample,bootstrap results are passed as inputs, so we only run the discovery step
    --alignment  "$ALI_DIR"
    --tree        "$TREE_FILE"
    --cycles      "$CYCLES"
    --accumulation_n_randomizations "$N_RANDOMIZATIONS"
    --ct_disambig_asr_mode      "precomputed"
    --ct_disambig_asr_cache_dir "${ASR_CACHE_DIR}"
    --discovery_out "${SOURCE_DISCOVERY}"
    --resample_out "${SOURCE_RESAMPLE_DIR}"
    --bootstrap_input "${SOURCE_BOOTSTRAP}"
    --background_input "${SOURCE_BACKGROUND}"
)

# ============================================================
# HELPER: run pipeline
# ============================================================
run_pipeline() {
    local outdir="$1"
    local workdir="$2"
    shift 2
    local extra=("$@")

    mkdir -p "$outdir" "$workdir"

    nextflow run "${REPO_DIR}/main.nf" \
        "${COMMON_NF_FLAGS[@]}" \
        --caas_postproc_mode "filter" \
        -w "$workdir" \
        --outdir "$outdir" \
        "${extra[@]}"
}

# ============================================================
# HEADER
# ============================================================
echo "=========================================="
echo " PHYLOPHERE SINGLE-PHENOTYPE RUNNER"
echo " CLASS      : $CLASS"
echo " TRAIT      : $TRAIT"
echo " Timestamp  : $timestamp"
echo " NF run name: ${TRAIT}_${NXF_RUN_ID}"
echo " IS_TOY     : $IS_TOY  (tag: '${TAG}')"
echo " Cycles     : $CYCLES"
echo " N_Rand     : $N_RANDOMIZATIONS"
echo " RUN_FADE   : $RUN_FADE  (mode=${FADE_MODE})"
echo " RUN_MOLERATE: $RUN_MOLERATE  (mode=${MOLERATE_MODE})"
echo " RUN_RER    : $RUN_RER  (tool=${RER_TOOL})"
echo "=========================================="

RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${timestamp}"
WORK_DIR="${WORK_BASE}/${TRAIT}${TAG}/${timestamp}"

# ============================================================
# CLASS 1 — Cancer / Pruned-Secondary
# ============================================================
if [ "$CLASS" = "1" ]; then

    PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
    PRUNE_SECONDARY_LIST="${PRUNE_DIR}/${PRUNE_SECONDARY_FILE}"

    for pf in "$PRUNE_LIST" "$PRUNE_SECONDARY_LIST"; do
        if [ ! -f "$pf" ]; then
            echo "Error: prune list not found: $pf"
            exit 1
        fi
    done

    echo ""
    echo "------------------------------------------"
    echo " Running CLASS 1: $TRAIT"
    echo "   Secondary     : $SECONDARY_TRAIT"
    echo "   n_trait       : $N_TRAIT_PRUNED"
    echo "   c_trait       : $C_TRAIT"
    echo "   Branch trait  : $BRANCH_TRAIT"
    echo "   Prune list    : $PRUNE_LIST"
    echo "   Prune (sec.)  : $PRUNE_SECONDARY_LIST"
    echo "   Output        : ${RESULTS_BASE}/filter"
    echo "------------------------------------------"

    TRAIT_FLAGS=(
        --my_traits            "$TRAIT_FILE"
        --traitname            "$TRAIT"
        --secondary_trait      "$SECONDARY_TRAIT"
        --n_trait              "$N_TRAIT_PRUNED"
        --c_trait              "$C_TRAIT"
        --branch_trait         "$BRANCH_TRAIT"
        --caas_config          "$SOURCE_CONFIG"
        --prune_data
        --prune_list           "$PRUNE_LIST"
        --prune_list_secondary "$PRUNE_SECONDARY_LIST"
        --reporting false
        --contrast_selection false
        --ct_postproc
        --ct_signification
        --ct_disambiguation
        --ora
        --string
        --ct_accumulation
        "${FADE_NF_FLAGS[@]:-}"
        "${MOLERATE_NF_FLAGS[@]:-}"
        "${RER_NF_FLAGS[@]:-}"
    )

    run_pipeline \
        "${RESULTS_BASE}/filter" \
        "${WORK_DIR}/filter" \
        "${TRAIT_FLAGS[@]}"

# ============================================================
# CLASS 2 — Diet / Simple (María Sánchez Bermúdez)
# ============================================================
elif [ "$CLASS" = "2" ]; then

    echo ""
    echo "------------------------------------------"
    echo " Running CLASS 2: $TRAIT"
    echo "   Trait file    : $SIMPLE_TRAIT_FILE"
    echo "   Discrete meth : $DISCRETE_METHOD"
    echo "   (no pruning, no secondary/n/c traits)"
    echo "   Output        : ${RESULTS_BASE}/filter"
    echo "------------------------------------------"

    TRAIT_FLAGS=(
        --my_traits       "$SIMPLE_TRAIT_FILE"
        --traitname       "$TRAIT"
        --branch_trait    "$BRANCH_TRAIT"
        --discrete_method "$DISCRETE_METHOD"
        --n_trait         ""
        --c_trait         ""
        --secondary_trait ""
        --reporting
        --contrast_selection
        --ct_postproc
        --ct_signification
        --ct_disambiguation
        --ora
        --string
        --ct_accumulation
        "${FADE_NF_FLAGS[@]:-}"
        "${MOLERATE_NF_FLAGS[@]:-}"
        "${RER_NF_FLAGS[@]:-}"
    )

    run_pipeline \
        "${RESULTS_BASE}/filter" \
        "${WORK_DIR}/filter" \
        "${TRAIT_FLAGS[@]}"

else
    echo "Error: CLASS must be 1 or 2, got: $CLASS"
    exit 1
fi

# ── Low-contrast gate ────────────────────────────────────────────────────────
if [ -f "${RESULTS_BASE}/filter/low_contrasts.skip" ]; then
    echo ""
    echo " ⚠ Skipped '${TRAIT}': too few foreground contrasts."
    echo "   Detail: $(cat "${RESULTS_BASE}/filter/low_contrasts.skip")"
    echo ""
    exit 0
fi

echo ""
echo " ✓ Completed: $TRAIT  [CLASS ${CLASS}]"
echo " Results    : ${RESULTS_BASE}/filter"
echo "=========================================="
