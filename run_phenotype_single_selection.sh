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
# PHYLOPHERE: Single-Phenotype Selection-Only Runner
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Runs only the selection tools (FADE, MoleRate, and optionally RERConverge)
# for a single phenotype, re-using gene sets produced by a prior completed
# full run.  All CT, signification, disambiguation, postproc, ORA, string, and
# accumulation steps are skipped.
#
# Gene sets are read directly from SOURCE_RUN_SUBDIR (set via env var or inside
# this script).  The resolved path is:
#   CAAS_OUTBASE/<TRAIT>[_toy]/<SOURCE_RUN_SUBDIR>/filter/
#
# Called by SBATCH_run_phenotypes_selection.sh (one SLURM array task per phenotype),
# or run directly:
#
# Usage (direct):
#   run_phenotype_single_selection.sh <CLASS> <TRAIT> [SECONDARY_TRAIT] [C_TRAIT] \
#                                     [PRUNE_FILE] [PRUNE_SECONDARY_FILE] [DISCRETE_METHOD]
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
# File: run_phenotype_single_selection.sh
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
CAAS_OUTBASE="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS"
WORK_BASE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs"

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
else
    TAG=""
    ALI_DIR="${DATADIR}/2.Alignments/Primate_alignments"
fi

# ── Selection tool toggles ────────────────────────────────────────────────────
RUN_FADE="${RUN_FADE:-true}"
RUN_MOLERATE="${RUN_MOLERATE:-true}"
FADE_MODE="${FADE_MODE:-gene_set}"
MOLERATE_MODE="${MOLERATE_MODE:-gene_set}"

# ── RERConverge toggles ───────────────────────────────────────────────────────
# Kept opt-in: RERConverge requires gene trees whose species names match the
# trait file.  Enable once species-name consistency has been verified.
RUN_RER="${RUN_RER:-false}"
RER_TOOL="${RER_TOOL:-build_trait,build_tree,build_matrix,continuous}"
RER_GENE_SET_MODE="${RER_GENE_SET_MODE:-gene_set}"
GENE_TREES="${GENE_TREES:-${DATADIR}/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt}"

# ============================================================
# ENVIRONMENT SETUP
# ============================================================
module load Nextflow
module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate phylophere

timestamp=$(date +%Y%m%d_%H%M%S)

# ── Source-run derived paths ───────────────────────────────────────────────────
# Gene sets and accumulation outputs are read from the prior completed run that
# lives under CAAS_OUTBASE/<TRAIT>[_toy]/<SOURCE_RUN_SUBDIR>/filter/.
# SOURCE_RUN_SUBDIR must be set (exported by SBATCH_run_phenotypes_selection.sh
# or manually before calling this script).
SOURCE_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${SOURCE_RUN_SUBDIR}/filter"

# Gene-set inputs for FADE / MoleRate / RER (precomputed by a prior full run)
SOURCE_POSTPROC_TOP="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
SOURCE_POSTPROC_BOTTOM="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"
SOURCE_ACCUMULATION_CSV="${SOURCE_BASE}/accumulation/aggregation/accumulation_global.csv"

# Short random hex ID used to name the Nextflow run (8 chars)
NXF_RUN_ID="$(cat /proc/sys/kernel/random/uuid 2>/dev/null | tr -d '-' | cut -c1-8 \
              || printf '%08x' $RANDOM)"

# ============================================================
# INPUT VALIDATION
# ============================================================
for f in "$TREE_FILE" "$SOURCE_POSTPROC_TOP" "$SOURCE_POSTPROC_BOTTOM" "$SOURCE_ACCUMULATION_CSV"; do
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
# Gene sets come directly from the source run — no CT needed.
# ============================================================
FADE_NF_FLAGS=()
if [ "$RUN_FADE" = true ]; then
    FADE_NF_FLAGS=(
        --fade
        --fade_mode             "$FADE_MODE"
        --fade_postproc_top     "$SOURCE_POSTPROC_TOP"
        --fade_postproc_bottom  "$SOURCE_POSTPROC_BOTTOM"
        --fade_accumulation_top    "$SOURCE_ACCUMULATION_CSV"
        --fade_accumulation_bottom "$SOURCE_ACCUMULATION_CSV"
    )
fi

MOLERATE_NF_FLAGS=()
if [ "$RUN_MOLERATE" = true ]; then
    MOLERATE_NF_FLAGS=(
        --molerate
        --molerate_mode                "$MOLERATE_MODE"
        --molerate_postproc_top        "$SOURCE_POSTPROC_TOP"
        --molerate_postproc_bottom     "$SOURCE_POSTPROC_BOTTOM"
        --molerate_accumulation_top    "$SOURCE_ACCUMULATION_CSV"
        --molerate_accumulation_bottom "$SOURCE_ACCUMULATION_CSV"
    )
fi

# ── RERConverge flags ─────────────────────────────────────────────────────────
RER_NF_FLAGS=()
if [ "$RUN_RER" = true ]; then
    RER_NF_FLAGS=(
        --rer_tool                "$RER_TOOL"
        --rer_gene_set_mode       "$RER_GENE_SET_MODE"
        --gene_trees              "$GENE_TREES"
        --rer_postproc_top        "$SOURCE_POSTPROC_TOP"
        --rer_postproc_bottom     "$SOURCE_POSTPROC_BOTTOM"
        --rer_accumulation_top    "$SOURCE_ACCUMULATION_CSV"
        --rer_accumulation_bottom "$SOURCE_ACCUMULATION_CSV"
    )
fi

# ============================================================
# COMMON NEXTFLOW FLAGS
# CT is not run: --ct_tool is intentionally absent to avoid triggering it.
# ============================================================
COMMON_NF_FLAGS=(
    -with-tower
    -name "${TRAIT}_sel_${NXF_RUN_ID}"
    -profile local
    --alignment "$ALI_DIR"
    --tree      "$TREE_FILE"
    --reporting false
    --contrast_selection false
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
        -w "$workdir" \
        --outdir "$outdir" \
        "${extra[@]}"
}

# ============================================================
# HEADER
# ============================================================
echo "=========================================="
echo " PHYLOPHERE — SELECTION-ONLY RUNNER"
echo " CLASS        : $CLASS"
echo " TRAIT        : $TRAIT"
echo " Timestamp    : $timestamp"
echo " NF run name  : ${TRAIT}_sel_${NXF_RUN_ID}"
echo " IS_TOY       : $IS_TOY  (tag: '${TAG}')"
echo " SOURCE_SUBDIR: $SOURCE_RUN_SUBDIR"
echo " RUN_FADE     : $RUN_FADE  (mode=${FADE_MODE})"
echo " RUN_MOLERATE : $RUN_MOLERATE  (mode=${MOLERATE_MODE})"
echo " RUN_RER      : $RUN_RER  (tool=${RER_TOOL})"
echo " Gene set TOP : $SOURCE_POSTPROC_TOP"
echo " Gene set BOT : $SOURCE_POSTPROC_BOTTOM"
echo " Accum. CSV   : $SOURCE_ACCUMULATION_CSV"
echo "=========================================="

RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}_selection/${timestamp}"
WORK_DIR="${WORK_BASE}/${TRAIT}${TAG}_selection/${timestamp}"

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
    echo " Running CLASS 1 (selection only): $TRAIT"
    echo "   Secondary   : $SECONDARY_TRAIT"
    echo "   n_trait     : $N_TRAIT_PRUNED"
    echo "   c_trait     : $C_TRAIT"
    echo "   Branch trait: $BRANCH_TRAIT"
    echo "   Prune list  : $PRUNE_LIST"
    echo "   Output      : ${RESULTS_BASE}"
    echo "------------------------------------------"

    TRAIT_FLAGS=(
        --my_traits            "$TRAIT_FILE"
        --traitname            "$TRAIT"
        --secondary_trait      "$SECONDARY_TRAIT"
        --n_trait              "$N_TRAIT_PRUNED"
        --c_trait              "$C_TRAIT"
        --branch_trait         "$BRANCH_TRAIT"
        --prune_data false
        --prune_list           "$PRUNE_LIST"
        --prune_list_secondary "$PRUNE_SECONDARY_LIST"
        "${FADE_NF_FLAGS[@]:-}"
        "${MOLERATE_NF_FLAGS[@]:-}"
        "${RER_NF_FLAGS[@]:-}"
    )

    run_pipeline \
        "${RESULTS_BASE}" \
        "${WORK_DIR}" \
        "${TRAIT_FLAGS[@]}"

# ============================================================
# CLASS 2 — Diet / Simple (María Sánchez Bermúdez)
# ============================================================
elif [ "$CLASS" = "2" ]; then

    echo ""
    echo "------------------------------------------"
    echo " Running CLASS 2 (selection only): $TRAIT"
    echo "   Trait file    : $SIMPLE_TRAIT_FILE"
    echo "   Discrete meth : $DISCRETE_METHOD"
    echo "   Output        : ${RESULTS_BASE}"
    echo "------------------------------------------"

    TRAIT_FLAGS=(
        --my_traits       "$SIMPLE_TRAIT_FILE"
        --traitname       "$TRAIT"
        --branch_trait    "$BRANCH_TRAIT"
        --discrete_method "$DISCRETE_METHOD"
        --n_trait         ""
        --c_trait         ""
        --secondary_trait ""
        "${FADE_NF_FLAGS[@]:-}"
        "${MOLERATE_NF_FLAGS[@]:-}"
        "${RER_NF_FLAGS[@]:-}"
    )

    run_pipeline \
        "${RESULTS_BASE}" \
        "${WORK_DIR}" \
        "${TRAIT_FLAGS[@]}"

else
    echo "Error: CLASS must be 1 or 2, got: $CLASS"
    exit 1
fi

echo ""
echo " ✓ Completed: $TRAIT  [CLASS ${CLASS}] — selection only"
echo " Results    : ${RESULTS_BASE}"
echo "=========================================="
