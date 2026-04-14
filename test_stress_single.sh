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
# PHYLOPHERE: Single-Trait Stress-Test Runner
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Trait-agnostic stress-test runner. Exercises every pipeline module in both
# integrated (end-to-end) and standalone (per-module, file-input) modes for a
# single phenotype. Called by SBATCH_test_stress.sh (one SLURM array task per
# phenotype) or directly from the command line.
#
# All trait-specific parameters are read from environment variables so the
# SBATCH controller can inject them without touching this file.
#
# ─── USAGE ────────────────────────────────────────────────────────────────────
# Direct (CLASS 1 example):
#   export TRAIT="neoplasia_prevalence"
#   export SECONDARY_TRAIT="malignant_prevalence"
#   export N_TRAIT="adult_necropsy_count"
#   export C_TRAIT="neoplasia_necropsy"
#   export PRUNE_FILE="neoplasia_exclude.txt"
#   export PRUNE_FILE_SECONDARY="malignant_exclude.txt"
#   export TRAIT_CLASS=1
#   bash test_stress_single.sh
#
# Direct (CLASS 2 example):
#   export TRAIT="frug_idx"
#   export TRAIT_CLASS=2
#   bash test_stress_single.sh
#
# Toggle variables (set before calling or export from SBATCH controller):
#   See MASTER TOGGLES section below — each has a safe default.
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: test_stress_single.sh
#

# set -Eeuo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# ENVIRONMENT ACTIVATION
# ─────────────────────────────────────────────────────────────────────────────
source ~/.bashrc
conda deactivate
conda activate phylophere

# ─────────────────────────────────────────────────────────────────────────────
# TIMESTAMP + REPO ROOT
# ─────────────────────────────────────────────────────────────────────────────
timestamp=$(date +%Y%m%d_%H%M%S)
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   TRAIT CONFIGURATION  —  injected from SBATCH controller or set by caller
# ══════════════════════════════════════════════════════════════════════════════

# CLASS 1 (cancer / pruned-secondary) or CLASS 2 (diet / simple)
TRAIT_CLASS="${TRAIT_CLASS:-1}"

# Primary trait column name in the trait CSV
TRAIT="${TRAIT:?TRAIT env var must be set}"

# CLASS 1 optional columns (leave empty for CLASS 2)
SECONDARY_TRAIT="${SECONDARY_TRAIT:-}"
N_TRAIT="${N_TRAIT:-adult_necropsy_count}"
C_TRAIT="${C_TRAIT:-}"
BRANCH_TRAIT="${BRANCH_TRAIT:-LQ}"

# Prune-list basenames (CLASS 1 only; resolved against PRUNE_DIR below)
PRUNE_FILE="${PRUNE_FILE:-}"
PRUNE_FILE_SECONDARY="${PRUNE_FILE_SECONDARY:-}"

# Sub-directory of a completed prior run whose filter/ outputs are used as
# standalone-module inputs.  "runtime" = the stable symlink / current run.
SOURCE_RUN_SUBDIR="${SOURCE_RUN_SUBDIR:-runtime}"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   MASTER TOGGLES  —  edit this section or export from the SBATCH controller
# ══════════════════════════════════════════════════════════════════════════════

# ── Infrastructure ──────────────────────────────────────────────────────────
CLEAN_WORK="${CLEAN_WORK:-false}"

# ── Integrated run toggles ──────────────────────────────────────────────────
RUN_INT_FILTER="${RUN_INT_FILTER:-false}"
RUN_INT_EXPLORATORY="${RUN_INT_EXPLORATORY:-false}"
INT_RESUME="${INT_RESUME:-false}"

INT_PRUNE_DATA="${INT_PRUNE_DATA:-false}"
INT_REPORTING="${INT_REPORTING:-false}"
INT_CONTRAST_SELECTION="${INT_CONTRAST_SELECTION:-false}"
INT_CT_SIGNIFICATION="${INT_CT_SIGNIFICATION:-false}"
INT_CT_DISAMBIGUATION="${INT_CT_DISAMBIGUATION:-false}"
INT_ORA="${INT_ORA:-true}"
INT_STRING="${INT_STRING:-true}"
INT_CT_ACCUMULATION="${INT_CT_ACCUMULATION:-true}"
INT_VEP="${INT_VEP:-true}"
INT_FADE="${INT_FADE:-false}"
INT_MOLERATE="${INT_MOLERATE:-false}"
INT_RER="${INT_RER:-true}"
INT_SCORING="${INT_SCORING:-true}"
INT_SCORING_STRESS="${INT_SCORING_STRESS:-true}"
INT_SCORING_STRESS_TOP_N="${INT_SCORING_STRESS_TOP_N:-25}"

INT_USE_SECONDARY_TRAIT="${INT_USE_SECONDARY_TRAIT:-true}"
INT_USE_BRANCH_TRAIT="${INT_USE_BRANCH_TRAIT:-true}"
INT_USE_N_TRAIT="${INT_USE_N_TRAIT:-true}"
INT_USE_C_TRAIT="${INT_USE_C_TRAIT:-true}"

INT_DISAMBIG_ASR_MODE="${INT_DISAMBIG_ASR_MODE:-precomputed}"
INT_DISAMBIG_ASR_CACHE_DIR="${INT_DISAMBIG_ASR_CACHE_DIR:-/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/asr}"

# ── Standalone run toggles ───────────────────────────────────────────────────
RUN_SA_CT="${RUN_SA_CT:-false}"
RUN_SA_SIGNIFICATION="${RUN_SA_SIGNIFICATION:-false}"
RUN_SA_DISAMBIGUATION="${RUN_SA_DISAMBIGUATION:-false}"
SA_DISAMBIG_COMPUTE="${SA_DISAMBIG_COMPUTE:-false}"
SA_DISAMBIG_PRECOMPUTED="${SA_DISAMBIG_PRECOMPUTED:-true}"
RUN_SA_POSTPROC_FILTER="${RUN_SA_POSTPROC_FILTER:-true}"
RUN_SA_POSTPROC_EXPLORATORY="${RUN_SA_POSTPROC_EXPLORATORY:-false}"
RUN_SA_ORA="${RUN_SA_ORA:-false}"
RUN_SA_STRING="${RUN_SA_STRING:-false}"
RUN_SA_ACCUMULATION="${RUN_SA_ACCUMULATION:-false}"
RUN_SA_REPORTING="${RUN_SA_REPORTING:-false}"
RUN_SA_CONTRAST_SELECTION="${RUN_SA_CONTRAST_SELECTION:-false}"
RUN_SA_VEP="${RUN_SA_VEP:-false}"
RUN_SA_SCORING="${RUN_SA_SCORING:-false}"
SA_SCORING_STRESS="${SA_SCORING_STRESS:-true}"
SA_SCORING_STRESS_TOP_N="${SA_SCORING_STRESS_TOP_N:-25}"
RUN_SA_FADE="${RUN_SA_FADE:-false}"
RUN_SA_MOLERATE="${RUN_SA_MOLERATE:-false}"
RUN_SA_RER="${RUN_SA_RER:-false}"

SA_SELECTION_MODE_ALL="${SA_SELECTION_MODE_ALL:-true}"
SA_SELECTION_MODE_GENESET_PP="${SA_SELECTION_MODE_GENESET_PP:-false}"

SA_RER_TOOL="${SA_RER_TOOL:-build_trait,build_tree,build_matrix,continuous}"
SA_RER_CONTINUOUS_ONLY="${SA_RER_CONTINUOUS_ONLY:-false}"
SA_RER_PERM_BATCHES="${SA_RER_PERM_BATCHES:-100}"
SA_RER_PERMS_PER_BATCH="${SA_RER_PERMS_PER_BATCH:-100}"
SA_RER_PERM_MODE="${SA_RER_PERM_MODE:-cc}"
SA_RER_GMT_FILE="${SA_RER_GMT_FILE:-/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere/subworkflows/RERCONVERGE/dat/c2.cp.pid.v2026.1.Hs.symbols.gmt}"

# ── Toy mode ─────────────────────────────────────────────────────────────────
TOY_MODE="${TOY_MODE:-false}"
TOY_N="${TOY_N:-1000}"

# ══════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# CLUSTER / ENVIRONMENT PATHS
# ─────────────────────────────────────────────────────────────────────────────
DATADIR="${DATADIR:-/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data}"
CAAS_OUTBASE="${CAAS_OUTBASE:-/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/final}"
WORK_BASE="${WORK_BASE:-/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs/final}"
ALI_ARCHIVE_ROOT="${ALI_ARCHIVE_ROOT:-/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons/alignments/Primates_BMGE}"

ALI_DIR="${ALI_DIR:-${ALI_ARCHIVE_ROOT}/PROT}"
INPUT_TRAITS="${INPUT_TRAITS:-${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv}"
SIMPLE_TRAIT_FILE="${SIMPLE_TRAIT_FILE:-${DATADIR}/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv}"
INPUT_TREE="${INPUT_TREE:-${DATADIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree}"
PRUNE_DIR="${PRUNE_DIR:-${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP}"
INPUT_GENE_TREES="${INPUT_GENE_TREES:-${DATADIR}/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt}"
INPUT_ASR_CACHE="${INPUT_ASR_CACHE:-${DATADIR}/asr_cache}"
INPUT_ASR_CACHED="${INPUT_ASR_CACHED:-${DATADIR}/asr}"

INPUT_VEP_CDS_DIR="${INPUT_VEP_CDS_DIR:-${ALI_ARCHIVE_ROOT}/TRIM}"
INPUT_VEP_TRACK_DIR="${INPUT_VEP_TRACK_DIR:-${ALI_ARCHIVE_ROOT}/TRACK}"
INPUT_VEP_PRIMATEAI_DB="${INPUT_VEP_PRIMATEAI_DB:-}"

CYCLES="${CYCLES:-1000000}"

# ── Resolve per-trait paths (depend on TRAIT being set above) ─────────────────
# Select the right trait CSV depending on class
if [ "$TRAIT_CLASS" = "2" ]; then
    INPUT_TRAITS="$SIMPLE_TRAIT_FILE"
fi

PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
PRUNE_LIST_SECONDARY="${PRUNE_DIR}/${PRUNE_FILE_SECONDARY}"

TEST_DIR="${CAAS_OUTBASE}/${TRAIT}/${SOURCE_RUN_SUBDIR}/standalone"
INT_RUNTIME="${CAAS_OUTBASE}/${TRAIT}/${SOURCE_RUN_SUBDIR}/filter"

INPUT_TRAITFILE="${INT_RUNTIME}/data_exploration/2.CT/1.Traitfiles/traitfile.tab"
INPUT_BOOT_TRAITFILE="${INT_RUNTIME}/data_exploration/2.CT/2.Bootstrap_traitfiles/boot_traitfile.tab"
INPUT_DISCOVERY="${INT_RUNTIME}/caastools/discovery.tab"
INPUT_BACKGROUND="${INT_RUNTIME}/caastools/background.output"
INPUT_BACKGROUND_GENES="${INT_RUNTIME}/caastools/background_genes.output"
INPUT_RESAMPLE="${INT_RUNTIME}/caastools/resample.tab"
INPUT_BOOTSTRAP="${INT_RUNTIME}/caastools/bootstrap.tab"
INPUT_GLOBAL_META_CAAS="${INT_RUNTIME}/signification/meta_caas/global_meta_caas.tsv"
INPUT_MASTER_CSV="${INT_RUNTIME}/ct_disambiguation/caas_convergence_master.csv"
INPUT_FILTERED_DISCOVERY="${INT_RUNTIME}/postproc/gene_filtering/filtered_discovery.tsv"
INPUT_CLEANED_BG="${INT_RUNTIME}/postproc/cleaned_backgrounds/cleaned_background_main.txt"
INPUT_ORA_GENE_LISTS_DIR="${INT_RUNTIME}/gene_lists/position"
INPUT_ACCUMULATION_GENE_LISTS_DIR="${INT_RUNTIME}/accumulation/randomization/gene_lists"
SELECTION_INPUTS_DIR="${INT_RUNTIME}/postproc/disambiguation_characterization/gene_relation_analysis/txt/"
INPUT_SELECTION_PP_TOP="${SELECTION_INPUTS_DIR}/all_top.txt"
INPUT_SELECTION_PP_BOTTOM="${SELECTION_INPUTS_DIR}/all_bottom.txt"
RER_STAGED_DIR="${INT_RUNTIME}/RERConverge"
INPUT_RER_TRAIT_OUT="${RER_STAGED_DIR}/RER_Traits/${TRAIT}.polished.output"
INPUT_RER_TREES_OUT="${RER_STAGED_DIR}/RER_Objects/ALL_FEB23_geneTrees.txt.masterTree.output"
INPUT_RER_MATRIX_OUT="${RER_STAGED_DIR}/RER_Objects/${TRAIT}.RERmatrix.output"

SCORING_POSTPROC_INPUT="${INPUT_FILTERED_DISCOVERY}"
SCORING_BACKGROUND_INPUT="${INPUT_CLEANED_BG}"
SCORING_FADE_SUMMARY_TOP="${INT_RUNTIME}/selection/fade/top/fade_summary_top.tsv"
SCORING_FADE_SUMMARY_BOTTOM="${INT_RUNTIME}/selection/fade/bottom/fade_summary_bottom.tsv"
SCORING_FADE_SITE_TOP="${INT_RUNTIME}/selection/fade/top/fade_site_bf_top.tsv"
SCORING_FADE_SITE_BOTTOM="${INT_RUNTIME}/selection/fade/bottom/fade_site_bf_bottom.tsv"
SCORING_MOLERATE_SUMMARY_TOP="${INT_RUNTIME}/selection/molerate/top/molerate_summary_top.tsv"
SCORING_MOLERATE_SUMMARY_BOTTOM="${INT_RUNTIME}/selection/molerate/bottom/molerate_summary_bottom.tsv"
SCORING_RER_INPUT="${INT_RUNTIME}/RERConverge/RER_Results/rerconverge_summary_${TRAIT}.tsv"
SCORING_ACCUM_DIR="${INT_RUNTIME}/accumulation/aggregation"
SCORING_ACCUM_TOP_DIR="${INT_RUNTIME}/accumulation/top/randomization"
SCORING_ACCUM_BOTTOM_DIR="${INT_RUNTIME}/accumulation/bottom/randomization"
SCORING_VEP_TRANSVAR="${INT_RUNTIME}/characterization/vep/transvar.tsv"
SCORING_VEP_PRIMATEAI="${INT_RUNTIME}/characterization/vep/primateai_mapped.tsv"
SCORING_VEP_AA2PROT="${INT_RUNTIME}/characterization/vep/aa2prot_global.csv"

# ─────────────────────────────────────────────────────────────────────────────
# APPTAINER / SINGULARITY HOME MOUNT
# ─────────────────────────────────────────────────────────────────────────────
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true

# Toy-mode flag array (empty when TOY_MODE=false)
SA_TOY_FLAGS=()
[ "$TOY_MODE" = true ] && SA_TOY_FLAGS=(--toy_mode --toy_n "$TOY_N")

# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

section() {
    echo ""
    echo "══════════════════════════════════════════════════════════════"
    echo "  $*"
    echo "══════════════════════════════════════════════════════════════"
}

step() { echo ""; echo "── $* ──"; }
ok()   { echo "  ✓ $*"; }
warn() { echo "  ⚠ $*"; }
fail() { echo "  ✗ ERROR: $*" >&2; exit 1; }

require_file() {
    local f="$1"
    [ -e "$f" ] || fail "Required input not found: $f"
}

find_latest_file_by_name() {
    local name="$1"
    shift || true
    local root latest="" cand t_latest t_cand
    for root in "$@"; do
        [ -d "$root" ] || continue
        cand="$(find "$root" -type f -name "$name" -printf '%T@ %p\n' 2>/dev/null | sort -nr | head -n1 | cut -d' ' -f2-)"
        if [ -n "$cand" ]; then
            if [ -z "$latest" ]; then
                latest="$cand"
            else
                t_latest="$(stat -c %Y "$latest" 2>/dev/null || echo 0)"
                t_cand="$(stat -c %Y "$cand" 2>/dev/null || echo 0)"
                [ "$t_cand" -gt "$t_latest" ] && latest="$cand"
            fi
        fi
    done
    printf '%s\n' "$latest"
}

maybe_clean_work() {
    local workdir="$1"
    if [ "$CLEAN_WORK" = true ] && [ -d "$workdir" ]; then
        rm -rf "$workdir"
        ok "Cleaned work dir: $workdir"
    fi
}

prepare_runtime_dir() {
    local outdir="$1"
    local parent backup
    mkdir -p "$outdir"
    parent="$(dirname "$outdir")"
    shopt -s nullglob dotglob
    local existing=("$outdir"/*)
    shopt -u nullglob dotglob
    if [ "${#existing[@]}" -gt 0 ]; then
        backup="${parent}/${timestamp}_bak"
        mkdir -p "$backup"
        mv "${existing[@]}" "$backup"/
        ok "Moved previous runtime contents to ${backup}"
    fi
}

check_low_contrasts() {
    local outdir="$1"
    local label="$2"
    local sentinel="${outdir}/low_contrasts.skip"
    if [ -f "$sentinel" ]; then
        warn "${label}: Fewer than ${MIN_CONTRASTS:-3} foreground contrasts — CT pipeline was skipped."
        warn "  Detail: $(cat "$sentinel")"
        return 0
    fi
    return 1
}

# ─────────────────────────────────────────────────────────────────────────────
# DEPENDENCY GUARDRAILS (integrated mode)
# ─────────────────────────────────────────────────────────────────────────────
check_integrated_deps() {
    if [ "$INT_CT_DISAMBIGUATION" = true ] && [ "$INT_CT_SIGNIFICATION" != true ]; then
        fail "INT_CT_DISAMBIGUATION requires INT_CT_SIGNIFICATION=true"
    fi
    if [ "$INT_ORA" = true ] && [ "$INT_CT_SIGNIFICATION" != true ]; then
        fail "INT_ORA in integrated mode requires INT_CT_SIGNIFICATION=true (feeds postproc)"
    fi
    if [ "$INT_STRING" = true ] && [ "$INT_ORA" != true ]; then
        fail "INT_STRING requires INT_ORA=true (shared channel)"
    fi
    if [ "$INT_CT_ACCUMULATION" = true ]; then
        [ "$INT_CT_SIGNIFICATION" = true ] || fail "INT_CT_ACCUMULATION requires INT_CT_SIGNIFICATION=true"
        [ "$INT_CT_DISAMBIGUATION" = true ] || fail "INT_CT_ACCUMULATION requires INT_CT_DISAMBIGUATION=true"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# BUILD SHARED FLAG ARRAYS (integrated mode)
# ─────────────────────────────────────────────────────────────────────────────
build_integrated_flags() {
    COMMON_FLAGS=(
        -profile slurm
        --ct_tool ""
        --alignment "$ALI_DIR"
        --ali_format "fasta"
        --my_traits "$INPUT_TRAITS"
        --traitname "$TRAIT"
        --cycles "$CYCLES"
        --caas_config "$INPUT_TRAITFILE"
        --tree "$INPUT_TREE"
        --traitvalues "$INPUT_BOOT_TRAITFILE"
    )

    TOY_FLAGS=()
    [ "$TOY_MODE" = true ] && TOY_FLAGS=(--toy_mode --toy_n "$TOY_N")

    OPTIONAL_BOOL_FLAGS=()

    if [ "$INT_PRUNE_DATA" = true ]; then
        OPTIONAL_BOOL_FLAGS+=(--prune_data --prune_list "$PRUNE_LIST" --prune_list_secondary "$PRUNE_LIST_SECONDARY")
    else
        OPTIONAL_BOOL_FLAGS+=(--prune_data false)
    fi

    [ "$INT_REPORTING" = true ]          && OPTIONAL_BOOL_FLAGS+=(--reporting)          || OPTIONAL_BOOL_FLAGS+=(--reporting false)
    [ "$INT_CONTRAST_SELECTION" = true ] && OPTIONAL_BOOL_FLAGS+=(--contrast_selection) || OPTIONAL_BOOL_FLAGS+=(--contrast_selection false)

    OPTIONAL_BOOL_FLAGS+=(--ct_postproc)

    [ "$INT_CT_SIGNIFICATION" = true ] && OPTIONAL_BOOL_FLAGS+=(--ct_signification) || OPTIONAL_BOOL_FLAGS+=(--ct_signification false)

    if [ "$INT_CT_DISAMBIGUATION" = true ]; then
        OPTIONAL_BOOL_FLAGS+=(
            --ct_disambiguation --ct_signification
            --ct_disambig_asr_mode "$INT_DISAMBIG_ASR_MODE"
            --ct_disambig_asr_cache_dir "$INT_DISAMBIG_ASR_CACHE_DIR"
            --ct_disambig_posterior_threshold 0.1
        )
    fi

    [ "$INT_ORA"    = true ] && OPTIONAL_BOOL_FLAGS+=(--ora)    || OPTIONAL_BOOL_FLAGS+=(--ora false)
    [ "$INT_STRING" = true ] && OPTIONAL_BOOL_FLAGS+=(--string) || OPTIONAL_BOOL_FLAGS+=(--string false)

    if [ "$INT_CT_ACCUMULATION" = true ]; then
        OPTIONAL_BOOL_FLAGS+=(--ct_accumulation --accumulation_n_randomizations "$CYCLES")
    else
        OPTIONAL_BOOL_FLAGS+=(--ct_accumulation false)
    fi

    [ "$INT_FADE"     = true ] && OPTIONAL_BOOL_FLAGS+=(--fade --fade_mode "all")         || OPTIONAL_BOOL_FLAGS+=(--fade false)
    [ "$INT_MOLERATE" = true ] && OPTIONAL_BOOL_FLAGS+=(--molerate --molerate_mode "all") || OPTIONAL_BOOL_FLAGS+=(--molerate false)

    if [ "$INT_RER" = true ]; then
        if [ "$SA_RER_CONTINUOUS_ONLY" = true ]; then
            OPTIONAL_BOOL_FLAGS+=(
                --rer_tool "continuous"
                --rer_gene_set_mode "all"
                --rer_perm_batches "$SA_RER_PERM_BATCHES"
                --rer_perms_per_batch "$SA_RER_PERMS_PER_BATCH"
                --rer_perm_mode "$SA_RER_PERM_MODE"
                --trait_out  "$INPUT_RER_TRAIT_OUT"
                --trees_out  "$INPUT_RER_TREES_OUT"
                --matrix_out "$INPUT_RER_MATRIX_OUT"
            )
        else
            OPTIONAL_BOOL_FLAGS+=(
                --rer_tool "$SA_RER_TOOL"
                --rer_gene_set_mode "all"
                --rer_perm_batches "$SA_RER_PERM_BATCHES"
                --rer_perms_per_batch "$SA_RER_PERMS_PER_BATCH"
                --rer_perm_mode "$SA_RER_PERM_MODE"
                --gene_trees "$INPUT_GENE_TREES"
            )
        fi
    fi

    OPTIONAL_TRAIT_FLAGS=()
    [ "$INT_USE_SECONDARY_TRAIT" = true ] && [ -n "$SECONDARY_TRAIT" ] && OPTIONAL_TRAIT_FLAGS+=(--secondary_trait "$SECONDARY_TRAIT")
    [ "$INT_USE_BRANCH_TRAIT"    = true ] && [ -n "$BRANCH_TRAIT"    ] && OPTIONAL_TRAIT_FLAGS+=(--branch_trait "$BRANCH_TRAIT")
    [ "$INT_USE_N_TRAIT"         = true ] && [ -n "$N_TRAIT"         ] && OPTIONAL_TRAIT_FLAGS+=(--n_trait "$N_TRAIT")
    [ "$INT_USE_C_TRAIT"         = true ] && [ -n "$C_TRAIT"         ] && OPTIONAL_TRAIT_FLAGS+=(--c_trait "$C_TRAIT")

    if [ "$INT_VEP" = true ]; then
        OPTIONAL_BOOL_FLAGS+=(
            --vep
            --vep_cds_dir   "$INPUT_VEP_CDS_DIR"
            --vep_track_dir "$INPUT_VEP_TRACK_DIR"
        )
        [ -n "$INPUT_VEP_PRIMATEAI_DB" ] && OPTIONAL_BOOL_FLAGS+=(--vep_primateai_db "$INPUT_VEP_PRIMATEAI_DB")
    else
        OPTIONAL_BOOL_FLAGS+=(--vep false)
    fi

    if [ "$INT_SCORING" = true ]; then
        OPTIONAL_BOOL_FLAGS+=(--scoring)
        if [ "$INT_SCORING_STRESS" = true ]; then
            OPTIONAL_BOOL_FLAGS+=(--scoring_stress --scoring_stress_top_n "$INT_SCORING_STRESS_TOP_N")
        else
            OPTIONAL_BOOL_FLAGS+=(--scoring_stress false)
        fi
    else
        OPTIONAL_BOOL_FLAGS+=(--scoring false)
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# INTEGRATED RUN
# ─────────────────────────────────────────────────────────────────────────────
run_integrated() {
    local mode="$1"
    local outdir
    [ "$mode" = "filter" ] && outdir="$INT_RUNTIME" || outdir="${TEST_DIR}/integrated_${mode}/runtime"
    local workdir="${outdir}/work"
    local resume_flags=()

    if [ "$INT_RESUME" = true ]; then
        resume_flags=(-resume)
    else
        prepare_runtime_dir "$outdir"
    fi
    mkdir -p "$outdir" "$workdir"
    build_integrated_flags

    echo ""
    echo "Running integrated pipeline — mode: ${mode}"
    echo "  Output: ${outdir}"
    echo ""

    nextflow run main.nf \
        -with-tower \
        "${resume_flags[@]+"${resume_flags[@]}"}" \
        "${COMMON_FLAGS[@]}" \
        "${OPTIONAL_BOOL_FLAGS[@]}" \
        "${OPTIONAL_TRAIT_FLAGS[@]}" \
        "${TOY_FLAGS[@]+"${TOY_FLAGS[@]}"}" \
        -w "$workdir" \
        --outdir "$outdir" \
        --caas_postproc_mode "$mode" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Integrated ${mode} completed → ${outdir}"
    check_low_contrasts "$outdir" "integrated_${mode}" || true
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE CT
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_ct() {
    local outdir="${TEST_DIR}/standalone_ct/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$ALI_DIR"; require_file "$INPUT_TRAITS"; require_file "$INPUT_TREE"; require_file "$INPUT_TRAITFILE"

    echo "Running standalone CT (discovery + resample + bootstrap)"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --ct_tool "discovery,resample,bootstrap" \
        --alignment "$ALI_DIR" --ali_format "fasta" \
        --my_traits "$INPUT_TRAITS" \
        --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        ${BRANCH_TRAIT:+--branch_trait "$BRANCH_TRAIT"} \
        ${N_TRAIT:+--n_trait "$N_TRAIT"} \
        ${C_TRAIT:+--c_trait "$C_TRAIT"} \
        --caas_config "$INPUT_TRAITFILE" \
        --tree "$INPUT_TREE" \
        --traitvalues "$INPUT_BOOT_TRAITFILE" \
        --cycles "$CYCLES" \
        "${SA_TOY_FLAGS[@]+"${SA_TOY_FLAGS[@]}"}" \
        --ct_postproc false --ct_signification false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false \
        --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone CT completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE SIGNIFICATION
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_signification() {
    local outdir="${TEST_DIR}/standalone_signification/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_DISCOVERY"; require_file "$INPUT_RESAMPLE"
    require_file "$INPUT_BOOTSTRAP"; require_file "$INPUT_BACKGROUND_GENES"

    echo "Running standalone Signification from file inputs"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --ct_signification \
        --discovery_out "$INPUT_DISCOVERY" \
        --bootstrap_input "$INPUT_BOOTSTRAP" \
        --background_input "$INPUT_BACKGROUND_GENES" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        --tree "$INPUT_TREE" --caas_config "$INPUT_TRAITFILE" \
        --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false \
        --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone Signification completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE DISAMBIGUATION
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_disambiguation() {
    require_file "$INPUT_GLOBAL_META_CAAS"; require_file "$INPUT_TRAITFILE"; require_file "$INPUT_TREE"

    if [ "$SA_DISAMBIG_COMPUTE" = true ]; then
        local outdir="${TEST_DIR}/standalone_disambiguation_compute/${timestamp}"
        local workdir="${outdir}/work"
        mkdir -p "$outdir" "$workdir"
        echo "  ASR mode: compute  →  cache dir: $INPUT_ASR_CACHE"

        nextflow run main.nf \
            -with-tower \
            -profile slurm \
            --ct_signification false --ct_disambiguation \
            --ct_disambig_caas_metadata "$INPUT_GLOBAL_META_CAAS" \
            --alignment "$ALI_DIR" --ali_format "fasta" \
            --caas_config "$INPUT_TRAITFILE" \
            --tree "$INPUT_TREE" --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
            --ct_disambig_asr_mode "compute" \
            --ct_disambig_asr_cache_dir "$INPUT_ASR_CACHE" \
            --ct_disambig_posterior_threshold 0.1 \
            --ct_postproc false --ora false --string false \
            --ct_accumulation false --reporting false --contrast_selection false \
            -w "$workdir" --outdir "$outdir" \
            -with-report "${outdir}/pipeline_info/execution_report.html" \
            -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

        ok "Standalone Disambiguation (compute) completed → ${outdir}"
        maybe_clean_work "$workdir"
    fi

    if [ "$SA_DISAMBIG_PRECOMPUTED" = true ]; then
        local outdir="${TEST_DIR}/standalone_disambiguation_precomputed/${timestamp}"
        local workdir="${outdir}/work"
        mkdir -p "$outdir" "$workdir"
        echo "  ASR mode: precomputed  →  cache dir: $INPUT_ASR_CACHED"

        nextflow run main.nf \
            -with-tower \
            -profile slurm \
            --ct_signification false --ct_disambiguation \
            --ct_disambig_caas_metadata "$INPUT_GLOBAL_META_CAAS" \
            --alignment "$ALI_DIR" --ali_format "fasta" \
            --caas_config "$INPUT_TRAITFILE" \
            --tree "$INPUT_TREE" --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
            --ct_disambig_asr_mode "precomputed" \
            --ct_disambig_asr_cache_dir "$INPUT_ASR_CACHED" \
            --ct_disambig_posterior_threshold 0.1 \
            --ct_postproc false --ora false --string false \
            --ct_accumulation false --reporting false --contrast_selection false \
            -w "$workdir" --outdir "$outdir" \
            -with-report "${outdir}/pipeline_info/execution_report.html" \
            -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

        ok "Standalone Disambiguation (precomputed) completed → ${outdir}"
        maybe_clean_work "$workdir"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE POSTPROC
# ─────────────────────────────────────────────────────────────────────────────
_run_standalone_postproc() {
    local mode="$1"
    local outdir="${TEST_DIR}/standalone_postproc_${mode}/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_MASTER_CSV"; require_file "$INPUT_BACKGROUND_GENES"

    echo "Running standalone CT Post-processing (${mode} mode)"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --ct_signification false --ct_disambiguation false --ct_postproc \
        --caas_postproc_mode "$mode" \
        --disambiguation_input "$INPUT_MASTER_CSV" \
        --background_input "$INPUT_BACKGROUND_GENES" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        --tree "$INPUT_TREE" --caas_config "$INPUT_TRAITFILE" \
        --ora false --string false --ct_accumulation false \
        --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone PostProc (${mode}) completed → ${outdir}"
    maybe_clean_work "$workdir"
}

run_standalone_postproc_filter()      { _run_standalone_postproc filter;      }
run_standalone_postproc_exploratory() { _run_standalone_postproc exploratory; }

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE ORA (+ optional STRING)
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_ora() {
    local outdir="${TEST_DIR}/standalone_ora/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_ORA_GENE_LISTS_DIR"; require_file "$INPUT_CLEANED_BG"

    local string_flag="--string false"
    [ "$RUN_SA_STRING" = true ] && string_flag="--string"

    echo "Running standalone ORA  (STRING: $RUN_SA_STRING)"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --ct_signification false --ct_disambiguation false --ct_postproc false \
        --ora $string_flag \
        --ora_gene_lists_input "$INPUT_ORA_GENE_LISTS_DIR" \
        --ora_background_input "$INPUT_CLEANED_BG" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        --tree "$INPUT_TREE" --caas_config "$INPUT_TRAITFILE" \
        --ct_accumulation false --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone ORA completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE CT ACCUMULATION
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_accumulation() {
    local outdir="${TEST_DIR}/standalone_accumulation/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_FILTERED_DISCOVERY"; require_file "$INPUT_CLEANED_BG"
    require_file "$INPUT_TRAITFILE"; require_file "$ALI_DIR"

    echo "Running standalone CT Accumulation"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --ct_signification false --ct_disambiguation false --ct_postproc false \
        --ora false --string false --ct_accumulation \
        --accumulation_n_randomizations "$CYCLES" \
        --accumulation_caas_input       "$INPUT_FILTERED_DISCOVERY" \
        --accumulation_background_input "$INPUT_CLEANED_BG" \
        --alignment "$ALI_DIR" --ali_format "fasta" \
        --caas_config "$INPUT_TRAITFILE" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        --tree "$INPUT_TREE" \
        "${SA_TOY_FLAGS[@]+"${SA_TOY_FLAGS[@]}"}" \
        --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone CT Accumulation completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE REPORTING
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_reporting() {
    local outdir="${TEST_DIR}/standalone_reporting/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_TRAITS"; require_file "$INPUT_TREE"

    echo "Running standalone Reporting"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --reporting \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        ${BRANCH_TRAIT:+--branch_trait "$BRANCH_TRAIT"} \
        ${N_TRAIT:+--n_trait "$N_TRAIT"} \
        ${C_TRAIT:+--c_trait "$C_TRAIT"} \
        --tree "$INPUT_TREE" \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone Reporting completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE CONTRAST SELECTION
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_contrast_selection() {
    local outdir="${TEST_DIR}/standalone_contrast_selection/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_TRAITS"; require_file "$INPUT_TREE"

    echo "Running standalone Contrast Selection"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --contrast_selection \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        ${BRANCH_TRAIT:+--branch_trait "$BRANCH_TRAIT"} \
        ${N_TRAIT:+--n_trait "$N_TRAIT"} \
        ${C_TRAIT:+--c_trait "$C_TRAIT"} \
        --tree "$INPUT_TREE" \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --reporting false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone Contrast Selection completed → ${outdir}"
    check_low_contrasts "$outdir" "standalone_contrast_selection" || true
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE VEP
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_vep() {
    local outdir="${TEST_DIR}/standalone_vep/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"
    require_file "$INPUT_FILTERED_DISCOVERY"
    require_file "$INPUT_VEP_CDS_DIR"
    require_file "$INPUT_VEP_TRACK_DIR"

    local primateai_db="${INPUT_VEP_PRIMATEAI_DB:-}"

    echo "Running standalone VEP"
    echo "  Output: ${outdir}"

    local vep_flags=(
        --vep
        --vep_caas_input "$INPUT_FILTERED_DISCOVERY"
        --vep_cds_dir    "$INPUT_VEP_CDS_DIR"
        --vep_track_dir  "$INPUT_VEP_TRACK_DIR"
    )
    [ -n "$primateai_db" ] && vep_flags+=(--vep_primateai_db "$primateai_db")

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        "${vep_flags[@]}" \
        --prune_data false --reporting false \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone VEP completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE SCORING
# ─────────────────────────────────────────────────────────────────────────────
run_standalone_scoring() {
    local outdir="${TEST_DIR}/standalone_scoring/${timestamp}"
    local workdir="${outdir}/work"
    local accum_stage_dir="${outdir}/inputs_accumulation"
    mkdir -p "$outdir" "$workdir" "$accum_stage_dir"
    require_file "$SCORING_POSTPROC_INPUT"; require_file "$SCORING_BACKGROUND_INPUT"

    # Auto-discover FADE per-site files when preferred paths are absent
    local scoring_fade_site_top="$SCORING_FADE_SITE_TOP"
    local scoring_fade_site_bottom="$SCORING_FADE_SITE_BOTTOM"
    [ ! -f "$scoring_fade_site_top"    ] && scoring_fade_site_top="$(find_latest_file_by_name "fade_site_bf_top.tsv"    "$INT_RUNTIME" "${TEST_DIR}/runs")"
    [ ! -f "$scoring_fade_site_bottom" ] && scoring_fade_site_bottom="$(find_latest_file_by_name "fade_site_bf_bottom.tsv" "$INT_RUNTIME" "${TEST_DIR}/runs")"

    echo "Running standalone Scoring"
    echo "  Output: ${outdir}"

    # Stage accumulation files
    local accum_staged_n=0 scheme srcf dstf
    for scheme in us gs4 gs3 gs2 gs1 full_pool gs0; do
        for dir in top bottom; do
            dstf="${accum_stage_dir}/accumulation_${dir}_${scheme}_aggregated_results.csv"
            srcf="${SCORING_ACCUM_TOP_DIR}/accumulation_${dir}_${scheme}_aggregated_results.csv"
            [ "$dir" = "bottom" ] && srcf="${SCORING_ACCUM_BOTTOM_DIR}/accumulation_${dir}_${scheme}_aggregated_results.csv"
            if [ -f "$srcf" ]; then
                cp -f "$srcf" "$dstf"
                accum_staged_n=$((accum_staged_n + 1))
            else
                srcf="${srcf/_${dir}_/_}"   # legacy filename without dir tag
                [ -f "$srcf" ] && cp -f "$srcf" "$dstf" && accum_staged_n=$((accum_staged_n + 1))
            fi
        done
    done
    [ "$accum_staged_n" -eq 0 ] \
        && warn "No accumulation CSVs staged — gene_rand_score will be skipped." \
        || ok "Staged ${accum_staged_n} accumulation CSV(s)."

    local optional_scoring_flags=()
    [ -f "$SCORING_FADE_SUMMARY_TOP"    ] && optional_scoring_flags+=(--scoring_fade_summary_top    "$SCORING_FADE_SUMMARY_TOP")
    [ -f "$SCORING_FADE_SUMMARY_BOTTOM" ] && optional_scoring_flags+=(--scoring_fade_summary_bottom "$SCORING_FADE_SUMMARY_BOTTOM")
    [ -n "$scoring_fade_site_top"    ] && [ -f "$scoring_fade_site_top"    ] && optional_scoring_flags+=(--scoring_fade_site_top    "$scoring_fade_site_top")
    [ -n "$scoring_fade_site_bottom" ] && [ -f "$scoring_fade_site_bottom" ] && optional_scoring_flags+=(--scoring_fade_site_bottom "$scoring_fade_site_bottom")
    [ -f "$SCORING_VEP_TRANSVAR"        ] && optional_scoring_flags+=(--scoring_vep_transvar   "$SCORING_VEP_TRANSVAR")
    [ -f "$SCORING_VEP_PRIMATEAI"       ] && optional_scoring_flags+=(--scoring_vep_primateai  "$SCORING_VEP_PRIMATEAI")
    [ -f "$SCORING_VEP_AA2PROT"         ] && optional_scoring_flags+=(--scoring_vep_aa2prot    "$SCORING_VEP_AA2PROT")
    [ "$SA_SCORING_STRESS" = true       ] && optional_scoring_flags+=(--scoring_stress --scoring_stress_top_n "$SA_SCORING_STRESS_TOP_N")

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --scoring \
        --scoring_postproc_input            "$SCORING_POSTPROC_INPUT" \
        --scoring_background_input          "$SCORING_BACKGROUND_INPUT" \
        --scoring_molerate_summary_top      "$SCORING_MOLERATE_SUMMARY_TOP" \
        --scoring_molerate_summary_bottom   "$SCORING_MOLERATE_SUMMARY_BOTTOM" \
        --scoring_rer_input                 "$SCORING_RER_INPUT" \
        --scoring_accum_dir                 "$accum_stage_dir" \
        "${optional_scoring_flags[@]}" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        --tree "$INPUT_TREE" \
        --prune_data false --reporting false \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone Scoring completed → ${outdir}"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE SELECTION HELPERS (FADE / MoleRate / RERConverge)
# ─────────────────────────────────────────────────────────────────────────────
_selection_mode_label() {
    case "$1" in all) echo "all" ;; gene_set_pp) echo "geneset_pp" ;; *) fail "Unknown selection mode: $1" ;; esac
}
_selection_mode_value() {
    case "$1" in all) echo "all" ;; gene_set_pp) echo "gene_set" ;; *) fail "Unknown selection mode: $1" ;; esac
}
_require_selection_mode_inputs() {
    case "$1" in
        all) ;;
        gene_set_pp)
            require_file "$INPUT_SELECTION_PP_TOP"
            require_file "$INPUT_SELECTION_PP_BOTTOM"
            ;;
        *) fail "Unknown selection mode: $1" ;;
    esac
}
_selection_source_flags() {
    local tool="$1" mode="$2"
    SELECTION_SOURCE_FLAGS=()
    [ "$mode" = "gene_set_pp" ] && SELECTION_SOURCE_FLAGS=(
        "--postproc_top"    "$INPUT_SELECTION_PP_TOP"
        "--postproc_bottom" "$INPUT_SELECTION_PP_BOTTOM"
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE FADE
# ─────────────────────────────────────────────────────────────────────────────
_run_standalone_fade_mode() {
    local mode="$1"
    local mode_label; mode_label="$(_selection_mode_label "$mode")"
    local mode_value; mode_value="$(_selection_mode_value "$mode")"
    _require_selection_mode_inputs "$mode"
    _selection_source_flags "fade" "$mode"

    local outdir="${TEST_DIR}/standalone_fade_${mode_label}/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"

    echo "Running standalone FADE (mode: ${mode_label})"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --fade --fade_mode "$mode_value" \
        "${SELECTION_SOURCE_FLAGS[@]+"${SELECTION_SOURCE_FLAGS[@]}"}" \
        --alignment "$ALI_DIR" --ali_format "fasta" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        --tree "$INPUT_TREE" --reporting \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone FADE (${mode_label}) completed → ${outdir}"
    maybe_clean_work "$workdir"
}

run_standalone_fade() {
    [ "$SA_SELECTION_MODE_ALL"      = true ] && _run_standalone_fade_mode all
    [ "$SA_SELECTION_MODE_GENESET_PP" = true ] && _run_standalone_fade_mode gene_set_pp
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE MOLERATE
# ─────────────────────────────────────────────────────────────────────────────
_run_standalone_molerate_mode() {
    local mode="$1"
    local mode_label; mode_label="$(_selection_mode_label "$mode")"
    local mode_value; mode_value="$(_selection_mode_value "$mode")"
    _require_selection_mode_inputs "$mode"
    _selection_source_flags "molerate" "$mode"

    local outdir="${TEST_DIR}/standalone_molerate_${mode_label}/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"

    echo "Running standalone MoleRate (mode: ${mode_label})"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        --molerate --molerate_mode "$mode_value" \
        "${SELECTION_SOURCE_FLAGS[@]+"${SELECTION_SOURCE_FLAGS[@]}"}" \
        --alignment "$ALI_DIR" --ali_format "fasta" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        --tree "$INPUT_TREE" --reporting \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone MoleRate (${mode_label}) completed → ${outdir}"
    maybe_clean_work "$workdir"
}

run_standalone_molerate() {
    [ "$SA_SELECTION_MODE_ALL"        = true ] && _run_standalone_molerate_mode all
    [ "$SA_SELECTION_MODE_GENESET_PP" = true ] && _run_standalone_molerate_mode gene_set_pp
}

# ─────────────────────────────────────────────────────────────────────────────
# STANDALONE RERConverge
# ─────────────────────────────────────────────────────────────────────────────
_run_standalone_rer_mode() {
    local mode="$1"
    local mode_label; mode_label="$(_selection_mode_label "$mode")"
    local mode_value; mode_value="$(_selection_mode_value "$mode")"
    _require_selection_mode_inputs "$mode"
    _selection_source_flags "rer" "$mode"

    local rer_tool_flags=(
        --rer_gene_set_mode "$mode_value"
        --rer_perm_batches "$SA_RER_PERM_BATCHES"
        --rer_perms_per_batch "$SA_RER_PERMS_PER_BATCH"
        --rer_perm_mode "$SA_RER_PERM_MODE"
    )
    [ -n "${SA_RER_GMT_FILE:-}" ] && rer_tool_flags+=(--rer_gmt_file "$SA_RER_GMT_FILE")

    if [ "$SA_RER_CONTINUOUS_ONLY" = true ]; then
        rer_tool_flags+=(
            --rer_tool "continuous"
            --trait_out  "$INPUT_RER_TRAIT_OUT"
            --trees_out  "$INPUT_RER_TREES_OUT"
            --matrix_out "$INPUT_RER_MATRIX_OUT"
        )
    else
        rer_tool_flags+=(
            --rer_tool "$SA_RER_TOOL"
            --gene_trees "$INPUT_GENE_TREES"
        )
    fi

    local outdir="${TEST_DIR}/standalone_rer_${mode_label}/${timestamp}"
    local workdir="${outdir}/work"
    mkdir -p "$outdir" "$workdir"

    echo "Running standalone RERConverge (mode: ${mode_label})"
    echo "  Output: ${outdir}"

    nextflow run main.nf \
        -with-tower \
        -profile slurm \
        "${rer_tool_flags[@]}" \
        "${SELECTION_SOURCE_FLAGS[@]+"${SELECTION_SOURCE_FLAGS[@]}"}" \
        --my_traits "$INPUT_TRAITS" --traitname "$TRAIT" \
        ${SECONDARY_TRAIT:+--secondary_trait "$SECONDARY_TRAIT"} \
        ${BRANCH_TRAIT:+--branch_trait "$BRANCH_TRAIT"} \
        ${N_TRAIT:+--n_trait "$N_TRAIT"} \
        ${C_TRAIT:+--c_trait "$C_TRAIT"} \
        --tree "$INPUT_TREE" \
        --ct_signification false --ct_postproc false --ct_disambiguation false \
        --ora false --string false --ct_accumulation false \
        --reporting false --contrast_selection false \
        -w "$workdir" --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Standalone RERConverge (${mode_label}) completed → ${outdir}"
    maybe_clean_work "$workdir"
}

run_standalone_rer() {
    [ "$SA_SELECTION_MODE_ALL"        = true ] && _run_standalone_rer_mode all
    [ "$SA_SELECTION_MODE_GENESET_PP" = true ] && _run_standalone_rer_mode gene_set_pp
}

# ─────────────────────────────────────────────────────────────────────────────
# PRINT CONFIGURATION SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print_config() {
    section "PHYLOPHERE STRESS-TEST CONFIGURATION"
    echo ""
    echo "Timestamp          : $timestamp"
    echo "TRAIT              : $TRAIT  [CLASS ${TRAIT_CLASS}]"
    echo "SECONDARY_TRAIT    : ${SECONDARY_TRAIT:-<none>}"
    echo "N_TRAIT            : ${N_TRAIT:-<none>}"
    echo "C_TRAIT            : ${C_TRAIT:-<none>}"
    echo "BRANCH_TRAIT       : ${BRANCH_TRAIT:-<none>}"
    echo "SOURCE_RUN_SUBDIR  : $SOURCE_RUN_SUBDIR"
    echo "INT_RUNTIME        : $INT_RUNTIME"
    echo "Alignment dir      : $ALI_DIR"
    echo "Bootstrap cycles   : $CYCLES"
    echo "Toy mode           : $TOY_MODE  (N=$TOY_N)"
    echo ""
    echo "── Infrastructure ──────────────────────────"
    echo "  CLEAN_WORK              : $CLEAN_WORK"
    echo ""
    echo "── Integrated runs ─────────────────────────"
    echo "  RUN_INT_FILTER          : $RUN_INT_FILTER"
    echo "  RUN_INT_EXPLORATORY     : $RUN_INT_EXPLORATORY"
    echo "  sub: INT_PRUNE_DATA           : $INT_PRUNE_DATA"
    echo "  sub: INT_REPORTING            : $INT_REPORTING"
    echo "  sub: INT_CONTRAST_SELECTION   : $INT_CONTRAST_SELECTION"
    echo "  sub: INT_CT_SIGNIFICATION     : $INT_CT_SIGNIFICATION"
    echo "  sub: INT_CT_DISAMBIGUATION    : $INT_CT_DISAMBIGUATION"
    echo "  sub: INT_ORA                  : $INT_ORA"
    echo "  sub: INT_STRING               : $INT_STRING"
    echo "  sub: INT_CT_ACCUMULATION      : $INT_CT_ACCUMULATION"
    echo "  sub: INT_VEP                  : $INT_VEP"
    echo "  sub: INT_FADE                 : $INT_FADE"
    echo "  sub: INT_MOLERATE             : $INT_MOLERATE"
    echo "  sub: INT_RER                  : $INT_RER"
    echo "  sub: INT_SCORING              : $INT_SCORING"
    echo "    (+ INT_SCORING_STRESS       : $INT_SCORING_STRESS)"
    echo "    (+ INT_SCORING_STRESS_TOP_N : $INT_SCORING_STRESS_TOP_N)"
    echo ""
    echo "── Standalone runs ─────────────────────────"
    echo "  RUN_SA_CT                     : $RUN_SA_CT"
    echo "  RUN_SA_SIGNIFICATION          : $RUN_SA_SIGNIFICATION"
    echo "  RUN_SA_DISAMBIGUATION         : $RUN_SA_DISAMBIGUATION"
    echo "    (+ SA_DISAMBIG_COMPUTE      : $SA_DISAMBIG_COMPUTE)"
    echo "    (+ SA_DISAMBIG_PRECOMPUTED  : $SA_DISAMBIG_PRECOMPUTED)"
    echo "  RUN_SA_POSTPROC_FILTER        : $RUN_SA_POSTPROC_FILTER"
    echo "  RUN_SA_POSTPROC_EXPLORATORY   : $RUN_SA_POSTPROC_EXPLORATORY"
    echo "  RUN_SA_ORA                    : $RUN_SA_ORA"
    echo "    (+ STRING via RUN_SA_STRING : $RUN_SA_STRING)"
    echo "  RUN_SA_ACCUMULATION           : $RUN_SA_ACCUMULATION"
    echo "  RUN_SA_REPORTING              : $RUN_SA_REPORTING"
    echo "  RUN_SA_CONTRAST_SELECTION     : $RUN_SA_CONTRAST_SELECTION"
    echo "  RUN_SA_VEP                    : $RUN_SA_VEP"
    echo "  RUN_SA_SCORING                : $RUN_SA_SCORING"
    echo "    (+ SA_SCORING_STRESS        : $SA_SCORING_STRESS)"
    echo "    (+ SA_SCORING_STRESS_TOP_N  : $SA_SCORING_STRESS_TOP_N)"
    echo "  RUN_SA_FADE                   : $RUN_SA_FADE"
    echo "  RUN_SA_MOLERATE               : $RUN_SA_MOLERATE"
    echo "  RUN_SA_RER                    : $RUN_SA_RER"
    echo "    (+ SA_SELECTION_MODE_ALL        : $SA_SELECTION_MODE_ALL)"
    echo "    (+ SA_SELECTION_MODE_GENESET_PP : $SA_SELECTION_MODE_GENESET_PP)"
    echo "    (+ SA_RER_CONTINUOUS_ONLY       : $SA_RER_CONTINUOUS_ONLY)"
    echo "    (+ SA_RER_TOOL                  : $SA_RER_TOOL)"
    echo "    (+ SA_RER_PERM_BATCHES          : $SA_RER_PERM_BATCHES)"
    echo "    (+ SA_RER_PERMS_PER_BATCH       : $SA_RER_PERMS_PER_BATCH)"
    echo "    (+ SA_RER_PERM_MODE             : $SA_RER_PERM_MODE)"
    echo "    (+ SA_RER_GMT_FILE              : ${SA_RER_GMT_FILE:-<none>})"
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# RESULT SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print_summary() {
    section "STRESS-TEST COMPLETE — ${TRAIT}"
    echo ""
    echo "Results root: ${TEST_DIR}/"
    echo ""
    echo "Run directories created:"
    find "${TEST_DIR}/runs" -maxdepth 2 -mindepth 2 -type d 2>/dev/null | sort | while read -r d; do
        echo "  $d"
    done
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
main() {
    print_config

    if [ "$RUN_INT_FILTER" = true ] || [ "$RUN_INT_EXPLORATORY" = true ]; then
        check_integrated_deps
    fi

    [ "$RUN_INT_FILTER"      = true ] && { section "INTEGRATED PIPELINE — filter mode";      run_integrated filter; }
    [ "$RUN_INT_EXPLORATORY" = true ] && { section "INTEGRATED PIPELINE — exploratory mode"; run_integrated exploratory; }

    [ "$RUN_SA_CT"                  = true ] && { section "STANDALONE: CT";                               run_standalone_ct; }
    [ "$RUN_SA_SIGNIFICATION"       = true ] && { section "STANDALONE: Signification";                    run_standalone_signification; }
    [ "$RUN_SA_DISAMBIGUATION"      = true ] && { section "STANDALONE: Disambiguation";                   run_standalone_disambiguation; }
    [ "$RUN_SA_POSTPROC_FILTER"     = true ] && { section "STANDALONE: PostProc — filter mode";           run_standalone_postproc_filter; }
    [ "$RUN_SA_POSTPROC_EXPLORATORY" = true ] && { section "STANDALONE: PostProc — exploratory mode";      run_standalone_postproc_exploratory; }
    [ "$RUN_SA_ORA"                 = true ] && { section "STANDALONE: ORA enrichment";                   run_standalone_ora; }
    [ "$RUN_SA_ACCUMULATION"        = true ] && { section "STANDALONE: CT Accumulation";                  run_standalone_accumulation; }
    [ "$RUN_SA_REPORTING"           = true ] && { section "STANDALONE: Reporting";                        run_standalone_reporting; }
    [ "$RUN_SA_CONTRAST_SELECTION"  = true ] && { section "STANDALONE: Contrast Selection";               run_standalone_contrast_selection; }
    [ "$RUN_SA_VEP"                 = true ] && { section "STANDALONE: VEP";                              run_standalone_vep; }
    [ "$RUN_SA_SCORING"             = true ] && { section "STANDALONE: Composite CAAS Scoring";           run_standalone_scoring; }
    [ "$RUN_SA_FADE"                = true ] && { section "STANDALONE: FADE";                             run_standalone_fade; }
    [ "$RUN_SA_MOLERATE"            = true ] && { section "STANDALONE: MoleRate";                         run_standalone_molerate; }
    [ "$RUN_SA_RER"                 = true ] && { section "STANDALONE: RERConverge";                      run_standalone_rer; }

    print_summary
}

main "$@"
