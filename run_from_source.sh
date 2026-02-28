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
# PHYLOPHERE: Source-Run Module Runner
#
# PURPOSE
# -------
# Run any individual pipeline module (or combination of modules) using
# pre-existing outputs from a completed source run as inputs.  Everything is
# derived from a single SOURCE_RUN directory; individual paths can be
# overridden in the "MANUAL OVERRIDES" block when the source layout differs.
#
# USAGE
# -----
# 1. Set SOURCE_RUN (and optionally SOURCE_FILTER_SUBDIR) to point at the
#    finished run you want to re-use.
# 2. Fill in "PHENOTYPE / DATA SETTINGS" for your trait.
# 3. Flip the module toggles ON/OFF in "MODULE TOGGLES".
# 4. Adjust any module-specific settings in "MODULE SETTINGS".
# 5. Run:  bash run_from_source.sh
# 6. Results appear under OUTPUT_BASE/<timestamp>/
#
# MODULE TOGGLES (quick reference)
# ---------------------------------
#   RUN_CT              – CT discovery + resample + bootstrap from alignments
#   RUN_SIGNIFICATION   – signification from CT outputs
#   RUN_DISAMBIGUATION  – ASR disambiguation (compute or precomputed)
#   RUN_POSTPROC        – postproc in filter or exploratory mode
#   RUN_ORA             – WebGestalt ORA from postproc gene lists
#   RUN_STRING          – STRING PPI enrichment (alongside ORA)
#   RUN_ACCUMULATION    – CT accumulation from filtered_discovery + background
#   RUN_ACCUMULATION_ORA– ORA on accumulation gene lists
#   RUN_FADE            – FADE selection from gene sets
#   RUN_MOLERATE        – MoleRate selection from gene sets
#   RUN_RER             – RERConverge from gene sets + gene trees
#   RUN_REPORTING       – Rmarkdown dataset / phenotype reports
#   RUN_CONTRAST_SELECTION – foreground/background contrast selection
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: run_scripts/run_from_source.sh
#

set -Eeuo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# ENVIRONMENT ACTIVATION
# ─────────────────────────────────────────────────────────────────────────────
echo "Activating environment: phylophere"
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if command -v micromamba &>/dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
elif command -v conda &>/dev/null; then
    conda activate phylophere || echo "Warning: could not activate conda environment"
else
    echo "Warning: neither micromamba nor conda found — continuing without activation"
fi

export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true

# Repo root — this script lives inside run_scripts/
REPO_DIR="/home/miguel/IBE-UPF/PhD/PhyloPhere"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   SOURCE RUN  —  the completed run whose outputs you want to re-use
# ══════════════════════════════════════════════════════════════════════════════
# Set this to the run-root (the timestamped directory).  The sub-directory
# named by SOURCE_FILTER_SUBDIR is appended automatically to locate the
# standard pipeline outputs (caastools/, postproc/, accumulation/, …).
# The full resulting path is SOURCE_BASE.
# ─────────────────────────────────────────────────────────────────────────────

# External drive mount point
EXT_DRIVE="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92"

# Path to the timestamped run root
SOURCE_RUN_ROOT="${EXT_DRIVE}/CAAS_2.0/Results/MARVIN_RESULTS/Ethanol/20260224_195959"

# Sub-directory under SOURCE_RUN_ROOT where the structured outputs live.
# Typical values: "filter" (default) or "exploratory".  Leave empty ("") if
# the run root IS the outputs directory.
SOURCE_FILTER_SUBDIR="filter"

# Derived path — do not edit unless you override paths individually below
if [ -n "$SOURCE_FILTER_SUBDIR" ]; then
    SOURCE_BASE="${SOURCE_RUN_ROOT}/${SOURCE_FILTER_SUBDIR}"
else
    SOURCE_BASE="${SOURCE_RUN_ROOT}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   PHENOTYPE / DATA SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

# Primary trait name (must match a column in the trait CSV)
TRAIT="Ethanol"

# Optional traits — set to "" to disable; they are passed only when non-empty
SECONDARY_TRAIT=""
N_TRAIT=""
C_TRAIT=""
BRANCH_TRAIT=""

# Main trait CSV (columns: species + trait columns)
TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv"

# Phylogenetic tree (Newick / NEXUS)
TREE_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree"

# ── Toy mode ──────────────────────────────────────────────────────────────
# IS_TOY=true  → uses toy alignments + small cycle counts for quick testing.
# IS_TOY=false → uses full paths and full cycle counts (production).
IS_TOY=true

# Toy alignment directory (small gene subset, fast runs)
ALI_DIR_TOY="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Ali_toy"

# Toy bootstrap / randomization counts
CYCLES_TOY="100"
N_RANDOMIZATIONS_TOY="1000"

# Full alignment directory (used when IS_TOY=false)
ALI_DIR_FULL="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Primate_alignments"

# Resolved at runtime by the toy-resolution block below — do not edit directly
ALI_DIR="$ALI_DIR_FULL"

# Gene trees file (used by RERConverge only)
GENE_TREES="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"

# Species prune lists — set to "" to disable pruning
PRUNE_LIST=""
PRUNE_LIST_SECONDARY=""

# Bootstrap / randomization settings (full run; overridden by IS_TOY)
CYCLES="10000"
N_RANDOMIZATIONS="1000000"

# ASR cache directory — used by disambiguation modules
ASR_CACHE_DIR="${EXT_DRIVE}/asr"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   AUTO-DERIVED INPUT PATHS  (from SOURCE_BASE)
#   Override any individual variable below if your source layout differs.
# ══════════════════════════════════════════════════════════════════════════════

# CT traitfiles (discovery config + bootstrap traitfile)
CAAS_CONFIG="${SOURCE_BASE}/data_exploration/2.CT/1.Traitfiles/traitfile.tab"
BOOT_CONFIG="${SOURCE_BASE}/data_exploration/2.CT/2.Bootstrap_traitfiles/boot_traitfile.tab"

# Raw caastools outputs
DISCOVERY_INPUT="${SOURCE_BASE}/caastools/discovery.tab"
DISCOVERY_OUT="${SOURCE_BASE}/caastools/discovery.tab"          # alias used by some modules
BACKGROUND_INPUT="${SOURCE_BASE}/caastools/background_genes.output"
BACKGROUND_RAW="${SOURCE_BASE}/caastools/background.output"
RESAMPLE_INPUT="${SOURCE_BASE}/caastools/resample.tab"
BOOTSTRAP_INPUT="${SOURCE_BASE}/caastools/bootstrap.tab"

# Signification outputs
GLOBAL_META_CAAS="${SOURCE_BASE}/signification/meta_caas/global_meta_caas.tsv"
ORA_GENE_LISTS_DIR="${SOURCE_BASE}/signification/gene_lists/global"

# Disambiguation output
MASTER_CSV="${SOURCE_BASE}/ct_disambiguation/ct_disambiguation/caas_convergence_master.csv"

# Gene filtering / postproc outputs
FILTERED_DISCOVERY="${SOURCE_BASE}/postproc/gene_filtering/filtered_discovery.tsv"
CLEANED_BACKGROUND="${SOURCE_BASE}/postproc/cleaned_backgrounds/cleaned_background_main.txt"
CLEANED_BACKGROUND_US="${SOURCE_BASE}/postproc/cleaned_backgrounds/cleaned_background_US.txt"

# Accumulation outputs
ACCUMULATION_CSV="${SOURCE_BASE}/accumulation/aggregation/accumulation_global.csv"
ACCUMULATION_GENE_LISTS_DIR="${SOURCE_BASE}/accumulation/randomization/gene_lists"

# Postproc union gene sets (pre-compiled, one gene per line)
POSTPROC_TOP="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
POSTPROC_BOTTOM="${SOURCE_BASE}/postproc/disambiguation_characterization/us_gs_relations/exports/txt/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   MODULE TOGGLES  —  set true / false for every module you want to run
# ══════════════════════════════════════════════════════════════════════════════

# CT discovery + resample + bootstrap (needs ALI_DIR, TRAIT_FILE, TREE_FILE)
RUN_CT=false

# Signification (needs DISCOVERY_INPUT, RESAMPLE_INPUT, BOOTSTRAP_INPUT)
RUN_SIGNIFICATION=true

# Disambiguation — set the sub-mode(s) to run
RUN_DISAMBIGUATION=true
DISAMBIG_ASR_MODE="precomputed"      # "compute" | "precomputed"
# If DISAMBIG_ASR_MODE="compute" the ASR is run fresh and written to ASR_CACHE_DIR
# If "precomputed" the cache is read from ASR_CACHE_DIR (must already exist)

# Post-processing
RUN_POSTPROC=true
POSTPROC_MODE="filter"               # "filter" | "exploratory"

# ORA enrichment from postproc gene lists
RUN_ORA=true

# STRING PPI enrichment (runs alongside ORA; requires RUN_ORA=true)
RUN_STRING=true

# CT accumulation (needs FILTERED_DISCOVERY, CLEANED_BACKGROUND, ALI_DIR)
RUN_ACCUMULATION=true

# ORA on accumulation gene lists
RUN_ACCUMULATION_ORA=true

# FADE positive selection (gene_set mode uses source-run gene sets)
RUN_FADE=true
FADE_MODE="gene_set"                 # "gene_set" | "all"

# MoleRate molecular rate shift detection
RUN_MOLERATE=true
MOLERATE_MODE="gene_set"             # "gene_set" | "all"

# RERConverge (needs GENE_TREES)
RUN_RER=false
RER_TOOL="build_trait,build_tree,build_matrix,continuous"
RER_GENE_SET_MODE="gene_set"         # "gene_set" | "all"

# Rmarkdown phenotype / dataset reports
RUN_REPORTING=true

# Foreground / background contrast selection
RUN_CONTRAST_SELECTION=false

# ── Execution mode ────────────────────────────────────────────────────────────
# "atomic"   – each enabled module runs as a separate nextflow invocation
#              (independent output dirs, easy to re-run one module at a time)
# "pipeline" – all enabled modules combined into a single nextflow run
#              (uses Nextflow channel chaining; mirrors a full-pipeline run)
RUN_MODE="pipeline"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   MODULE SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

# Nextflow Tower (set to true to stream run to tower.nf.io)
USE_TOWER=true

# Postproc posterior threshold for disambiguation filtering
CT_DISAMBIG_POSTERIOR_THRESHOLD=0.7

# ORA background: override if you want a different background file for ORA
ORA_BACKGROUND="${CLEANED_BACKGROUND}"

# Accumulation randomization N (may differ from CT bootstrap N)
ACCUMULATION_N_RANDOMIZATIONS="${N_RANDOMIZATIONS}"

# Whether to prune species from the trait file before analysis
PRUNE_DATA=false

# Reporting (Rmarkdown): enable/disable
REPORTING=true

# Contrast selection: enable/disable in integrated-style runs
CONTRAST_SELECTION=false

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#   OUTPUT SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

# Base directory for all output runs produced by this script.
# Each module run creates its own sub-directory: <label>/<timestamp>/
OUTPUT_BASE="${EXT_DRIVE}/CAAS_2.0/Results/${TRAIT}_from_source"

# Work directory — keeps Nextflow intermediate files.
# Uses a local filesystem path to avoid Apptainer nosuid/nodev issues.
WORK_BASE="${REPO_DIR}/work/${TRAIT}_from_source"

# Remove each run's work/ dir after completion (saves disk space)
CLEAN_WORK=true

# ─────────────────────────────────────────────────────────────────────────────
# Apptainer / Singularity cache/tmp on external drive
# ─────────────────────────────────────────────────────────────────────────────
export APPTAINER_TMPDIR="${EXT_DRIVE}/apptainer_tmp"
export APPTAINER_CACHEDIR="${EXT_DRIVE}/apptainer_cache"
export SINGULARITY_TMPDIR="${APPTAINER_TMPDIR}"
export SINGULARITY_CACHEDIR="${APPTAINER_CACHEDIR}"
mkdir -p "${APPTAINER_TMPDIR}" "${APPTAINER_CACHEDIR}"

# ─────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
#  (end of user-editable section)
# ══════════════════════════════════════════════════════════════════════════════
# ─────────────────────────────────────────────────────────────────────────────

timestamp=$(date +%Y%m%d_%H%M%S)

# ─────────────────────────────────────────────────────────────────────────────
# TOY-MODE RESOLUTION  — override paths/counts when IS_TOY=true
# ─────────────────────────────────────────────────────────────────────────────
if [ "$IS_TOY" = true ]; then
    ALI_DIR="$ALI_DIR_TOY"
    CYCLES="$CYCLES_TOY"
    N_RANDOMIZATIONS="$N_RANDOMIZATIONS_TOY"
    ACCUMULATION_N_RANDOMIZATIONS="$N_RANDOMIZATIONS_TOY"
fi

# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

section() {
    echo ""
    echo "══════════════════════════════════════════════════════════════"
    echo "  $*"
    echo "══════════════════════════════════════════════════════════════"
}

step()  { echo ""; echo "── $* ──"; }
ok()    { echo "  ✓ $*"; }
warn()  { echo "  ⚠ $*"; }
fail()  { echo "  ✗ ERROR: $*" >&2; exit 1; }

require_file() {
    local f="$1"
    [ -e "$f" ] || fail "Required input not found: $f"
}

maybe_clean_work() {
    local workdir="$1"
    if [ "$CLEAN_WORK" = true ] && [ -d "$workdir" ]; then
        rm -rf "$workdir"
        ok "Cleaned work dir: $workdir"
    fi
}

check_low_contrasts() {
    local outdir="$1" label="$2"
    local sentinel="${outdir}/low_contrasts.skip"
    if [ -f "$sentinel" ]; then
        warn "${label}: Fewer than 3 foreground contrasts — pipeline was skipped."
        warn "  Detail: $(cat "$sentinel")"
        return 0
    fi
    return 1
}

# Build the -with-tower flag (or nothing) based on USE_TOWER
tower_flag() {
    [ "$USE_TOWER" = true ] && echo "-with-tower" || echo ""
}

# Build --prune_data / --prune_list flags
prune_flags() {
    local flags=()
    if [ "$PRUNE_DATA" = true ] && [ -n "$PRUNE_LIST" ]; then
        flags+=(--prune_data --prune_list "$PRUNE_LIST")
        [ -n "$PRUNE_LIST_SECONDARY" ] && flags+=(--prune_list_secondary "$PRUNE_LIST_SECONDARY")
    else
        flags+=(--prune_data false)
    fi
    echo "${flags[@]+"${flags[@]}"}"
}

# Build optional trait flags (secondary, n, c, branch) — only when non-empty
optional_trait_flags() {
    local flags=()
    [ -n "$SECONDARY_TRAIT" ] && flags+=(--secondary_trait "$SECONDARY_TRAIT")
    [ -n "$N_TRAIT"          ] && flags+=(--n_trait "$N_TRAIT")
    [ -n "$C_TRAIT"          ] && flags+=(--c_trait "$C_TRAIT")
    [ -n "$BRANCH_TRAIT"     ] && flags+=(--branch_trait "$BRANCH_TRAIT")
    echo "${flags[@]+"${flags[@]}"}"
}

# Shared "everything off" disabled flags — added to every run that doesn't
# use a given module, so Nextflow params are explicit and unambiguous.
all_off_flags() {
    echo "--ct_signification false --ct_disambiguation false --ct_postproc false \
--ora false --string false --ct_accumulation false --reporting false \
--contrast_selection false"
}

# ─────────────────────────────────────────────────────────────────────────────
# PRE-FLIGHT VALIDATION
# ─────────────────────────────────────────────────────────────────────────────
validate_inputs() {
    section "VALIDATING INPUTS"

    # Always required
    require_file "$TRAIT_FILE"
    require_file "$TREE_FILE"

    # Required by module
    if [ "$RUN_CT" = true ]; then
        require_file "$ALI_DIR"
        require_file "$CAAS_CONFIG"
    fi

    if [ "$RUN_SIGNIFICATION" = true ]; then
        require_file "$DISCOVERY_INPUT"
        require_file "$RESAMPLE_INPUT"
        require_file "$BOOTSTRAP_INPUT"
        require_file "$BACKGROUND_INPUT"
    fi

    if [ "$RUN_DISAMBIGUATION" = true ]; then
        require_file "$GLOBAL_META_CAAS"
        require_file "$CAAS_CONFIG"
        require_file "$ASR_CACHE_DIR"
    fi

    if [ "$RUN_POSTPROC" = true ]; then
        require_file "$MASTER_CSV"
        require_file "$BACKGROUND_INPUT"
    fi

    if [ "$RUN_ORA" = true ]; then
        require_file "$ORA_GENE_LISTS_DIR"
        require_file "$CLEANED_BACKGROUND"
    fi

    if [ "$RUN_ACCUMULATION" = true ]; then
        require_file "$FILTERED_DISCOVERY"
        require_file "$CLEANED_BACKGROUND"
        require_file "$ALI_DIR"
        require_file "$CAAS_CONFIG"
    fi

    if [ "$RUN_ACCUMULATION_ORA" = true ]; then
        require_file "$ACCUMULATION_GENE_LISTS_DIR"
        require_file "$CLEANED_BACKGROUND"
    fi

    if [ "$RUN_FADE" = true ] || [ "$RUN_MOLERATE" = true ] || [ "$RUN_RER" = true ]; then
        require_file "$ALI_DIR"
        require_file "$CAAS_CONFIG"
        require_file "$POSTPROC_TOP"
        require_file "$POSTPROC_BOTTOM"
        require_file "$ACCUMULATION_CSV"
    fi

    if [ "$RUN_RER" = true ]; then
        require_file "$GENE_TREES"
    fi

    if [ "$RUN_STRING" = true ] && [ "$RUN_ORA" != true ]; then
        fail "RUN_STRING=true requires RUN_ORA=true (shared channel)"
    fi

    ok "All required inputs found."
}

# ─────────────────────────────────────────────────────────────────────────────
# PRINT CONFIGURATION SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print_config() {
    section "PHYLOPHERE MODULE RUNNER — CONFIGURATION"
    echo ""
    echo "Timestamp          : $timestamp"
    echo "Source run root    : $SOURCE_RUN_ROOT"
    echo "Source base        : $SOURCE_BASE"
    echo "Trait              : $TRAIT"
    echo "Trait file         : $TRAIT_FILE"
    echo "Tree               : $TREE_FILE"
    echo "Alignment dir      : $ALI_DIR"
    echo "Output base        : $OUTPUT_BASE"
    echo "Work base          : $WORK_BASE"
    echo "Use Tower          : $USE_TOWER"
    echo "Clean work dirs    : $CLEAN_WORK"
    echo ""
    echo "── Module toggles ──────────────────────────"
    echo "  RUN_CT                  : $RUN_CT"
    echo "  RUN_SIGNIFICATION       : $RUN_SIGNIFICATION"
    echo "  RUN_DISAMBIGUATION      : $RUN_DISAMBIGUATION  (asr_mode=${DISAMBIG_ASR_MODE})"
    echo "  RUN_POSTPROC            : $RUN_POSTPROC  (mode=${POSTPROC_MODE})"
    echo "  RUN_ORA                 : $RUN_ORA"
    echo "  RUN_STRING              : $RUN_STRING"
    echo "  RUN_ACCUMULATION        : $RUN_ACCUMULATION"
    echo "  RUN_ACCUMULATION_ORA    : $RUN_ACCUMULATION_ORA"
    echo "  RUN_FADE                : $RUN_FADE  (mode=${FADE_MODE})"
    echo "  RUN_MOLERATE            : $RUN_MOLERATE  (mode=${MOLERATE_MODE})"
    echo "  RUN_RER                 : $RUN_RER  (tool=${RER_TOOL}, mode=${RER_GENE_SET_MODE})"
    echo "  RUN_REPORTING           : $RUN_REPORTING"
    echo "  RUN_CONTRAST_SELECTION  : $RUN_CONTRAST_SELECTION"
    echo ""
    echo "── Execution ───────────────────────────────"
    echo "  RUN_MODE  : $RUN_MODE"
    echo "  IS_TOY    : $IS_TOY  (ali=${ALI_DIR}, cycles=${CYCLES}, rand=${N_RANDOMIZATIONS})"
    echo ""
    echo "── Optional traits ─────────────────────────"
    echo "  SECONDARY_TRAIT : ${SECONDARY_TRAIT:-(not set)}"
    echo "  N_TRAIT         : ${N_TRAIT:-(not set)}"
    echo "  C_TRAIT         : ${C_TRAIT:-(not set)}"
    echo "  BRANCH_TRAIT    : ${BRANCH_TRAIT:-(not set)}"
    echo ""
    echo "── Data pruning ────────────────────────────"
    echo "  PRUNE_DATA      : $PRUNE_DATA"
    echo "  PRUNE_LIST      : ${PRUNE_LIST:-(not set)}"
    echo "  PRUNE_LIST_SECONDARY : ${PRUNE_LIST_SECONDARY:-(not set)}"
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: CT discovery + resample + bootstrap
# ─────────────────────────────────────────────────────────────────────────────
run_ct() {
    local label="ct"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    require_file "$ALI_DIR"
    require_file "$CAAS_CONFIG"
    require_file "$BOOT_CONFIG"

    echo "Running CT (discovery + resample + bootstrap)"
    echo "  Alignment : $ALI_DIR"
    echo "  Output    : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_tool "discovery,resample,bootstrap" \
        --alignment "$ALI_DIR" \
        --my_traits "$TRAIT_FILE" \
        --traitname "$TRAIT" \
        $(optional_trait_flags) \
        --caas_config "$CAAS_CONFIG" \
        --traitvalues "$BOOT_CONFIG" \
        --tree "$TREE_FILE" \
        --cycles "$CYCLES" \
        $(prune_flags) \
        --ct_postproc false \
        --ct_signification false \
        --ct_disambiguation false \
        --ora false \
        --string false \
        --ct_accumulation false \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "CT completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: Signification
# ─────────────────────────────────────────────────────────────────────────────
run_signification() {
    local label="signification"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running Signification"
    echo "  discovery : $DISCOVERY_INPUT"
    echo "  bootstrap : $BOOTSTRAP_INPUT"
    echo "  background: $BACKGROUND_INPUT"
    echo "  Output    : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification \
        --discovery_input  "$DISCOVERY_INPUT" \
        --bootstrap_input  "$BOOTSTRAP_INPUT" \
        --background_input "$BACKGROUND_INPUT" \
        --my_traits  "$TRAIT_FILE" \
        --traitname  "$TRAIT" \
        $(optional_trait_flags) \
        --tree       "$TREE_FILE" \
        --caas_config "$CAAS_CONFIG" \
        --ct_postproc false \
        --ct_disambiguation false \
        --ora false \
        --string false \
        --ct_accumulation false \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Signification completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: Disambiguation
# ─────────────────────────────────────────────────────────────────────────────
run_disambiguation() {
    local label="disambiguation_${DISAMBIG_ASR_MODE}"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running Disambiguation (asr_mode=${DISAMBIG_ASR_MODE})"
    echo "  meta_caas : $GLOBAL_META_CAAS"
    echo "  asr cache : $ASR_CACHE_DIR"
    echo "  Output    : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification false \
        --ct_disambiguation \
        --ct_disambig_caas_metadata  "$GLOBAL_META_CAAS" \
        --caas_config                "$CAAS_CONFIG" \
        --tree                       "$TREE_FILE" \
        --my_traits                  "$TRAIT_FILE" \
        --traitname                  "$TRAIT" \
        $(optional_trait_flags) \
        --ct_disambig_asr_mode       "$DISAMBIG_ASR_MODE" \
        --ct_disambig_asr_cache_dir  "$ASR_CACHE_DIR" \
        --ct_disambig_posterior_threshold "$CT_DISAMBIG_POSTERIOR_THRESHOLD" \
        --ct_postproc false \
        --ora false \
        --string false \
        --ct_accumulation false \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Disambiguation (${DISAMBIG_ASR_MODE}) completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: Post-processing
# ─────────────────────────────────────────────────────────────────────────────
run_postproc() {
    local label="postproc_${POSTPROC_MODE}"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running PostProc (mode=${POSTPROC_MODE})"
    echo "  master CSV : $MASTER_CSV"
    echo "  background : $BACKGROUND_INPUT"
    echo "  Output     : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification false \
        --ct_disambiguation false \
        --ct_postproc \
        --caas_postproc_mode "$POSTPROC_MODE" \
        --discovery_input    "$MASTER_CSV" \
        --background_input   "$BACKGROUND_INPUT" \
        --my_traits          "$TRAIT_FILE" \
        --traitname          "$TRAIT" \
        $(optional_trait_flags) \
        --tree               "$TREE_FILE" \
        --caas_config        "$CAAS_CONFIG" \
        --ora false \
        --string false \
        --ct_accumulation false \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "PostProc (${POSTPROC_MODE}) completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: ORA enrichment (+ optional STRING)
# ─────────────────────────────────────────────────────────────────────────────
run_ora() {
    local label="ora"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    local string_flag="--string false"
    [ "$RUN_STRING" = true ] && string_flag="--string"

    echo "Running ORA enrichment"
    echo "  gene lists : $ORA_GENE_LISTS_DIR"
    echo "  background : $ORA_BACKGROUND"
    echo "  STRING     : $RUN_STRING"
    echo "  Output     : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification false \
        --ct_disambiguation false \
        --ct_postproc false \
        --ora \
        $string_flag \
        --ora_gene_lists_input "$ORA_GENE_LISTS_DIR" \
        --ora_background_input "$ORA_BACKGROUND" \
        --my_traits   "$TRAIT_FILE" \
        --traitname   "$TRAIT" \
        $(optional_trait_flags) \
        --tree        "$TREE_FILE" \
        --caas_config "$CAAS_CONFIG" \
        --ct_accumulation false \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "ORA completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: CT Accumulation
# ─────────────────────────────────────────────────────────────────────────────
run_accumulation() {
    local label="accumulation"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running CT Accumulation"
    echo "  filtered_discovery : $FILTERED_DISCOVERY"
    echo "  background         : $CLEANED_BACKGROUND"
    echo "  alignment          : $ALI_DIR"
    echo "  Output             : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification false \
        --ct_disambiguation false \
        --ct_postproc false \
        --ora false \
        --string false \
        --ct_accumulation \
        --accumulation_caas_input       "$FILTERED_DISCOVERY" \
        --accumulation_background_input "$CLEANED_BACKGROUND" \
        --accumulation_n_randomizations "$ACCUMULATION_N_RANDOMIZATIONS" \
        --alignment   "$ALI_DIR" \
        --caas_config "$CAAS_CONFIG" \
        --my_traits   "$TRAIT_FILE" \
        --traitname   "$TRAIT" \
        $(optional_trait_flags) \
        --tree        "$TREE_FILE" \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "CT Accumulation completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: ORA on Accumulation gene lists (+ optional STRING)
# ─────────────────────────────────────────────────────────────────────────────
run_accumulation_ora() {
    local label="accumulation_ora"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    local string_flag="--string false"
    [ "$RUN_STRING" = true ] && string_flag="--string"

    echo "Running ORA on Accumulation gene lists"
    echo "  gene lists : $ACCUMULATION_GENE_LISTS_DIR"
    echo "  background : $CLEANED_BACKGROUND"
    echo "  STRING     : $RUN_STRING"
    echo "  Output     : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --ct_signification false \
        --ct_disambiguation false \
        --ct_postproc false \
        --ct_accumulation false \
        --ora \
        $string_flag \
        --accumulation_caas_input             "$FILTERED_DISCOVERY" \
        --accumulation_background_input       "$CLEANED_BACKGROUND" \
        --accumulation_ora_gene_lists_input   "$ACCUMULATION_GENE_LISTS_DIR" \
        --ora_gene_lists_input                "$ACCUMULATION_GENE_LISTS_DIR" \
        --ora_background_input                "$CLEANED_BACKGROUND" \
        --alignment   "$ALI_DIR" \
        --caas_config "$CAAS_CONFIG" \
        --my_traits   "$TRAIT_FILE" \
        --traitname   "$TRAIT" \
        $(optional_trait_flags) \
        --tree        "$TREE_FILE" \
        --reporting false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "ORA Accumulation completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: FADE + MoleRate + RERConverge (selection suite)
# All three share the same gene-set inputs so they run in a single NF invocation
# when all are active, or individually when only some are on.
# ─────────────────────────────────────────────────────────────────────────────
run_selection() {
    local label="selection"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running Selection suite"
    echo "  FADE      : $RUN_FADE  (mode=${FADE_MODE})"
    echo "  MoleRate  : $RUN_MOLERATE  (mode=${MOLERATE_MODE})"
    echo "  RER       : $RUN_RER  (tool=${RER_TOOL})"
    echo "  postproc_top    : $POSTPROC_TOP"
    echo "  postproc_bottom : $POSTPROC_BOTTOM"
    echo "  accumulation    : $ACCUMULATION_CSV"
    echo "  Output    : $outdir"

    local nf_flags=(
        $(tower_flag)
        -profile local
        --my_traits  "$TRAIT_FILE"
        --traitname  "$TRAIT"
        --tree       "$TREE_FILE"
        --alignment  "$ALI_DIR"
        --caas_config  "$CAAS_CONFIG"
        --traitvalues  "$BOOT_CONFIG"
        --discovery_input  "$DISCOVERY_INPUT"
        --background_input "$BACKGROUND_INPUT"
        --reporting        false
        --contrast_selection false
    )

    # Disable all CT / ORA / accumulation unless explicitly requested here
    nf_flags+=(--ct_signification false --ct_disambiguation false --ct_postproc false)
    nf_flags+=(--ct_accumulation false)

    # Optional traits
    read -ra _opt_traits <<<"$(optional_trait_flags)" 2>/dev/null || true
    nf_flags+=("${_opt_traits[@]+"${_opt_traits[@]}"}")

    # Accumulation + ORA (pass-through so downstream annotation still works)
    if [ "$RUN_ACCUMULATION" = true ]; then
        nf_flags+=(
            --ct_accumulation
            --accumulation_caas_input       "$FILTERED_DISCOVERY"
            --accumulation_background_input "$CLEANED_BACKGROUND"
        )
    fi

    if [ "$RUN_ORA" = true ]; then
        nf_flags+=(--ora --ora_background_input "$ORA_BACKGROUND")
    else
        nf_flags+=(--ora false)
    fi
    if [ "$RUN_STRING" = true ]; then
        nf_flags+=(--string)
    else
        nf_flags+=(--string false)
    fi

    # FADE
    if [ "$RUN_FADE" = true ]; then
        nf_flags+=(
            --fade
            --fade_mode                 "$FADE_MODE"
            --fade_postproc_top         "$POSTPROC_TOP"
            --fade_postproc_bottom      "$POSTPROC_BOTTOM"
            --fade_accumulation_top     "$ACCUMULATION_CSV"
            --fade_accumulation_bottom  "$ACCUMULATION_CSV"
        )
    fi

    # MoleRate
    if [ "$RUN_MOLERATE" = true ]; then
        nf_flags+=(
            --molerate
            --molerate_mode                "$MOLERATE_MODE"
            --molerate_postproc_top        "$POSTPROC_TOP"
            --molerate_postproc_bottom     "$POSTPROC_BOTTOM"
            --molerate_accumulation_top    "$ACCUMULATION_CSV"
            --molerate_accumulation_bottom "$ACCUMULATION_CSV"
        )
    fi

    # RERConverge
    if [ "$RUN_RER" = true ]; then
        nf_flags+=(
            --rer_tool               "$RER_TOOL"
            --rer_gene_set_mode      "$RER_GENE_SET_MODE"
            --gene_trees             "$GENE_TREES"
            --rer_postproc_top       "$POSTPROC_TOP"
            --rer_postproc_bottom    "$POSTPROC_BOTTOM"
            --rer_accumulation_top   "$ACCUMULATION_CSV"
            --rer_accumulation_bottom "$ACCUMULATION_CSV"
        )
    fi

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        "${nf_flags[@]}" \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Selection suite completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# PIPELINE MODE — single nextflow invocation with all enabled modules
# ─────────────────────────────────────────────────────────────────────────────
run_pipeline() {
    local label="pipeline"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running PIPELINE mode (single Nextflow invocation)"
    echo "  Output : $outdir"

    local nf_flags=(
        -profile local
        --my_traits   "$TRAIT_FILE"
        --traitname   "$TRAIT"
        --tree        "$TREE_FILE"
        --alignment   "$ALI_DIR"
        --caas_config "$CAAS_CONFIG"
        --traitvalues "$BOOT_CONFIG"
        # Fallback discovery/background for modules that need them
        --discovery_input  "$DISCOVERY_INPUT"
        --background_input "$BACKGROUND_INPUT"
    )

    # Tower
    [ "$USE_TOWER" = true ] && nf_flags+=(-with-tower)

    # Optional traits
    [ -n "$SECONDARY_TRAIT" ] && nf_flags+=(--secondary_trait "$SECONDARY_TRAIT")
    [ -n "$N_TRAIT"          ] && nf_flags+=(--n_trait "$N_TRAIT")
    [ -n "$C_TRAIT"          ] && nf_flags+=(--c_trait "$C_TRAIT")
    [ -n "$BRANCH_TRAIT"     ] && nf_flags+=(--branch_trait "$BRANCH_TRAIT")

    # Prune
    if [ "$PRUNE_DATA" = true ] && [ -n "$PRUNE_LIST" ]; then
        nf_flags+=(--prune_data --prune_list "$PRUNE_LIST")
        [ -n "$PRUNE_LIST_SECONDARY" ] && nf_flags+=(--prune_list_secondary "$PRUNE_LIST_SECONDARY")
    else
        nf_flags+=(--prune_data false)
    fi

    # CT
    if [ "$RUN_CT" = true ]; then
        nf_flags+=(--ct_tool "discovery,resample,bootstrap" --cycles "$CYCLES")
    fi

    # Signification
    if [ "$RUN_SIGNIFICATION" = true ]; then
        nf_flags+=(
            --ct_signification
            --discovery_input  "$DISCOVERY_INPUT"
            --bootstrap_input  "$BOOTSTRAP_INPUT"
            --background_input "$BACKGROUND_INPUT"
        )
    else
        nf_flags+=(--ct_signification false)
    fi

    # Disambiguation
    if [ "$RUN_DISAMBIGUATION" = true ]; then
        nf_flags+=(
            --ct_disambiguation
            --ct_disambig_caas_metadata       "$GLOBAL_META_CAAS"
            --ct_disambig_asr_mode            "$DISAMBIG_ASR_MODE"
            --ct_disambig_asr_cache_dir       "$ASR_CACHE_DIR"
            --ct_disambig_posterior_threshold "$CT_DISAMBIG_POSTERIOR_THRESHOLD"
        )
    else
        nf_flags+=(--ct_disambiguation false)
    fi

    # PostProc
    if [ "$RUN_POSTPROC" = true ]; then
        nf_flags+=(
            --ct_postproc
            --caas_postproc_mode "$POSTPROC_MODE"
        )
    else
        nf_flags+=(--ct_postproc false)
    fi

    # ORA
    if [ "$RUN_ORA" = true ]; then
        nf_flags+=(
            --ora
            --ora_gene_lists_input "$ORA_GENE_LISTS_DIR"
            --ora_background_input "$ORA_BACKGROUND"
        )
    else
        nf_flags+=(--ora false)
    fi

    # STRING
    if [ "$RUN_STRING" = true ]; then
        nf_flags+=(--string)
    else
        nf_flags+=(--string false)
    fi

    # CT Accumulation
    if [ "$RUN_ACCUMULATION" = true ]; then
        nf_flags+=(
            --ct_accumulation
            --accumulation_caas_input       "$FILTERED_DISCOVERY"
            --accumulation_background_input "$CLEANED_BACKGROUND"
            --accumulation_n_randomizations "$ACCUMULATION_N_RANDOMIZATIONS"
        )
    else
        nf_flags+=(--ct_accumulation false)
    fi

    # Accumulation ORA gene-list pass-through
    if [ "$RUN_ACCUMULATION_ORA" = true ]; then
        nf_flags+=(--accumulation_ora_gene_lists_input "$ACCUMULATION_GENE_LISTS_DIR")
    fi

    # FADE
    if [ "$RUN_FADE" = true ]; then
        nf_flags+=(
            --fade
            --fade_mode                "$FADE_MODE"
            --fade_postproc_top        "$POSTPROC_TOP"
            --fade_postproc_bottom     "$POSTPROC_BOTTOM"
            --fade_accumulation_top    "$ACCUMULATION_CSV"
            --fade_accumulation_bottom "$ACCUMULATION_CSV"
        )
    fi

    # MoleRate
    if [ "$RUN_MOLERATE" = true ]; then
        nf_flags+=(
            --molerate
            --molerate_mode                "$MOLERATE_MODE"
            --molerate_postproc_top        "$POSTPROC_TOP"
            --molerate_postproc_bottom     "$POSTPROC_BOTTOM"
            --molerate_accumulation_top    "$ACCUMULATION_CSV"
            --molerate_accumulation_bottom "$ACCUMULATION_CSV"
        )
    fi

    # RERConverge
    if [ "$RUN_RER" = true ]; then
        nf_flags+=(
            --rer_tool                "$RER_TOOL"
            --rer_gene_set_mode       "$RER_GENE_SET_MODE"
            --gene_trees              "$GENE_TREES"
            --rer_postproc_top        "$POSTPROC_TOP"
            --rer_postproc_bottom     "$POSTPROC_BOTTOM"
            --rer_accumulation_top    "$ACCUMULATION_CSV"
            --rer_accumulation_bottom "$ACCUMULATION_CSV"
        )
    fi

    # Reporting
    if [ "$RUN_REPORTING" = true ]; then
        nf_flags+=(--reporting)
    else
        nf_flags+=(--reporting false)
    fi

    # Contrast selection
    if [ "$RUN_CONTRAST_SELECTION" = true ]; then
        nf_flags+=(--contrast_selection)
    else
        nf_flags+=(--contrast_selection false)
    fi

    nextflow run "${REPO_DIR}/main.nf" \
        "${nf_flags[@]}" \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Pipeline run completed → $outdir"
    check_low_contrasts "$outdir" "pipeline" || true
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: Reporting
# ─────────────────────────────────────────────────────────────────────────────
run_reporting() {
    local label="reporting"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running Reporting"
    echo "  traits : $TRAIT_FILE"
    echo "  tree   : $TREE_FILE"
    echo "  Output : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --reporting \
        --my_traits "$TRAIT_FILE" \
        --traitname "$TRAIT" \
        $(optional_trait_flags) \
        --tree "$TREE_FILE" \
        --ct_signification false \
        --ct_postproc false \
        --ct_disambiguation false \
        --ora false \
        --string false \
        --ct_accumulation false \
        --contrast_selection false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Reporting completed → $outdir"
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: Contrast Selection
# ─────────────────────────────────────────────────────────────────────────────
run_contrast_selection() {
    local label="contrast_selection"
    local outdir="${OUTPUT_BASE}/${label}/${timestamp}"
    local workdir="${WORK_BASE}/${label}/${timestamp}"
    mkdir -p "$outdir" "$workdir"

    echo "Running Contrast Selection"
    echo "  traits : $TRAIT_FILE"
    echo "  tree   : $TREE_FILE"
    echo "  Output : $outdir"

    # shellcheck disable=SC2046
    nextflow run "${REPO_DIR}/main.nf" \
        $(tower_flag) \
        -profile local \
        --contrast_selection \
        --my_traits "$TRAIT_FILE" \
        --traitname "$TRAIT" \
        $(optional_trait_flags) \
        --tree "$TREE_FILE" \
        --ct_signification false \
        --ct_postproc false \
        --ct_disambiguation false \
        --ora false \
        --string false \
        --ct_accumulation false \
        --reporting false \
        -w "$workdir" \
        --outdir "$outdir" \
        -with-report "${outdir}/pipeline_info/execution_report.html" \
        -with-trace  "${outdir}/pipeline_info/execution_trace.txt"

    ok "Contrast Selection completed → $outdir"
    check_low_contrasts "$outdir" "contrast_selection" || true
    maybe_clean_work "$workdir"
}

# ─────────────────────────────────────────────────────────────────────────────
# RESULT SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print_summary() {
    section "RUN COMPLETE"
    echo ""
    echo "Results root: ${OUTPUT_BASE}/"
    echo ""
    echo "Run directories created:"
    find "${OUTPUT_BASE}" -maxdepth 2 -mindepth 2 -type d 2>/dev/null | sort | while read -r d; do
        echo "  $d"
    done
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
main() {
    print_config
    validate_inputs

    # ── Pipeline mode: single NF invocation ─────────────────────────────────
    if [ "$RUN_MODE" = "pipeline" ]; then
        section "PIPELINE MODE — all enabled modules in one Nextflow run"
        run_pipeline
        print_summary
        return
    fi

    # ── Atomic mode: one NF invocation per module ────────────────────────────

    # Check whether any selection module is active
    local any_selection=false
    { [ "$RUN_FADE" = true ] || [ "$RUN_MOLERATE" = true ] || [ "$RUN_RER" = true ]; } \
        && any_selection=true

    # ── Individual modules ───────────────────────────────────────────────────
    if [ "$RUN_CT" = true ]; then
        section "CT — discovery + resample + bootstrap"
        run_ct
    fi

    if [ "$RUN_SIGNIFICATION" = true ]; then
        section "Signification"
        run_signification
    fi

    if [ "$RUN_DISAMBIGUATION" = true ]; then
        section "Disambiguation (${DISAMBIG_ASR_MODE})"
        run_disambiguation
    fi

    if [ "$RUN_POSTPROC" = true ]; then
        section "PostProc (${POSTPROC_MODE})"
        run_postproc
    fi

    if [ "$RUN_ORA" = true ] && [ "$any_selection" != true ]; then
        section "ORA enrichment"
        run_ora
    fi

    if [ "$RUN_ACCUMULATION" = true ] && [ "$any_selection" != true ]; then
        section "CT Accumulation"
        run_accumulation
    fi

    if [ "$RUN_ACCUMULATION_ORA" = true ]; then
        section "ORA on Accumulation gene lists"
        run_accumulation_ora
    fi

    # Selection suite (FADE + MoleRate + RER) runs together in one NF call
    # along with optional accumulation/ORA pass-through
    if [ "$any_selection" = true ]; then
        section "Selection suite (FADE=${RUN_FADE}, MoleRate=${RUN_MOLERATE}, RER=${RUN_RER})"
        run_selection
    fi

    if [ "$RUN_REPORTING" = true ]; then
        section "Reporting"
        run_reporting
    fi

    if [ "$RUN_CONTRAST_SELECTION" = true ]; then
        section "Contrast Selection"
        run_contrast_selection
    fi

    print_summary
}

main "$@"
