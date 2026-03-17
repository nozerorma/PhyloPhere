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
# PHYLOPHERE: Multi-Phenotype Test Runner (TOY MODE)
#
# ─── OVERVIEW ─────────────────────────────────────────────────────────────────
#
# Drop-in counterpart of run_phenotypes.sh for local integration testing.
# Runs in TOY MODE exclusively (IS_TOY=true), using the bundled Data/test/
# inputs so that the full analytical chain can be verified quickly without
# needing the large production datasets or the external drive.
#
# Covers all pipeline modes presently available, including modes not yet
# exercised in the production run scripts:
#
#   ✓ CLASS 1 · PRUNED-SECONDARY  (one cancer phenotype: neoplasia_prevalence)
#   ✓ CLASS 2 · SIMPLE            (one diet phenotype:   frug_idx)
#   ✓ FADE                        (directional AA selection — gene_set mode)
#   ✓ MOLERATE                    (relative evolutionary rate — gene_set mode)
#   ✓ RERconverge                  (RERconverge association — gene_set mode)
#
# ORA on excluded genes (ORA_EXCLUDED) is automatic when --ora + --ct_postproc
# are both enabled; no separate toggle is required.
#
# Results land in:
#   ${REPO_DIR}/test_results/<trait><_toy>/<timestamp>/filter/
#
# Work directories are kept inside test_results so they do not clutter the
# production output base.  Re-runs automatically resume via -resume unless
# CLEAN_WORK=true is set below.
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: run_phenotypes_test.sh
#

set -Eeuo pipefail

# ============================================================
# Repo root — this script lives at the top of PhyloPhere
# ============================================================
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# Environment activation
# ============================================================
echo "Activating micromamba environment: phylophere"
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
else
    echo "Warning: micromamba not found, attempting conda activation"
    if command -v conda &> /dev/null; then
        conda activate phylophere || echo "Warning: could not activate conda environment"
    else
        echo "Warning: neither micromamba nor conda found, proceeding without activation"
    fi
fi

export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true

# Re-use the Apptainer tmp/cache on the external drive when available;
# fall back gracefully to /tmp for CI environments.
_EXT_DRIVE="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92"
if [ -d "${_EXT_DRIVE}" ]; then
    export APPTAINER_TMPDIR="${_EXT_DRIVE}/apptainer_tmp"
    export APPTAINER_CACHEDIR="${_EXT_DRIVE}/apptainer_cache"
else
    export APPTAINER_TMPDIR="/tmp/apptainer_tmp"
    export APPTAINER_CACHEDIR="/tmp/apptainer_cache"
fi
export SINGULARITY_TMPDIR="${APPTAINER_TMPDIR}"
export SINGULARITY_CACHEDIR="${APPTAINER_CACHEDIR}"
mkdir -p "${APPTAINER_TMPDIR}" "${APPTAINER_CACHEDIR}"

# Pre-pull the Apptainer SIF if absent
_APPTAINER_SIF="${REPO_DIR}/apptainer/miralnso-phylophere-latest.img"
if [ ! -f "${_APPTAINER_SIF}" ]; then
    echo "[INFO] Pre-pulling Apptainer image to ${_APPTAINER_SIF} ..."
    mkdir -p "${REPO_DIR}/apptainer"
    apptainer pull "${_APPTAINER_SIF}" docker://miralnso/phylophere:latest
    echo "[INFO] Apptainer image ready."
fi
unset _EXT_DRIVE _APPTAINER_SIF

timestamp=$(date +%Y%m%d_%H%M%S)

# ============================================================
# TOP-LEVEL TOGGLES
# ============================================================

# Which phenotype class(es) to test
RUN_PRUNED_SECONDARY=true   # CLASS 1: neoplasia_prevalence (pruning + secondary + CI)
RUN_SIMPLE=true             # CLASS 2: frug_idx (no pruning, no secondary/n/c traits)

# ── Selection analysis toggles ───────────────────────────────────────────────
# Toy inputs for gene_set mode live in Data/test/selection/.
RUN_FADE=true               # HyPhy FADE (directional AA selection)
RUN_MOLERATE=true           # HyPhy MoleRate (relative evolutionary rate)
RUN_RER=false                # RERconverge (trait-gene-rate association)

FADE_MODE="gene_set"
MOLERATE_MODE="gene_set"
RER_MODE="gene_set"

# Set to true to wipe work/ directories before running (forces full re-run).
CLEAN_WORK=false

# ============================================================
# TOY MODE — always on in this test script
# ============================================================
IS_TOY=true
TAG="_toy"
CYCLES="100"
N_RANDOMIZATIONS="1000"

# ============================================================
# PATHS (all self-contained inside PhyloPhere; no external drive needed)
# ============================================================

TEST_DATA="${REPO_DIR}/Data"
ASR_CACHE_DIR="${REPO_DIR}/Data/asr"
ALI_DIR="${REPO_DIR}/Data/2.Alignments/Ali_toy"

# Gene Ensembl mapping file (gene ID → gene symbol + coordinates).
# The default in conf/common.config points to the cluster path; override for
# local testing with the copy present in NEOPLASY_PRIMATES.
GENE_ENSEMBL_FILE="${TEST_DATA}/2.Alignments/ensembl_genes.output"

# Output base: results land inside the repo under test_results/
CAAS_OUTBASE="${REPO_DIR}/test_results"

# ============================================================
# CLASS 1 — PRUNED-SECONDARY INPUTS
# ============================================================
TRAIT_FILE="${TEST_DATA}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"
TREE_FILE="${TEST_DATA}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
PRUNE_DIR="${TEST_DATA}/1.Cancer_data/Neoplasia_species360/prune_lists"
BRANCH_TRAIT="LQ"
N_TRAIT_PRUNED="adult_necropsy_count"

# Single phenotype for a fast test pass (add the second entry to exercise both)
declare -a PRUNED_SECONDARY_DEFS=(
    # "neoplasia_prevalence:malignant_prevalence:neoplasia_necropsy:neoplasia_exclude.txt:malignant_exclude.txt"
    # "malignant_prevalence:neoplasia_prevalence:malignant_count:malignant_exclude.txt:neoplasia_exclude.txt"
)

# ============================================================
# CLASS 2 — SIMPLE INPUTS  (María Sánchez Bermúdez — Diet)
# ============================================================
SIMPLE_TRAIT_FILE="${TEST_DATA}/11.Diet_phenotypes/diet_traitfile_comma.csv"

# Single trait for a fast test pass
declare -a SIMPLE_PHENOTYPES=(
#     "frug_idx"
    # "fol_idx"
    # "ins_idx"
    "Ethanol"
)

# ============================================================
# INPUT VALIDATION
# ============================================================
for f in "$TRAIT_FILE" "$TREE_FILE"; do
    if [ ! -f "$f" ]; then
        echo "Error: required test input file not found: $f"
        exit 1
    fi
done

if [ "$RUN_SIMPLE" = true ] && [ ! -f "$SIMPLE_TRAIT_FILE" ]; then
    echo "Error: SIMPLE_TRAIT_FILE not found: $SIMPLE_TRAIT_FILE"
    exit 1
fi

if [ ! -d "$ALI_DIR" ]; then
    echo "Error: toy alignment directory not found: $ALI_DIR"
    exit 1
fi

if [ ! -f "$GENE_ENSEMBL_FILE" ]; then
    echo "Error: gene_ensembl_file not found: $GENE_ENSEMBL_FILE"
    exit 1
fi

# ============================================================
# SELECTION MODEL PRE-COMPUTED INPUTS (for gene_set mode)
# Toy versions live in Data/test/selection/
# ============================================================
SEL_TEST_DIR="${REPO_DIR}/Data/test/selection"
SEL_ACC_TOP="${SEL_TEST_DIR}/accumulation_convergent_caap_top_aggregated_results.csv"
SEL_ACC_BOTTOM="${SEL_TEST_DIR}/accumulation_convergent_caap_bottom_aggregated_results.csv"
SEL_PP_TOP="${SEL_TEST_DIR}/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
SEL_PP_BOTTOM="${SEL_TEST_DIR}/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"

for sel_f in "$SEL_ACC_TOP" "$SEL_ACC_BOTTOM" "$SEL_PP_TOP" "$SEL_PP_BOTTOM"; do
    if [ ! -f "$sel_f" ]; then
        echo "Warning: selection test input not found, disabling gene_set modes: $sel_f"
        RUN_FADE=false
        RUN_MOLERATE=false
        RUN_RER=false
        break
    fi
done

# ============================================================
# FADE / MOLERATE NEXTFLOW FLAGS
# ============================================================
FADE_NF_FLAGS=()
if [ "$RUN_FADE" = true ]; then
    FADE_NF_FLAGS=(
        --fade
        --fade_mode "$FADE_MODE"
        --fade_accumulation_top    "$SEL_ACC_TOP"
        --fade_accumulation_bottom "$SEL_ACC_BOTTOM"
        --fade_postproc_top        "$SEL_PP_TOP"
        --fade_postproc_bottom     "$SEL_PP_BOTTOM"
    )
fi

MOLERATE_NF_FLAGS=()
if [ "$RUN_MOLERATE" = true ]; then
    MOLERATE_NF_FLAGS=(
        --molerate
        --molerate_mode "$MOLERATE_MODE"
        --molerate_accumulation_top    "$SEL_ACC_TOP"
        --molerate_accumulation_bottom "$SEL_ACC_BOTTOM"
        --molerate_postproc_top        "$SEL_PP_TOP"
        --molerate_postproc_bottom     "$SEL_PP_BOTTOM"
    )
fi

# ── RERconverge flags ─────────────────────────────────────────────────────────
# RER_MAIN automatically reads --my_traits (passed via COMMON or per-class) and
# uses sel_pp_*/sel_acc_* channels wired from CT_ACCUMULATION / CT_POSTPROC.
# In standalone gene_set mode the postproc paths are forwarded explicitly.
RER_NF_FLAGS=()
if [ "$RUN_RER" = true ]; then
    RER_NF_FLAGS=(
        --rer_tool
        --rer_mode "$RER_MODE"
        --rer_postproc_top    "$SEL_PP_TOP"
        --rer_postproc_bottom "$SEL_PP_BOTTOM"
    )
fi

# ============================================================
# COMMON NEXTFLOW INFRASTRUCTURE FLAGS
# ============================================================
COMMON_NF_FLAGS=(
    -profile local
    --ct_tool "discovery,resample,bootstrap"
    --alignment  "$ALI_DIR"
    --tree        "$TREE_FILE"
    --cycles      "$CYCLES"
    --accumulation_n_randomizations "$N_RANDOMIZATIONS"
    --ct_disambig_asr_mode      "precomputed"
    --ct_disambig_asr_cache_dir "${ASR_CACHE_DIR}"
    --gene_ensembl_file         "$GENE_ENSEMBL_FILE"
    # Taxon-ID mapping required by CT_DISAMBIGUATION tip-level analysis
    --tax_id                    "${REPO_DIR}/Data/5.Phylogeny/taxid_species_family.tsv"
    --publish_intermediates     true
    --export_groups             "debug"
    --export_perm_discovery     "debug"
    # Tower is intentionally disabled for local testing; re-enable with:
    # -with-tower
)

# ============================================================
# HELPER: run_pipeline <outdir> <workdir> [extra nextflow flags...]
# ============================================================
run_pipeline() {
    local outdir="$1"
    local workdir="$2"
    shift 2
    local extra=("$@")

    if [ "$CLEAN_WORK" = true ] && [ -d "$workdir" ]; then
        echo "  [clean] Removing work directory: $workdir"
        rm -rf "$workdir"
    fi

    mkdir -p "$outdir" "$workdir"

    nextflow run "${REPO_DIR}/main.nf" \
        "${COMMON_NF_FLAGS[@]}" \
        --caas_postproc_mode "filter" \
        -resume \
        -w "$workdir" \
        --outdir "$outdir" \
        "${extra[@]}"
}

# ============================================================
# SUMMARY HEADER
# ============================================================
echo "=========================================="
echo " PHYLOPHERE TEST RUNNER  (TOY MODE)"
echo "=========================================="
echo " Timestamp    :  $timestamp"
echo " IS_TOY       :  true  (tag: '${TAG}')"
echo " Cycles       :  $CYCLES"
echo " N_Rand       :  $N_RANDOMIZATIONS"
echo " Alignment    :  $ALI_DIR"
echo " Output base  :  $CAAS_OUTBASE"
echo " Clean work   :  $CLEAN_WORK"
echo "------------------------------------------"
echo " Classes to run:"
echo "   Pruned-secondary : $RUN_PRUNED_SECONDARY"
echo "   Simple           : $RUN_SIMPLE"
echo "------------------------------------------"
echo " Selection analyses:"
echo "   FADE       : $RUN_FADE     (mode: $FADE_MODE)"
echo "   MoleRate   : $RUN_MOLERATE  (mode: $MOLERATE_MODE)"
echo "   RERconverge: $RUN_RER      (mode: $RER_MODE)"
echo "=========================================="
echo ""

# ============================================================
# CLASS 1: PRUNED-SECONDARY PHENOTYPES  (Cancer in Primates)
# ============================================================
if [ "$RUN_PRUNED_SECONDARY" = true ]; then

    echo "=========================================="
    echo " CLASS 1: CANCER — PRUNED-SECONDARY (TOY)"
    echo "=========================================="

    for def in "${PRUNED_SECONDARY_DEFS[@]}"; do

        IFS=':' read -r TRAIT SECONDARY_TRAIT C_TRAIT PRUNE_FILE PRUNE_SECONDARY_FILE <<< "$def"

        PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
        PRUNE_SECONDARY_LIST="${PRUNE_DIR}/${PRUNE_SECONDARY_FILE}"

        for pf in "$PRUNE_LIST" "$PRUNE_SECONDARY_LIST"; do
            if [ ! -f "$pf" ]; then
                echo "Warning: prune list not found (skipping '$TRAIT'): $pf"
                continue 2
            fi
        done

        RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${timestamp}"
        WORK_DIR="${RESULTS_BASE}/work"
        mkdir -p "$RESULTS_BASE" "$WORK_DIR"

        echo ""
        echo "------------------------------------------"
        echo " Running: $TRAIT"
        echo "   Secondary     : $SECONDARY_TRAIT"
        echo "   n_trait       : $N_TRAIT_PRUNED"
        echo "   c_trait       : $C_TRAIT"
        echo "   Branch trait  : $BRANCH_TRAIT"
        echo "   Prune list    : $PRUNE_LIST"
        echo "   Prune (sec.)  : $PRUNE_SECONDARY_LIST"
        echo "   Output        : $RESULTS_BASE"
        echo "------------------------------------------"

        PRUNED_FLAGS=(
            --my_traits         "$TRAIT_FILE"
            --traitname         "$TRAIT"
            --secondary_trait   "$SECONDARY_TRAIT"
            --n_trait           "$N_TRAIT_PRUNED"
            --c_trait           "$C_TRAIT"
            --branch_trait      "$BRANCH_TRAIT"
            --prune_data
            --prune_list        "$PRUNE_LIST"
            --prune_list_secondary "$PRUNE_SECONDARY_LIST"
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
            "${PRUNED_FLAGS[@]}"

        if [ -f "${RESULTS_BASE}/filter/low_contrasts.skip" ]; then
            echo ""
            echo " ⚠ Skipped '${TRAIT}': fewer than ${MIN_CONTRASTS:-3} foreground contrasts."
            echo "   Detail: $(cat "${RESULTS_BASE}/filter/low_contrasts.skip")"
            echo ""
            continue
        fi

        echo ""
        echo " ✓ Completed: $TRAIT"
        echo ""

    done

    echo "=========================================="
    echo " CLASS 1 complete."
    echo "=========================================="
    echo ""
fi

# ============================================================
# CLASS 2: SIMPLE PHENOTYPES  (María Sánchez Bermúdez — Diet)
# ============================================================
if [ "$RUN_SIMPLE" = true ]; then

    echo "=========================================="
    echo " CLASS 2: DIET — SIMPLE (TOY)"
    echo "=========================================="

    if [ "${#SIMPLE_PHENOTYPES[@]}" -eq 0 ]; then
        echo "Warning: SIMPLE_PHENOTYPES list is empty – nothing to run."
    fi

    for TRAIT in "${SIMPLE_PHENOTYPES[@]}"; do

        RESULTS_BASE="${CAAS_OUTBASE}/${TRAIT}${TAG}/${timestamp}"
        WORK_DIR="${RESULTS_BASE}/work"
        mkdir -p "$RESULTS_BASE" "$WORK_DIR"

        echo ""
        echo "------------------------------------------"
        echo " Running: $TRAIT"
        echo "   Trait file    : $SIMPLE_TRAIT_FILE"
        echo "   (no pruning, no secondary/n/c traits)"
        echo "   Output        : $RESULTS_BASE"
        echo "------------------------------------------"

        SIMPLE_DISCRETE_METHOD="decile"
        if [ "$TRAIT" = "Ethanol" ]; then
            SIMPLE_DISCRETE_METHOD="quintile"
        fi

        SIMPLE_FLAGS=(
            --my_traits         "$SIMPLE_TRAIT_FILE"
            --traitname         "$TRAIT"
            --branch_trait      "$BRANCH_TRAIT"
            --discrete_method   "$SIMPLE_DISCRETE_METHOD"
            --n_trait           ""
            --c_trait           ""
            --secondary_trait   ""
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
            "${SIMPLE_FLAGS[@]}"

        if [ -f "${RESULTS_BASE}/filter/low_contrasts.skip" ]; then
            echo ""
            echo " ⚠ Skipped '${TRAIT}': fewer than ${MIN_CONTRASTS:-3} foreground contrasts."
            echo "   Detail: $(cat "${RESULTS_BASE}/filter/low_contrasts.skip")"
            echo ""
            continue
        fi

        echo ""
        echo " ✓ Completed: $TRAIT"
        echo ""

    done

    echo "=========================================="
    echo " CLASS 2 complete."
    echo "=========================================="
    echo ""
fi

# ============================================================
# ASCII WORKFLOW OVERVIEW
# ============================================================
#
#  CLASS 1 — Cancer / Pruned-Secondary             CLASS 2 — Diet / Simple
#  ─────────────────────────────────────            ───────────────────────
#  Test traits (cancer_traits_processed-LQ.csv)    Diet traits (frug_idx …)
#  Test prune lists (neoplasia/malignant_exclude)   No pruning
#          │                                                 │
#          ▼                                                 ▼
#   ┌─────────────────┐                           ┌──────────────────┐
#   │ REPORTING        │                           │ REPORTING         │
#   └────────┬────────┘                           └────────┬─────────┘
#            │                                             │
#            ▼                                             ▼
#   ┌─────────────────┐                           ┌──────────────────┐
#   │ CONTRAST_SELECT  │ ← n_trait, c_trait (CI)  │ CONTRAST_SELECT  │ (no CI)
#   └────────┬────────┘                           └────────┬─────────┘
#            │                                             │
#            └────────────────────┬────────────────────────┘
#                                 ▼
#   ┌─────────────────────────────────────────────────────────────────┐
#   │              CT TOOL  (discovery → resample → bootstrap)        │
#   │                  --cycles 100  --alignment Ali_toy              │
#   └─────────────────────────────┬───────────────────────────────────┘
#                                 │
#                                 ▼
#   ┌─────────────────────────────────────────────────────────────────┐
#   │              CT SIGNIFICATION  (bootstrap p-values)             │
#   └─────────────────────────────┬───────────────────────────────────┘
#                                 │
#                                 ▼
#   ┌─────────────────────────────────────────────────────────────────┐
#   │         CT DISAMBIGUATION  (ASR convergence, precomputed)       │
#   └─────────────────────────────┬───────────────────────────────────┘
#                                 │
#                                 ▼
#   ┌─────────────────────────────────────────────────────────────────┐
#   │              CT POSTPROC  (filter mode, gene filtering)         │
#   └──────┬──────────────┬──────────────────┬────────────┬───────────┘
#          │              │                  │            │
#          ▼              ▼                  ▼            ▼
#   ┌───────────┐  ┌───────────┐   ┌──────────────┐ ┌───────────────┐
#   │    ORA     │  │ ORA       │   │ STRING        │ │ CT ACCUMUL.   │
#   │ (sig genes)│  │ EXCLUDED  │   │ (rbioapi)     │ │ (gene perm.)  │
#   └───────────┘  │(excl.lists│   └──────────────┘ └──────┬────────┘
#                  │ bg_ori)   │                            │
#                  └───────────┘                            │
#                          ┌─────────────────────────────── ┘
#                          │
#                          ▼
#             ┌────────────────────────────────────┐
#             │ FADE  /  MOLERATE  /  RERconverge   │
#             │  (gene_set mode, Data/test/selection│
#             │   pre-computed inputs)              │
#             └────────────────────────────────────┘
#
# ============================================================
# DONE
# ============================================================
echo "=========================================="
echo " TEST RUNNER FINISHED"
echo " Timestamp : $timestamp"
echo " Results   : $CAAS_OUTBASE"
echo "=========================================="
