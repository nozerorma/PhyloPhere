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
# PHYLOPHERE: Multi-Phenotype Runner
#
# ─── OVERVIEW ─────────────────────────────────────────────────────────────────
# Orchestrates PhyloPhere runs for two independent biological projects, each
# requiring a different analytical configuration:
#
#   CLASS 1 · PRUNED-SECONDARY  (RUN_PRUNED_SECONDARY)
#   ─────────────────────────────────────────────────────
#   Comparative Oncology in Primates — IBE-UPF / Comparative Oncology Alliance
#
#   This class covers the cancer phenotypes of the Primates & Cancer project,
#   an ongoing collaboration between the IBE-UPF group and the Comparative
#   Oncology Alliance, whose aim is to uncover the genomic determinants of
#   differential cancer susceptibility across the primate clade. Phenotypic
#   data were obtained from Species360 / ZIMS veterinary records spanning
#   captive primate populations worldwide.
#
#   Two complementary tumour-burden phenotypes are analysed as primary traits
#   (each becoming the other's secondary trait):
#
#     • neoplasia_prevalence  – prevalence of all neoplastic disease
#                               (benign + malignant combined)
#     • malignant_prevalence  – prevalence restricted to malignant tumours
#
#   Both phenotypes are frequency/prevalence metrics derived from necropsy
#   records. Consequently CONTRAST SELECTION runs with full confidence-
#   interval (CI) support, enabled by two auxiliary population-size traits:
#
#     • n_trait  = adult_necropsy_count  – total adult necropsies per species
#                                          (population denominator for CI)
#     • c_trait  = neoplasia_necropsy /   – observed case count per species
#                  malignant_count      (swapped to match the primary trait)
#
#   A branch-colouring trait (LQ, Longevity Quotient) is carried through for
#   downstream visualisation.
#
#   Species pruning is mandatory for this class: per-trait exclusion lists
#   (neoplasia_exclude.txt / malignant_exclude.txt) remove species whose
#   necropsy sample sizes fall below quality thresholds, ensuring that only
#   species with statistically meaningful prevalence estimates enter the
#   contrast groups. Secondary prune lists mirror the primary lists with
#   swapped roles.
#
#   Full downstream chain: contrast_selection → CT (discovery + resample +
#   bootstrap) → signification → disambiguation (precomputed ASR) →
#   postproc / filter → ORA → STRING → CT accumulation →
#   FADE (directional AA selection) / MoleRate (rel. evo. rate) in gene_set
#   mode; RERConverge in gene_set mode (pre-filtered gene trees).
#
#
#   CLASS 2 · SIMPLE  (RUN_SIMPLE)
#   ────────────────────────────────
#   Primate Diet & Ethanol — MSc Thesis of María Sánchez Bermúdez
#   IBE-UPF / UOC · Supervised by Dr. David de Juan Sopeña &
#   Dr. José Francisco Sánchez Herrero
#
#   This class continues the comparative-genomics work initiated in the
#   Master's thesis of María Sánchez Bermúdez (MU Bioinformàtica i
#   Bioestadística, UOC, defended January 2026), carried out in close
#   collaboration with the IBE-UPF group. The project studies the genomic
#   mechanisms associated with the diverse dietary strategies observed across
#   159 primate species (16 families, 504 species in the order), including
#   specialised consumption of meat, invertebrates (insects), leaves, fruit,
#   seeds, and fermented fruit (ethanol). Raw dietary percentages and trophic-
#   guild classifications were curated from the literature and integrated with
#   the 233-species primate genomic catalogue of Kuderna et al. (2023).
#
#   María defined and derived the following phenotypes / indices, which are
#   the traits run through PhyloPhere in this class:
#
#     Raw diet fractions (% of diet, from EltonTraits / primary literature):
#       Diet.Inv    – invertebrate (insect) consumption
#       Diet.Vert   – vertebrate (meat) consumption
#       Diet.Fruit  – fruit consumption
#       Diet.Nect   – nectar consumption
#       Diet.Seed   – seed consumption
#       Diet.PlantO – other plant material
#
#     Trophic guild binary flags:
#       Herb  – herbivore (0/1)
#       Carn  – carnivore (0/1)
#
#     Trophic specialisation indices (derived by María; higher = more
#     specialised towards that guild relative to the species mean):
#       frug_idx     – frugivory index
#       fol_idx      – folivory index
#       ins_idx      – insectivory index
#       herb_idx     – herbivory index
#       omn_idx      – omnivory index
#       omn_spec_idx – omnivory specialisation index
#       folfrug_idx  – folivore-frugivore composite index
#
#     Dietary diversity (entropy-based metrics):
#       Shannon_norm – Shannon diversity of diet composition (normalised)
#       Shannon_spec – Shannon diversity at species level
#
#     Ethanol exposure:
#       Ethanol     – level of ethanol consumption / exposure (scored from
#                     field and captive observations; 0–3 scale)
#
#   IMPORTANT: for this project NO population-size metrics (n_trait / c_trait)
#   are available, as the dietary data are observational estimates rather than
#   individual-level clinical records. Accordingly, contrast_selection runs
#   WITHOUT CI support, and no species pruning is applied. Traits are run
#   serially, one per PhyloPhere invocation.
#
#   Full downstream chain: contrast_selection (no CI) → CT (discovery +
#   resample + bootstrap) → signification → disambiguation (precomputed ASR)
#   → postproc / filter → ORA → STRING → CT accumulation →
#   FADE / MoleRate (gene_set mode); RERConverge (gene_set mode, pre-filtered
#   gene trees).
#
#
# In both classes: --caas_postproc_mode filter; --ct_disambig_asr_mode
# precomputed; toy/full mode controlled by IS_TOY (cycles & randomisations).
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: run_phenotypes.sh
#

set -Eeuo pipefail

# Repo root — this script lives in the PhyloPhere directory
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

# Redirect Apptainer's tmp and cache directories to the external drive.
# /home is ~95% full (13 GB free) and /tmp is a 7.5 GB tmpfs — both too small
# to unpack the full conda environment from the Docker image layers (~5-8 GB
# uncompressed).  The external drive has ~65 GB free and is adequate for builds.
_EXT_DRIVE="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92"
export APPTAINER_TMPDIR="${_EXT_DRIVE}/apptainer_tmp"
export APPTAINER_CACHEDIR="${_EXT_DRIVE}/apptainer_cache"
export SINGULARITY_TMPDIR="${APPTAINER_TMPDIR}"
export SINGULARITY_CACHEDIR="${APPTAINER_CACHEDIR}"
mkdir -p "${APPTAINER_TMPDIR}" "${APPTAINER_CACHEDIR}"

# Pre-pull the Apptainer SIF if absent (same guard as run_neoplasia script).
_REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_APPTAINER_SIF="${_REPO_DIR}/apptainer/miralnso-phylophere-latest.img"
if [ ! -f "${_APPTAINER_SIF}" ]; then
    echo "[INFO] Pre-pulling Apptainer image to ${_APPTAINER_SIF} ..."
    mkdir -p "${_REPO_DIR}/apptainer"
    apptainer pull "${_APPTAINER_SIF}" docker://miralnso/phylophere:latest
    echo "[INFO] Apptainer image ready."
fi
unset _EXT_DRIVE _REPO_DIR _APPTAINER_SIF

timestamp=$(date +%Y%m%d_%H%M%S)

# ============================================================
# TOP-LEVEL TOGGLES
# ============================================================

# Which phenotype class(es) to run
RUN_PRUNED_SECONDARY=true   # neoplasia/malignant pair (pruning + secondary trait + CI)
RUN_SIMPLE=true            # user-defined list of unpruned, secondary-free phenotypes

# ── Selection analysis toggles ───────────────────────────────────────────────
# Both default to false; enable independently of each other.
# Mode 'all'       : runs on the full alignment directory (no upstream dependency).
# Mode 'gene_set'  : runs only on genes from a prior accumulation / postproc run;
#                    you must set the *_accumulation_* / *_postproc_* path params
#                    (see FADE_NF_FLAGS / MOLERATE_NF_FLAGS below for details).
# NOTE: when run inline (same Nextflow invocation as CT_ACCUMULATION), use
#       mode 'all' or point the gene-set paths at outputs of a *previous* run.
RUN_FADE=false              # HyPhy FADE (directional amino-acid selection)
RUN_MOLERATE=false          # HyPhy MoleRate (relative evolutionary rate)
RUN_RER=false               # RERconverge (trait-gene evolutionary-rate association)
FADE_MODE="gene_set"             # 'all' or 'gene_set'
MOLERATE_MODE="gene_set"         # 'all' or 'gene_set'
RER_MODE="gene_set"              # 'all' or 'gene_set'

# ============================================================
# TOY / FULL MODE
# Toy mode: small alignment subset, fewer cycles and randomizations.
# ============================================================
IS_TOY=false

DATADIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data"

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

CAAS_OUTBASE="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/CAAS_2.0/Results"
ASR_CACHE_DIR="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/asr"

# ============================================================
# SHARED INPUTS (used by both phenotype classes)
# ============================================================
TRAIT_FILE="${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"
TREE_FILE="${DATADIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
PRUNE_DIR="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP"

# BRANCH_TRAIT (LQ — Longevity Quotient)
#   LQ is carried through CLASS 1 runs in two complementary roles:
#
#   1. Phenotypic observation / visualisation (branch_trait role):
#      In the phenotype exploration report (2.Phenotype_exploration.Rmd),
#      LQ is mapped onto the phylogenetic tree via ancestral state
#      reconstruction, producing coloured-branch visualisations that reveal
#      how longevity evolves alongside cancer prevalence across the clade.
#      This provides a first, qualitative readout of the LQ–cancer
#      relationship in the phenotypic space.
#
#   2. Direct correlation with cancer phenotypes:
#      The exploration report also renders scatter plots and correlation
#      statistics between the primary/secondary cancer traits and LQ,
#      giving a quantitative estimate of how strongly longevity tracks
#      neoplasia or malignancy prevalence across species. These correlations
#      are computed both raw (Pearson/Spearman) and phylogenetically
#      corrected (PGLS), serving as a quality-control and biological
#      interpretation layer before the CAAS analysis itself.
#
#   LQ is NOT passed as secondary_trait to avoid confounding the CAAS
#   contrast-group definition (secondary_trait already taken by the
#   cross-cancer phenotype). Its informational content is captured at
#   the reporting/exploration stage instead.
BRANCH_TRAIT="LQ"

# ============================================================
# PRUNED-SECONDARY PHENOTYPE DEFINITIONS  (CLASS 1 — Cancer)
#
# Format per entry (colon-separated):
#   PRIMARY_TRAIT : SECONDARY_TRAIT : C_TRAIT : PRUNE_LIST : PRUNE_SECONDARY_LIST
#
# neoplasia/malignant are each run as primary with the other as secondary.
# c_trait is swapped accordingly (neoplasia_necropsy vs malignant_count)
# so the CI calculation always uses the case count matching the primary.
# n_trait (adult_necropsy_count) is shared: it is the common denominator
# (total adult necropsies) for both prevalence estimates.
# Prune lists mirror each other: a species excluded from the neoplasia
# analysis may still be retained for malignancy and vice versa.
# ============================================================
declare -a PRUNED_SECONDARY_DEFS=(
    "neoplasia_prevalence:malignant_prevalence:neoplasia_necropsy:neoplasia_exclude.txt:malignant_exclude.txt"
    "malignant_prevalence:neoplasia_prevalence:malignant_count:malignant_exclude.txt:neoplasia_exclude.txt"
)
N_TRAIT_PRUNED="adult_necropsy_count"

# ============================================================
# SIMPLE PHENOTYPE DEFINITIONS  (CLASS 2 — María Sánchez Bermúdez)
#
# Trait CSV produced during María's MSc thesis: dietary composition
# data (raw fractions + derived indices) for 159 primate species.
# Source: EltonTraits, primary literature, curated by María Sánchez
# Bermúdez (UOC / IBE-UPF, 2026). Supervised by Dr. David de Juan
# Sopeña and Dr. José Francisco Sánchez Herrero.
#
# Available traits (columns in SIMPLE_TRAIT_FILE):
#   Raw fractions  : Diet.Inv Diet.Vert Diet.Fruit Diet.Nect Diet.Seed Diet.PlantO
#   Trophic flags  : Herb Carn
#   Trophic indices: frug_idx fol_idx ins_idx omn_idx omn_spec_idx herb_idx folfrug_idx
#   Diversity      : Shannon_norm Shannon_spec
#   Ethanol        : Ethanol
# No n_trait / c_trait available → CI disabled; no pruning.
# ============================================================
SIMPLE_TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/maria_caas/Datos_fenotipos/diet_traitfile_comma.csv"

declare -a SIMPLE_PHENOTYPES=(
    # ── Raw dietary fractions (% of diet) ──────────────────────
    # "Diet.Inv"     # Invertebrate (insect) consumption
    # "Diet.Vert"    # Vertebrate (meat) consumption
    # "Diet.Fruit"   # Fruit consumption
    # "Diet.Nect"    # Nectar consumption
    # "Diet.Seed"    # Seed consumption
    # "Diet.PlantO"  # Other plant material
    # ── Trophic binary flags ────────────────────────────────────
    # "Herb"         # Herbivore flag (binary: 0/1)
    # "Carn"         # Carnivore flag (binary: 0/1)
    # ── Trophic specialisation indices (María-derived) ─────────
    "frug_idx"     # Frugivory index
    "fol_idx"      # Folivory index
    "ins_idx"      # Insectivory index
    "omn_idx"      # Omnivory index
    "omn_spec_idx" # Omnivory specialisation index
    "herb_idx"     # Herbivory index
    "folfrug_idx"  # Folivore-frugivore composite index
    # ── Dietary diversity ──────────────────────────────────────
    # "Shannon_norm" # Normalised Shannon entropy of diet composition
    # "Shannon_spec" # Shannon diversity at species level
    # ── Ethanol exposure ───────────────────────────────────────
    "Ethanol"      # Ethanol consumption level (0–3 ordinal scale)
)

# ============================================================
# INPUT VALIDATION
# ============================================================
for f in "$TRAIT_FILE" "$TREE_FILE"; do
    if [ ! -f "$f" ]; then
        echo "Error: required file not found: $f"
        exit 1
    fi
done

if [ "$RUN_SIMPLE" = true ] && [ ! -f "$SIMPLE_TRAIT_FILE" ]; then
    echo "Error: SIMPLE_TRAIT_FILE not found: $SIMPLE_TRAIT_FILE"
    exit 1
fi

if [ ! -d "$ALI_DIR" ]; then
    echo "Error: alignment directory not found: $ALI_DIR"
    exit 1
fi

# ============================================================
# FADE / MOLERATE NEXTFLOW FLAGS
# Extend or override per trait by editing below.
# In gene_set mode set the *_accumulation_* / *_postproc_* paths to the
# output files from a completed CT_ACCUMULATION run, e.g.:
#   --fade_accumulation_top  "${CAAS_OUTBASE}/${TRAIT}${TAG}/<timestamp>/filter/accumulation/\
#     accumulation_convergent_caap_top_aggregated_results.csv"
# ============================================================
FADE_NF_FLAGS=()
if [ "$RUN_FADE" = true ]; then
    FADE_NF_FLAGS=(
        --fade
        --fade_mode "$FADE_MODE"
        # Uncomment and fill for gene_set mode:
        # --fade_accumulation_top    "/path/to/accumulation_convergent_caap_top_aggregated_results.csv"
        # --fade_accumulation_bottom "/path/to/accumulation_convergent_caap_bottom_aggregated_results.csv"
        # --fade_postproc_top        "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
        # --fade_postproc_bottom     "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"
    )
fi

MOLERATE_NF_FLAGS=()
if [ "$RUN_MOLERATE" = true ]; then
    MOLERATE_NF_FLAGS=(
        --molerate
        --molerate_mode "$MOLERATE_MODE"
        # Uncomment and fill for gene_set mode:
        # --molerate_accumulation_top    "/path/to/accumulation_convergent_caap_top_aggregated_results.csv"
        # --molerate_accumulation_bottom "/path/to/accumulation_convergent_caap_bottom_aggregated_results.csv"
        # --molerate_postproc_top        "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
        # --molerate_postproc_bottom     "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"
    )
fi

RER_NF_FLAGS=()
if [ "$RUN_RER" = true ]; then
    RER_NF_FLAGS=(
        --rer_tool
        --rer_mode "$RER_MODE"
        # Uncomment and fill for gene_set mode:
        # --rer_postproc_top    "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt"
        # --rer_postproc_bottom "/path/to/special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt"
    )
fi

# ============================================================
# COMMON NEXTFLOW INFRASTRUCTURE FLAGS (shared across all runs)
# ============================================================
COMMON_NF_FLAGS=(
    -with-tower
    -profile local
    --ct_tool "discovery,resample,bootstrap"
    --alignment  "$ALI_DIR"
    --tree        "$TREE_FILE"
    --cycles      "$CYCLES"
    --accumulation_n_randomizations "$N_RANDOMIZATIONS"
    --ct_disambig_asr_mode      "precomputed"
    --ct_disambig_asr_cache_dir "${ASR_CACHE_DIR}"
)

# ============================================================
# HELPER: run a single complete pipeline (filter mode)
# ============================================================
# Usage:
#   run_pipeline <outdir> <workdir> <flags...>
run_pipeline() {
    local outdir="$1"
    local workdir="$2"
    shift 2
    local extra=("$@")

    mkdir -p "$outdir" "$workdir"

    nextflow run main.nf \
        "${COMMON_NF_FLAGS[@]}" \
        --caas_postproc_mode "filter" \
        -w "$workdir" \
        --outdir "$outdir" \
        "${extra[@]}"
}

# ============================================================
# SUMMARY HEADER
# ============================================================
echo "=========================================="
echo " PHYLOPHERE MULTI-PHENOTYPE RUNNER"
echo "=========================================="
echo " Timestamp :  $timestamp"
echo " IS_TOY    :  $IS_TOY  (tag: '${TAG}')"
echo " Cycles    :  $CYCLES"
echo " N_Rand    :  $N_RANDOMIZATIONS"
echo " Alignment :  $ALI_DIR"
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
# Comparative Oncology Alliance · IBE-UPF
# ============================================================
if [ "$RUN_PRUNED_SECONDARY" = true ]; then

    echo "=========================================="
    echo " CLASS 1: CANCER — PRUNED-SECONDARY"
    echo " Comparative Oncology in Primates · IBE-UPF"
    echo "=========================================="

    for def in "${PRUNED_SECONDARY_DEFS[@]}"; do

        # Parse colon-separated definition
        IFS=':' read -r TRAIT SECONDARY_TRAIT C_TRAIT PRUNE_FILE PRUNE_SECONDARY_FILE <<< "$def"

        PRUNE_LIST="${PRUNE_DIR}/${PRUNE_FILE}"
        PRUNE_SECONDARY_LIST="${PRUNE_DIR}/${PRUNE_SECONDARY_FILE}"

        # Validate prune lists exist
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

        # ── Low-contrast gate ────────────────────────────────────────────────
        # If CHECK_MIN_CONTRASTS in the pipeline determined that the traitfile
        # has too few foreground pairs, it writes a sentinel and exits cleanly.
        # Detect that here and skip to the next phenotype instead of reporting
        # a spurious success.
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
# MU Bioinformàtica i Bioestadística · UOC / IBE-UPF · 2026
# Supervised by Dr. D. de Juan Sopeña & Dr. J.F. Sánchez Herrero
# ============================================================
if [ "$RUN_SIMPLE" = true ]; then

    echo "=========================================="
    echo " CLASS 2: DIET (María Sánchez Bermúdez MSc)"
    echo " Primate Dietary Phenotypes · IBE-UPF / UOC"
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

        # ── Low-contrast gate ────────────────────────────────────────────────
        # If CHECK_MIN_CONTRASTS in the pipeline determined that the traitfile
        # has too few foreground pairs, it writes a sentinel and exits cleanly.
        # Detect that here and skip to the next phenotype instead of reporting
        # a spurious success.
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
#  Trait CSV (cancer_traits_processed-LQ.csv)       Trait CSV (Datos_melted.csv)
#  Prune lists (per-trait exclude lists)            No pruning
#          │                                                 │
#          ▼                                                 ▼
#   ┌─────────────────┐                           ┌──────────────────┐
#   │ REPORTING        │                           │ REPORTING         │
#   │ (prune + report) │                           │ (no prune)        │
#   └────────┬────────┘                           └────────┬─────────┘
#            │                                             │
#            ▼                                             ▼
#   ┌─────────────────┐                           ┌──────────────────┐
#   │ CONTRAST_SELECT  │ ← n_trait, c_trait (CI)  │ CONTRAST_SELECT  │ (no CI)
#   │ primary + sec.   │   branch_trait (LQ)       │ traitname only   │
#   └────────┬────────┘                           └────────┬─────────┘
#            │                                             │
#            ▼                                             ▼
#   ┌─────────────────────────────────────────────────────────────────┐
#   │              CT TOOL  (discovery → resample → bootstrap)        │
#   │                  --cycles N  --alignment alignment/dir          │
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
#   └──────────┬──────────────────┬────────────────────────┬──────────┘
#              │                  │                         │
#              ▼                  ▼                         ▼
#   ┌───────────────┐   ┌─────────────────┐    ┌──────────────────────┐
#   │      ORA       │   │    STRING        │    │   CT ACCUMULATION    │
#   │  (WebGestalt)  │   │  (rbioapi)       │    │  (gene-level perm.)  │
#   └───────────────┘   └─────────────────┘    └──────────────────────┘
#
# ============================================================
# DONE
# ============================================================
echo "=========================================="
echo " MULTI-PHENOTYPE RUNNER FINISHED"
echo " Timestamp: $timestamp"
echo " Results:   $CAAS_OUTBASE"
echo "=========================================="
