#!/usr/bin/env bash
#
# Integrated 1000-gene test WITH the CAAS permulation-excess framework.
# Copy of run_integrated_1000.sh tuned for the perms feature (cycles=1000, full_perms=10).
#
# IMPORTANT — why this script does NOT use -resume by default:
#   Nextflow caches tasks by the process *command* (e.g. `ct discovery`), NOT by the
#   Python modules the `ct` binary loads at runtime (disco.py/boot.py/caas_id.py...).
#   After a change to those modules (e.g. the header canonicalization), -resume happily
#   serves STALE old-schema discovery → CT_SIGNIFICATION then fails
#   "Missing required columns: pvalue". So a clean run is required the first time after
#   any CAAStools-module edit. See docs/CAAS_PERMULATION_RUNTIME.md.
#
# Usage:
#   bash test/run_integrated_1000_perms.sh                 # clean run (safe default)
#   RESUME=1 bash test/run_integrated_1000_perms.sh        # -resume (only if modules unchanged)
#   CYCLES=200 FULL_PERMS=20 bash test/run_integrated_1000_perms.sh   # override knobs
#   ... any extra args are passed through to `nextflow run`.

set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="${REPO_DIR}/test"
INPUTS_DIR="${TEST_DIR}/inputs"
ALIGNMENT_DIR="${REPO_DIR}/test/inputs/alignments/Ali_1000"

# Tunable knobs (env-overridable)
CYCLES="${CYCLES:-1000}"
FULL_PERMS="${FULL_PERMS:-100}"
RESUME="${RESUME:-0}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="/media/miguel/phd/run/integrated_1000_perms/${TIMESTAMP}"
WORKDIR="${OUTDIR}/work"

TRAIT_FILE="${INPUTS_DIR}/traitfiles/traitfile.tab"
BOOT_TRAIT_FILE="${INPUTS_DIR}/traitfiles/boot_traitfile.tab"
TREE_FILE="${INPUTS_DIR}/tree/pruned_tree_file.nwk"
ASR_CACHE_DIR="${REPO_DIR}/test/inputs/asr"

require_path() {
    local path="$1"
    if [[ ! -e "$path" ]]; then
        echo "Error: required path not found: $path" >&2
        exit 1
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# ENVIRONMENT ACTIVATION
# ─────────────────────────────────────────────────────────────────────────────
echo "Activating environment: phylophere"
# Fallback so `micromamba shell hook` works under set -u even when .bashrc didn't set it.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
elif command -v conda &> /dev/null; then
    conda activate phylophere || echo "Warning: Could not activate conda environment"
else
    echo "Warning: Neither micromamba nor conda found — proceeding without activation"
fi

require_path "$ALIGNMENT_DIR"
require_path "$TRAIT_FILE"
require_path "$BOOT_TRAIT_FILE"
require_path "$TREE_FILE"
require_path "$ASR_CACHE_DIR"

mkdir -p "$OUTDIR" "$WORKDIR"

echo "=========================================="
echo "PHYLOPHERE 1000-GENE TEST — CAAS PERMULATION FRAMEWORK"
echo "=========================================="
echo "cycles            : $CYCLES"
echo "caas_full_perms   : $FULL_PERMS"
echo "resume            : $([ "$RESUME" = 1 ] && echo yes || echo 'no (clean run)')"
echo "Alignments        : $ALIGNMENT_DIR"
echo "Output dir        : $OUTDIR"
echo "Work dir          : $WORKDIR"
echo "=========================================="

NF_FLAGS=(
    -profile local
    --ali_format fasta
    --reporting false
    --ct_tool "discovery,resample,bootstrap"
    --ct_signification
    --ct_disambiguation
    --ct_postproc
    --ct_accumulation
    --traitname "neoplasia_prevalence"
    --scoring
    --cycles "$CYCLES"
    --caas_permulation_enrichment true
    --caas_full_perms "$FULL_PERMS"
    --gene_ensembl_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/ensembl_genes.output"
    --tax_id "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/taxid_species_family_primates.tsv"
    --ali_sp_names "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/ali_sp_names.txt"
    --vep_map_dir "${REPO_DIR}/test/inputs/map_1000/"
    --alignment "$ALIGNMENT_DIR"
    --caas_config "$TRAIT_FILE"
    --traitvalues "$BOOT_TRAIT_FILE"
    --tree "$TREE_FILE"
    --ct_disambig_asr_mode precomputed
    --ct_disambig_asr_cache_dir "$ASR_CACHE_DIR"
    --ct_disambig_posterior_threshold 0.1
    --caas_postproc_mode filter
    --gene_filter_mode both
    --generate_reports true
    --asr_robustness true
    --scoring_rer_input "${REPO_DIR}/test/inputs/rer/RERConverge/RER_Results/rerconverge_summary_neoplasia_prevalence.tsv"
    --vep
    --enrichment true
    --gmt_dir "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/GMTs"
    --posenrich true
    --domain_variability_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/PFAM/domain_variability.tsv"
    --ucr_positions_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/UCR/ucr_positions.tsv"
    --fubar_sites_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/FUBAR/fubar_sites.tsv"
    --vep_primateai_db "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/PAI3D/PrimateAI-3D.hg38.txt.gz"
    --egg_members_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/eggNOG/9443_members.tsv.gz"
    --egg_annotations_file "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/eggNOG/9443_annotations.tsv.gz"
    --cosmic_db "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/4.External_DBs/COSMIC/Cosmic_MutantCensus_v104_GRCh38.tsv.gz"
    --outdir "$OUTDIR"
)
# Opt-in resume only (default clean — see header note on the module-cache trap).
if [ "$RESUME" = 1 ]; then
    NF_FLAGS+=( -resume )
fi

nextflow run "${REPO_DIR}/main.nf" \
    -with-tower \
    "${NF_FLAGS[@]}" \
    -w "$WORKDIR" \
    "$@"

echo ""
echo "Permulation framework test finished."
echo "Results:          $OUTDIR"
echo "Null matrix:      $OUTDIR/caas_permulation/caas_perms.rds"
echo "CAAS FCS report:  $OUTDIR/fcs/FCS_scoring_neoplasia_prevalence.html  (p.perm on the *_asr rankings)"
# POSENRICH (position-wise XL-mHG): background = the in-run caastools background.output
# (discovery) restricted to cleaned_background — no manual background input needed.
echo "POSENRICH report: $OUTDIR/posenrich/POSENRICH_report_neoplasia_prevalence.html"
echo "POSENRICH results:$OUTDIR/posenrich/posenrich_results.tsv"
