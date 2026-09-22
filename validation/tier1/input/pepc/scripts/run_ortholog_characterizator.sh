#!/usr/bin/env bash
# run_ortholog_characterizator.sh — launcher for ortholog_characterizator
# (quality -> translation -> phylogeny -> positive_selection) against a
# Tier 1 fixture's codon-level CDS alignment (align_cds/*.fasta, built by
# each fixture's own build_cds.py).
#
# Shared by both fixtures currently wired up this way (PEPC) rather than
# hand-rolling a one-off invocation per fixture. Runs with quality/dedup off
# (see below) and stages .fa symlinks for the CDS dir, since
# ortholog_characterizator's own Channel.fromPath hardcodes that extension.
#
# Usage:
#   bash validation/tier1/run_ortholog_characterizator.sh pepc
#   bash validation/tier1/run_ortholog_characterizator.sh --cds-dir <dir> --species-tree <nwk> --out-dir <dir> --label <name>
#
# The first form is a shortcut for the two fixtures already set up this way;
# it resolves CDS_DIR/SPECIES_TREE/OUT_DIR relative to this file. The second
# form works for any future fixture built the same way (align_cds/*.fasta +
# a species tree.nwk), without editing this script.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OC_ROOT="/home/miguel/IBE-UPF/PhD/ortholog_characterizator"

# ── Resolve fixture shortcut or explicit flags ───────────────────────────────
CDS_DIR=""
SPECIES_TREE=""
OUT_DIR=""
LABEL=""
OUTGROUP_TAXA=""

case "${1:-}" in
  pepc)
    CDS_DIR="${HERE}/../align_cds"
    SPECIES_TREE="${HERE}/../tree.nwk"
    OUT_DIR="${HERE}/../oc_run"
    LABEL="pepc"
    shift
    ;;
esac

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cds-dir)      CDS_DIR="$2"; shift 2 ;;
    --species-tree) SPECIES_TREE="$2"; shift 2 ;;
    --out-dir)      OUT_DIR="$2"; shift 2 ;;
    --label)        LABEL="$2"; shift 2 ;;
    --outgroup)     OUTGROUP_TAXA="$2"; shift 2 ;;
    --run-quality)  RUN_QUALITY="$2"; shift 2 ;;
    --run-psel)     RUN_POSITIVE_SELECTION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -n "$CDS_DIR" && -d "$CDS_DIR" ]]      || { echo "ERROR: --cds-dir missing or not a directory: $CDS_DIR" >&2; exit 1; }
[[ -n "$SPECIES_TREE" && -f "$SPECIES_TREE" ]] || { echo "ERROR: --species-tree missing: $SPECIES_TREE" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: --out-dir required" >&2; exit 1; }
[[ -n "$LABEL" ]]   || LABEL="run"

RUN_TRANSLATION="${RUN_TRANSLATION:-true}"
# PHYLOGENY (GENE_TREES + ASTRAL SPECIES_TREE) is not wired to POSITIVE_SELECTION
# in main.nf -- psel_species_tree is always a separately-supplied external tree
# (the fixture's own tree.nwk here) -- so it's dead compute for this launcher.
RUN_PHYLOGENY="${RUN_PHYLOGENY:-false}"
RUN_VARIABILITY="${RUN_VARIABILITY:-false}"
RUN_POSITIVE_SELECTION="${RUN_POSITIVE_SELECTION:-true}"
FILTER_MODE="${FILTER_MODE:-none}"   # our CDS is already curated per-tip; no assembly-quality filtering to do
HYPHY_METHODS="${HYPHY_METHODS:-FUBAR}"   # FEL not needed downstream (FADE is the consumer, not FEL)
# MEME is triggered automatically by ortholog_characterizator as a second pass
# on any gene with >= meme_min_sites FUBAR-positive sites (conf/
# positive_selection.config default: 1) -- there's no plain on/off flag for
# it, so an unreachable threshold is how you turn it off. Off by default:
# MEME is far slower than FUBAR (process_medium's 48h*attempt time budget vs
# FUBAR's own, see resources.config) and PEPC's single-gene FUBAR-positive
# call already triggers it on the only gene there is.
MEME_MIN_SITES="${MEME_MIN_SITES:-999999}"

mkdir -p "$OUT_DIR"
META_DIR="${OUT_DIR}/oc_metadata"
mkdir -p "$META_DIR"

shopt -s nullglob
CDS_FILES=("${CDS_DIR}"/*.fa "${CDS_DIR}"/*.fasta)
shopt -u nullglob
[[ ${#CDS_FILES[@]} -gt 0 ]] || { echo "ERROR: no .fa/.fasta files in ${CDS_DIR}" >&2; exit 1; }

# ortholog_characterizator's own Channel.fromPath hardcodes "*.fa" (see
# main.nf); stage .fa symlinks alongside any .fasta-named build_cds.py
# output rather than rename the fixture's own committed file.
STAGED_CDS_DIR="${META_DIR}/cds_dir_fa"
mkdir -p "$STAGED_CDS_DIR"
rm -f "${STAGED_CDS_DIR}"/*.fa
for f in "${CDS_FILES[@]}"; do
  base="$(basename "$f")"
  ln -sf "$f" "${STAGED_CDS_DIR}/${base%.fasta}.fa"
done
N_TIPS=$(grep -ch "^>" "${CDS_FILES[@]}" | paste -sd+ - | bc 2>/dev/null || echo "?")

# The quality/dedup stage (assemblies_tsv/species_dedup_map/RENAME_ASSEMBLIES)
# is built for the mammal use case: multiple genome assembly versions per
# species, picked down to one by assembly quality. None of our fixtures have
# that -- every tip here is already a single, specific individual/accession --
# so quality is off by default and translation reads straight from the CDS
# dir via --dedup_dir (the same *.fa shape quality's own output would have
# produced, per make_fallback_channel('dedup_dir', ...) in main.nf).
RUN_QUALITY="${RUN_QUALITY:-false}"
CDS_DIR="$STAGED_CDS_DIR"

echo "════════════════════════════════════════════════════════"
echo " ortholog_characterizator — ${LABEL}"
echo " CDS dir      : ${CDS_DIR}  (${N_TIPS} tips)"
echo " Species tree : ${SPECIES_TREE}"
echo " Output       : ${OUT_DIR}"
echo " Stages       : quality=${RUN_QUALITY} translation=${RUN_TRANSLATION}"
echo "                phylogeny=${RUN_PHYLOGENY} positive_selection=${RUN_POSITIVE_SELECTION}"
echo "════════════════════════════════════════════════════════"

OPTIONAL_ARGS=()
[[ -n "$OUTGROUP_TAXA" ]] && OPTIONAL_ARGS+=(--outgroup_taxa "$OUTGROUP_TAXA")
# --cds_dir is read unconditionally at startup regardless of --run_quality;
# --dedup_dir is what translation actually reads from when quality is off.
[[ "$RUN_QUALITY" == "false" ]] && OPTIONAL_ARGS+=(--dedup_dir "$CDS_DIR")
OPTIONAL_ARGS+=(--meme_min_sites "$MEME_MIN_SITES")

# Run from ortholog_characterizator's own directory, not PhyloPhere's: Nextflow
# auto-loads a `nextflow.config` from the launch directory in addition to the
# pipeline's own, and PhyloPhere's repo-root config is full of process
# selectors (SCORING_COMPUTE, FADE_RUN, ...) that don't exist in this
# pipeline at all -- harmless "no process matching" noise at best, a real
# resource/executor clash at worst. An explicit -w also keeps this run's
# work dir under the fixture's own output tree instead of dropping a stray
# work/ into PhyloPhere's repo root.
cd "$OC_ROOT"

RUN_NAME="${LABEL}_$(date +%Y%m%d_%H%M%S)"
time nextflow run "${OC_ROOT}/main.nf" \
    -name "$RUN_NAME" \
    -w "${OUT_DIR}/work" \
    -with-report   "${OUT_DIR}/nextflow_report_${RUN_NAME}.html" \
    -with-trace    "${OUT_DIR}/nextflow_trace_${RUN_NAME}.tsv" \
    -with-timeline "${OUT_DIR}/nextflow_timeline_${RUN_NAME}.html" \
    --cds_dir                "$CDS_DIR" \
    --out_dir                "$OUT_DIR" \
    --run_quality            "$RUN_QUALITY" \
    --run_translation        "$RUN_TRANSLATION" \
    --run_phylogeny          "$RUN_PHYLOGENY" \
    --run_variability        "$RUN_VARIABILITY" \
    --run_positive_selection "$RUN_POSITIVE_SELECTION" \
    --filter_mode            "$FILTER_MODE" \
    --gene_tree_method       iqtree \
    --psel_species_tree      "$SPECIES_TREE" \
    --hyphy_methods          "$HYPHY_METHODS" \
    "${OPTIONAL_ARGS[@]}" \
    -profile local \
    -resume
