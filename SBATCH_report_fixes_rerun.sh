#!/bin/bash
#SBATCH --job-name=phylophere_plot_fixes
#SBATCH --partition=high-cpu
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -e Slurm/slurm-plotfix-%A_%a.err
#SBATCH -o Slurm/slurm-plotfix-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-2%2
#
# One-off manual re-render job (NOT part of the normal Nextflow pipeline).
# Copies SBATCH_posenrich_only.sh's pattern (isolated per-report work dir,
# reuse already-computed upstream inputs, only re-render the display code
# that changed) to cover the 2026-07-22 plot fixes across TWO reports:
#   1) posenrich_enrich.py + 16.Position_enrichment_report.Rmd (posenrich/)
#   2) 4.Independent_contrasts.Rmd -- "4b. Phylogenetic Tree of Selected
#      Contrast Pairs" (data_exploration/)
#
# Both fixes are display-only (no upstream recompute needed):
#   - both reports: geom_text_repel/geom_tiplab no longer duplicate the point
#     layer's legend key with an overlaid letter glyph (show.legend = FALSE)
#   - 4.Independent_contrasts.Rmd tree: tip labels offset off the tips
#     (0.03 x tree depth, was sitting directly on them), larger tip-label/
#     legend/title font sizes, 1.8x thicker branches
#
# Step 2 re-passes the SAME seed (1998) and the SAME already-pruned
# trait/tree files the original run produced, so contrast selection itself
# is reproduced identically -- only the plot code differs.

set -o pipefail
mkdir -p Slurm

source /etc/profile
source ~/.bashrc
conda deactivate 2>/dev/null
conda activate phylophere

REPO_DIR="/homes/users/mramon/scratch/0.Phylophere"

case $SLURM_ARRAY_TASK_ID in
    1) TRAIT="neoplasia_prevalence" ;;
    2) TRAIT="malignant_prevalence" ;;
esac

BASE="/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_true/${TRAIT}"
WORK="${BASE}/_manual_reports_20260722"

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT"
echo " BASE          : $BASE"
echo "======================================================"

fail() { echo "[FATAL] $1"; exit 1; }

# ── 1) POSENRICH: legend fix only, upstream Fisher computation unchanged ────
echo ""
echo "─── [1/2] POSENRICH ────────────────────────────────────────────"
WD="${WORK}/posenrich_v3"
rm -rf "$WD"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .

python3 "${REPO_DIR}/subworkflows/ENRICHMENT/local/src/posenrich_enrich.py" \
    --obs-scores      "${BASE}/scoring/position_scores.tsv" \
    --gmt-dir         "${BASE}/posenrich/gmts" \
    --characterization "${BASE}/posenrich/gmts/characterization_layers.tsv" \
    --universe        "${BASE}/postproc/cleaned_backgrounds/cleaned_background_main.txt" \
    --background      "${BASE}/caastools/background.output" \
    --annot-file      "${BASE}/scoring/fcs_stats.tsv" \
    --cosmic-coverage "${BASE}/posenrich/gmts/cosmic_coverage_genes.txt" \
    --pai3d-coverage  "${BASE}/posenrich/gmts/pai3d_coverage_genes.txt" \
    --min-size 5 \
    --max-size 0 \
    --char-fracs "0.25,0.10,0.05,0.01" \
    --padj-thr 0.15 \
    --fold-thr 1.5 \
    --output-dir . \
    || fail "posenrich_enrich.py failed"

[ -f posenrich_characterization.tsv ] || fail "posenrich_enrich.py produced no posenrich_characterization.tsv"

cat > render_posenrich.R <<RSCRIPT
rmarkdown::render(
    '16.Position_enrichment_report.Rmd',
    params = list(
        results_file       = 'posenrich_characterization.tsv',
        leading_edge_file  = 'posenrich_leading_edge.tsv',
        traitname           = '${TRAIT}',
        padj_thr            = 0.15,
        fold_thr            = 1.5
    ),
    output_file = '16.Position_enrichment_report_${TRAIT}.html'
)
RSCRIPT

Rscript render_posenrich.R || fail "POSENRICH report failed"
[ -f "16.Position_enrichment_report_${TRAIT}.html" ] || fail "POSENRICH report produced no html"

cp posenrich_characterization.tsv posenrich_leading_edge.tsv "${BASE}/posenrich/"
cp "16.Position_enrichment_report_${TRAIT}.html" "${BASE}/posenrich/"
cp "16.Position_enrichment_report_${TRAIT}.html" "${BASE}/html_reports/"
echo "[1/2] POSENRICH done"

# ── 2) Independent contrasts tree: re-render with the SAME seed/inputs so ───
# the contrast selection itself is bit-identical -- only the plot code changed.
echo ""
echo "─── [2/2] Independent contrasts tree ───────────────────────────"
WD="${WORK}/independent_contrasts_v1"
rm -rf "$WD"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/TRAIT_ANALYSIS/local/"* .

case "$TRAIT" in
    malignant_prevalence)
        TRAIT_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs_final/malignant_prevalence/64/018d741d99775bd5c5134ebee4f176/data_exploration/0.Data-pruning/pruned_trait_file.tsv"
        TREE_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs_final/malignant_prevalence/64/018d741d99775bd5c5134ebee4f176/data_exploration/0.Data-pruning/pruned_tree_file.nwk"
        C_TRAIT="malignant_count"
        SECONDARY_TRAIT="neoplasia_prevalence"
        ;;
    neoplasia_prevalence)
        TRAIT_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs_final/neoplasia_prevalence/30/3ce3eec44cae0513d9f6ee2f5ac383/data_exploration/0.Data-pruning/pruned_trait_file.tsv"
        TREE_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs_final/neoplasia_prevalence/30/3ce3eec44cae0513d9f6ee2f5ac383/data_exploration/0.Data-pruning/pruned_tree_file.nwk"
        C_TRAIT="neoplasia_necropsy"
        SECONDARY_TRAIT="malignant_prevalence"
        ;;
esac

[ -f "$TRAIT_FILE" ] || fail "pruned trait file not found: $TRAIT_FILE"
[ -f "$TREE_FILE" ] || fail "pruned tree file not found: $TREE_FILE"
cp "$TRAIT_FILE" pruned_trait_file.tsv
cp "$TREE_FILE" pruned_tree_file.nwk

# The Rmd isn't standalone: it reads upstream files already accumulated into
# data_exploration/ by reports 0-3 (e.g. trait_stats.csv from
# 1.Data-exploration/1.Species_distribution/) -- the real Nextflow process
# receives this whole tree as its staged `results_dir` input. Reproduce that
# here by seeding the work dir with the already-published copy.
[ -d "${BASE}/data_exploration" ] || fail "no existing data_exploration/ to seed from: ${BASE}/data_exploration"
mkdir -p data_exploration
cp -r "${BASE}/data_exploration"/. data_exploration/

cat > render_contrasts.R <<RSCRIPT
rmarkdown::render(
    '4.Independent_contrasts.Rmd',
    params = list(
        trait_file = 'pruned_trait_file.tsv',
        tree_file = 'pruned_tree_file.nwk',
        output_dir = 'data_exploration',
        seed = '1998',
        clade_name = 'primates',
        taxon_of_interest = 'family',
        traitname = '${TRAIT}',
        n_trait = 'adult_necropsy_count',
        c_trait = '${C_TRAIT}',
        tax_id = '/homes/users/mramon/scratch/2.Primates/1.Primates_data/5.Phylogeny/taxid_species_family_primates.tsv',
        secondary_trait = '${SECONDARY_TRAIT}',
        branch_trait = 'LQ',
        discrete_method = 'quintile',
        top_quantile = '0.90',
        bottom_quantile = '0.10',
        contrast_max_iter = '3'
    ),
    output_file = '4.Independent_contrasts_${TRAIT}.html',
    envir = new.env()
)
RSCRIPT

Rscript render_contrasts.R || fail "Independent contrasts report failed"
[ -f "4.Independent_contrasts_${TRAIT}.html" ] || fail "Independent contrasts report produced no html"

# Left in html_reports/ alongside the rest -- NOT overwriting the original
# data_exploration/ tree (this re-render's own copy lives in $WD/data_exploration
# for inspection, e.g. $WD/data_exploration/.../selected_pairs_tree.png).
cp "4.Independent_contrasts_${TRAIT}.html" "${BASE}/html_reports/4.Independent_contrasts_${TRAIT}_fixed.html"
echo "[2/2] Independent contrasts tree done"

echo ""
echo "======================================================"
echo " ALL DONE: $TRAIT"
echo "======================================================"
