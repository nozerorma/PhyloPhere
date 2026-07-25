#!/bin/bash
#SBATCH --job-name=phylophere_posenrich
#SBATCH --partition=high-cpu
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -e Slurm/slurm-posenrich-%A_%a.err
#SBATCH -o Slurm/slurm-posenrich-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-2%2
#
# One-off manual re-render job (NOT part of the normal Nextflow pipeline).
# Runs ONLY the POSENRICH step (posenrich_enrich.py + 16.Position_enrichment_report.Rmd)
# for both traits. Split out of SBATCH_manual_fcs_rerun.sh so it can run
# independently of / in parallel with the STRING and Comparison steps there.
#
# Reuses the already-built position GMTs (external-DB/alignment derived --
# untouched by any of the RER/FADE/scoring fixes, no need to rebuild them).
#
# Changes vs the previous posenrich run:
#   - 'global' ranking dropped (top/bottom-only, consistent with every other
#     report; also the new scatter below pairs top vs bottom directly, and
#     'global' has no partner axis in that pairing)
#   - top-fraction cutoffs expanded to 25/10/5/1% via posenrich_enrich.py's
#     new default (--char-fracs) -- was 10/5/1% only; 100%/50% considered and
#     dropped (100% is a no-op foreground = full background, 50% is too
#     coarse to be informative for a top-fraction cutoff)
#   - hard fold-enrichment gate removed from significance (p.adj alone
#     decides now; fold_enrichment stays as a displayed/plotted column)
#   - new "Top vs Bottom (log2 fold enrichment)" scatter: x = log2fold at the
#     top ranking, y = log2fold at the bottom ranking, one point per
#     (database, pathway, top-fraction), coloured by cutoff, all four cutoffs
#     combined onto one plot, significant points labelled description_percentage
#   - leading-edge gene/position lists now truncate to 5 with a <details>
#     toggle for the rest, instead of one long comma-separated string

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
WORK="${BASE}/_manual_reports_20260720"

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT"
echo " BASE          : $BASE"
echo "======================================================"

fail() { echo "[FATAL] $1"; exit 1; }

echo ""
echo "─── POSENRICH ────────────────────────────────────────────────"
WD="${WORK}/posenrich_v2"
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

mkdir -p "${BASE}/posenrich"
cp posenrich_characterization.tsv posenrich_leading_edge.tsv "${BASE}/posenrich/"
cp "16.Position_enrichment_report_${TRAIT}.html" "${BASE}/posenrich/"
cp "16.Position_enrichment_report_${TRAIT}.html" "${BASE}/html_reports/"
echo "POSENRICH done"

echo ""
echo "======================================================"
echo " ALL DONE: $TRAIT"
echo "======================================================"
