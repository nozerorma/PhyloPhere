#!/bin/bash
#SBATCH --job-name=phylophere_fcs_fix
#SBATCH --partition=high-cpu
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -e Slurm/slurm-fcsfix-%A_%a.err
#SBATCH -o Slurm/slurm-fcsfix-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-2%2
#
# One-off manual re-render job (NOT part of the normal Nextflow pipeline).
# Re-renders, for one trait, in order:
#   1) 6.FADE_report.Rmd for top and bottom (full gene_mode='all' JSON set)
#   2) scoring_compute.R (picks up the fresh fade_summary_{top,bottom}.tsv)
#   3) 13.FCS_general_report.Rmd for RER   (scoring/rer/)
#   4) 13.FCS_general_report.Rmd for FADE  (scoring/fade/)
#   5) 15.STRING_report.Rmd                (string/)
#   6) posenrich_enrich.py + 16.Position_enrichment_report.Rmd (posenrich/)
#   7) 12.Comparison_report.Rmd            (compare/) -- pulls in 3-6
#
# Background: scoring_compute.R's gene_scores join was fixed (left_join ->
# full_join, was silently truncating RER/FADE/accum genes to the CAAS gene
# set) and the FCS universe wiring was fixed (RER/FADE now use their own
# tested-gene background instead of cleaned_background_main.txt). FADE was
# also re-run in full 'all' mode (was gene_set-restricted before). STRING now
# always runs with PPI enrichment on and tabbed term tables/networks.
# POSENRICH's hard fold-enrichment gate was removed (p.adj alone decides
# significance now; fold_enrichment stays as a displayed column). The
# Comparison report is now per-direction throughout (was module-only/pooled
# before) and reuses STRING's + this run's fresh module outputs.

set -o pipefail
mkdir -p Slurm

source /etc/profile
source ~/.bashrc
conda deactivate 2>/dev/null
conda activate phylophere

REPO_DIR="/homes/users/mramon/scratch/0.Phylophere"
GMT_DIR="/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/GMTs"
STRING_DB_DIR="/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/STRING"
export FADE_REPORT_CORES="${SLURM_CPUS_PER_TASK:-8}"

# Per-step toggles -- override at submission to split the work across
# multiple concurrent jobs instead of one long sequential one. Steps 1-4 are
# already correct from a prior run and untouched by the latest report
# changes, so a typical split is:
#   J1=$(sbatch --export=ALL,RUN_FADE_REPORT=false,RUN_SCORING=false,RUN_RER_FCS=false,RUN_FADE_FCS=false,RUN_POSENRICH=false,RUN_COMPARE=false SBATCH_manual_fcs_rerun.sh | awk '{print $4}')
#   sbatch        --export=ALL,RUN_FADE_REPORT=false,RUN_SCORING=false,RUN_RER_FCS=false,RUN_FADE_FCS=false,RUN_STRING=false,RUN_COMPARE=false   SBATCH_manual_fcs_rerun.sh   # POSENRICH, runs in parallel with J1 -- independent of STRING
#   sbatch --dependency=afterok:$J1 --export=ALL,RUN_FADE_REPORT=false,RUN_SCORING=false,RUN_RER_FCS=false,RUN_FADE_FCS=false,RUN_STRING=false,RUN_POSENRICH=false SBATCH_manual_fcs_rerun.sh   # Comparison, needs STRING's fresh output -> waits on J1 only
# All toggles default to true so a plain `sbatch SBATCH_manual_fcs_rerun.sh`
# still runs the full 7-step chain in one job.
RUN_FADE_REPORT="${RUN_FADE_REPORT:-true}"
RUN_SCORING="${RUN_SCORING:-true}"
RUN_RER_FCS="${RUN_RER_FCS:-true}"
RUN_FADE_FCS="${RUN_FADE_FCS:-true}"
RUN_STRING="${RUN_STRING:-true}"
RUN_POSENRICH="${RUN_POSENRICH:-true}"
RUN_COMPARE="${RUN_COMPARE:-true}"

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

# ── 1) FADE individual reports (top + bottom, sequential) ───────────────────
# 6.FADE_report.Rmd now reduces each gene's JSON straight to one summary row
# instead of materialising a full (site x AA) table, so this no longer needs
# the huge memory the first attempt OOM'd on. Kept resilient regardless:
# only fail this step if fade_summary_*.tsv is missing (HTML failure is
# logged but not fatal), and skip re-rendering entirely if a prior attempt
# already produced the summary (idempotent re-run after a crash).
if [ "$RUN_FADE_REPORT" = "true" ]; then
for DIR in top bottom; do
    echo ""
    echo "─── [1/7] FADE report: $DIR ──────────────────────────────────"
    WD="${WORK}/fade_report_v2_${DIR}"
    mkdir -p "$WD"
    cd "$WD" || fail "cannot cd $WD"

    if [ -s "fade_summary_${DIR}.tsv" ]; then
        echo "[1/7] fade_summary_${DIR}.tsv already present from a prior attempt — skipping re-render"
    else
        cp -R "${REPO_DIR}/subworkflows/FADE/local/"* .

        cat > "render_fade_${DIR}.R" <<RSCRIPT
rmarkdown::render(
    '6.FADE_report.Rmd',
    params = list(
        json_dir       = '${BASE}/selection/fade/${DIR}',
        direction      = '${DIR}',
        traitname      = '${TRAIT}',
        bf_threshold   = 100,
        min_genes_hmap = 2,
        output_dir     = '${BASE}/selection/fade/${DIR}',
        fg_list_file   = '${BASE}/selection/species_sets/${DIR}_species.txt'
    ),
    output_file = '6.FADE_report_${DIR}.html'
)
RSCRIPT

        Rscript "render_fade_${DIR}.R"
        RC=$?
        if [ $RC -ne 0 ]; then
            echo "[WARN] FADE report ($DIR) Rscript exited $RC (possibly OOM in the HTML/heatmap tail) — checking whether the summary TSV still landed"
        fi
    fi

    [ -s "fade_summary_${DIR}.tsv" ] || fail "FADE report ($DIR): fade_summary_${DIR}.tsv missing — no usable output to continue with"

    mkdir -p "${BASE}/selection/fade/${DIR}"
    if [ -f "6.FADE_report_${DIR}.html" ]; then
        cp "6.FADE_report_${DIR}.html" "${BASE}/selection/fade/${DIR}/"
        cp "6.FADE_report_${DIR}.html" "${BASE}/html_reports/"
    else
        echo "[WARN] no 6.FADE_report_${DIR}.html produced — continuing without it (fcs_stats/comparison don't need it)"
    fi
    cp "fade_summary_${DIR}.tsv" "${BASE}/selection/fade/${DIR}/"
    echo "[1/7] FADE report ($DIR) done: $(wc -l < ${BASE}/selection/fade/${DIR}/fade_summary_${DIR}.tsv) rows in fade_summary_${DIR}.tsv"
done

# ── FADE's own universe: union of genes in the fresh fade_summary_{top,bottom}.tsv ──
mkdir -p "${WORK}/universes"
{ tail -n +2 "${BASE}/selection/fade/top/fade_summary_top.tsv" | cut -f1
  tail -n +2 "${BASE}/selection/fade/bottom/fade_summary_bottom.tsv" | cut -f1
} | sort -u > "${WORK}/universes/fade_background_${TRAIT}.txt"
echo "FADE universe: $(wc -l < ${WORK}/universes/fade_background_${TRAIT}.txt) genes"
else
    echo "[1/7] RUN_FADE_REPORT=false — skipping FADE report re-render (reusing existing selection/fade/*/fade_summary_*.tsv)"
    mkdir -p "${WORK}/universes"
    { tail -n +2 "${BASE}/selection/fade/top/fade_summary_top.tsv" | cut -f1
      tail -n +2 "${BASE}/selection/fade/bottom/fade_summary_bottom.tsv" | cut -f1
    } | sort -u > "${WORK}/universes/fade_background_${TRAIT}.txt"
fi

# ── 2) Re-run scoring_compute.R with the fresh fade summaries ───────────────
if [ "$RUN_SCORING" = "true" ]; then
echo ""
echo "─── [2/7] scoring_compute.R ───────────────────────────────────"
WD="${WORK}/scoring_compute_v2"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp "${REPO_DIR}/subworkflows/SCORING/local/src/scoring_compute.R" .

Rscript scoring_compute.R \
    --postproc       "${BASE}/postproc/gene_filtering/filtered_discovery.tsv" \
    --fade_top       "${BASE}/selection/fade/top/fade_summary_top.tsv" \
    --fade_bottom    "${BASE}/selection/fade/bottom/fade_summary_bottom.tsv" \
    --fade_site_top  'NO_FADE_SITE_TOP' \
    --fade_site_bot  'NO_FADE_SITE_BOT' \
    --rer            "${BASE}/rerconverge/rer_results/rerconverge_summary_${TRAIT}.tsv" \
    --accum_dir      "${BASE}/accumulation/all/randomization" \
    --stress               true \
    --stress_top_n         25 \
    --top_pct              0.10 \
    --gene_top_pct         0.10 \
    || fail "scoring_compute.R failed"

[ -f gene_scores.tsv ] || fail "scoring_compute.R produced no gene_scores.tsv"
cp *.tsv "${BASE}/scoring/"
rm -rf "${BASE}/scoring/gene_lists"
cp -R gene_lists "${BASE}/scoring/gene_lists"
echo "[2/7] scoring_compute.R done: $(wc -l < gene_scores.tsv) rows in gene_scores.tsv"
else
    echo "[2/7] RUN_SCORING=false — skipping, reusing existing scoring/gene_scores.tsv + fcs_stats*.tsv"
fi

# ── 3) RER FCS report ────────────────────────────────────────────────────────
if [ "$RUN_RER_FCS" = "true" ]; then
echo ""
echo "─── [3/7] RER FCS report ──────────────────────────────────────"
WD="${WORK}/rer_fcs_v3"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .

cat > render_rer_fcs.R <<RSCRIPT
rmarkdown::render(
    '13.FCS_general_report.Rmd',
    params = list(
        stats_file    = '${BASE}/scoring/fcs_stats_rer.tsv',
        universe_file = '${BASE}/rerconverge/gene_lists/background.txt',
        gmt_dir       = '${GMT_DIR}',
        project_name  = '13.FCS_rer_${TRAIT}',
        num_g         = 5,
        fdr_thr       = 0.15,
        pperm_thr     = 0.025,
        top_n         = 20,
        traitname     = '${TRAIT}',
        perms_file    = '${BASE}/rerconverge/rer_results/${TRAIT}.continuous.perms.rds',
        enrich        = TRUE,
        annot_file    = '${BASE}/scoring/fcs_stats.tsv'
    ),
    output_file = '13.FCS_rer_${TRAIT}.html'
)
RSCRIPT

Rscript render_rer_fcs.R || fail "RER FCS report failed"
[ -f "13.FCS_rer_${TRAIT}.html" ] || fail "RER FCS report produced no html"

mkdir -p "${BASE}/scoring/rer/fcs_results"
cp "13.FCS_rer_${TRAIT}.html" "${BASE}/scoring/rer/"
cp "13.FCS_rer_${TRAIT}.html" "${BASE}/html_reports/"
rm -rf "${BASE}/scoring/rer/fcs_results"
mkdir -p "${BASE}/scoring/rer/fcs_results"
cp -R fcs_results/* "${BASE}/scoring/rer/fcs_results/"
echo "[3/7] RER FCS done"
else
    echo "[3/7] RUN_RER_FCS=false — skipping, reusing existing scoring/rer/fcs_results/"
fi

# ── 4) FADE FCS report ───────────────────────────────────────────────────────
if [ "$RUN_FADE_FCS" = "true" ]; then
echo ""
echo "─── [4/7] FADE FCS report ─────────────────────────────────────"
WD="${WORK}/fade_fcs_v3"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .

cat > render_fade_fcs.R <<RSCRIPT
rmarkdown::render(
    '13.FCS_general_report.Rmd',
    params = list(
        stats_file    = '${BASE}/scoring/fcs_stats_fade.tsv',
        universe_file = '${WORK}/universes/fade_background_${TRAIT}.txt',
        gmt_dir       = '${GMT_DIR}',
        project_name  = '13.FCS_fade_${TRAIT}',
        num_g         = 5,
        fdr_thr       = 0.15,
        pperm_thr     = 0.025,
        top_n         = 20,
        traitname     = '${TRAIT}',
        enrich        = TRUE,
        annot_file    = '${BASE}/scoring/fcs_stats.tsv'
    ),
    output_file = '13.FCS_fade_${TRAIT}.html'
)
RSCRIPT

Rscript render_fade_fcs.R || fail "FADE FCS report failed"
[ -f "13.FCS_fade_${TRAIT}.html" ] || fail "FADE FCS report produced no html"

mkdir -p "${BASE}/scoring/fade/fcs_results"
cp "13.FCS_fade_${TRAIT}.html" "${BASE}/scoring/fade/"
cp "13.FCS_fade_${TRAIT}.html" "${BASE}/html_reports/"
rm -rf "${BASE}/scoring/fade/fcs_results"
mkdir -p "${BASE}/scoring/fade/fcs_results"
cp -R fcs_results/* "${BASE}/scoring/fade/fcs_results/"
echo "[4/7] FADE FCS done"
else
    echo "[4/7] RUN_FADE_FCS=false — skipping, reusing existing scoring/fade/fcs_results/"
fi

# ── 5) STRING report ─────────────────────────────────────────────────────────
if [ "$RUN_STRING" = "true" ]; then
echo ""
echo "─── [5/7] STRING report ───────────────────────────────────────"
WD="${WORK}/string_v1"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .

# Same slice_*.tsv -> *.txt conversion SCORING_STRING_REPORT does (first
# column, header dropped) so the Rmd's auto-detected gene_list_files picks
# them up from the working directory.
for f in "${BASE}/scoring/gene_lists"/slice_*.tsv; do
    [ -f "$f" ] || continue
    base="$(basename "$f" .tsv)"
    name="${base#slice_}"
    tail -n +2 "$f" | cut -f1 | { grep -v "^[[:space:]]*$" || true; } > "${name}.txt"
done

cat > render_string.R <<RSCRIPT
rmarkdown::render(
    '15.STRING_report.Rmd',
    params = list(
        background_file       = '${BASE}/postproc/cleaned_backgrounds/cleaned_background_main.txt',
        background_basename   = 'cleaned_background_main.txt',
        output_dir            = 'string_results',
        project_name          = '${TRAIT}',
        species                = 9606,
        required_score         = 400,
        network_score_thr      = 700,
        fdr_thr                = 0.15,
        top_thr                 = 15,
        enable_ppi_enrichment  = TRUE,
        gene_scores_file       = '${BASE}/scoring/gene_scores.tsv',
        string_db_dir          = '${STRING_DB_DIR}'
    ),
    output_file = '15.STRING_report_${TRAIT}.html'
)
RSCRIPT

Rscript render_string.R || fail "STRING report failed"
[ -f "15.STRING_report_${TRAIT}.html" ] || fail "STRING report produced no html"

mkdir -p "${BASE}/string"
cp "15.STRING_report_${TRAIT}.html" "${BASE}/string/"
cp "15.STRING_report_${TRAIT}.html" "${BASE}/html_reports/"
for d in string_results string_summary string_plots string_networks; do
    [ -d "$d" ] || continue
    mkdir -p "${BASE}/string/$d"
    cp -R "$d"/* "${BASE}/string/$d/" 2>/dev/null || true
done
# Compat rename (mirrors SCORING_STRING_REPORT's compat_cmd) so Comparison
# report's cmp_string staging below can pick these up.
mkdir -p string_results
if [ -d string_summary ]; then
    for f in string_summary/*_results.tsv; do
        [ -f "$f" ] || continue
        base="$(basename "$f" _results.tsv)"
        cp "$f" "string_results/${base}_enrichment.tsv"
    done
fi
mkdir -p "${BASE}/string/string_results"
cp string_results/*_enrichment.tsv "${BASE}/string/string_results/" 2>/dev/null || true
echo "[5/7] STRING report done"
else
    echo "[5/7] RUN_STRING=false — skipping, reusing existing string/string_results/ (if any) for the Comparison report"
fi

# ── 6) POSENRICH: re-run the Fisher computation (fold-gate removed), then report ──
if [ "$RUN_POSENRICH" = "true" ]; then
echo ""
echo "─── [6/7] POSENRICH ────────────────────────────────────────────"
WD="${WORK}/posenrich_v1"
mkdir -p "$WD"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .

# Reuse the already-built position GMTs (external-DB/alignment derived —
# untouched by any of this session's RER/FADE/scoring fixes, no need to
# rebuild them).
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
echo "[6/7] POSENRICH done"
else
    echo "[6/7] RUN_POSENRICH=false — skipping, reusing existing posenrich/posenrich_characterization.tsv"
fi

# ── 7) Comparison report ─────────────────────────────────────────────────────
# CAAS + accumulation fcs_all_results.tsv reused as-is (unaffected by the
# RER/FADE-specific fixes); RER + FADE come from steps 3-4, STRING from step 5.
if [ "$RUN_COMPARE" = "true" ]; then
echo ""
echo "─── [7/7] Comparison report ───────────────────────────────────"
WD="${WORK}/compare_v3"
mkdir -p "$WD/cmp_fcs" "$WD/cmp_string"
cd "$WD" || fail "cannot cd $WD"
cp -R "${REPO_DIR}/subworkflows/ENRICHMENT/local/"* .
mkdir -p cmp_fcs cmp_string

cp "${BASE}/fcs/fcs_results/fcs_results/fcs_all_results.tsv"        cmp_fcs/caas_fcs_all.tsv
cp "${BASE}/scoring/rer/fcs_results/fcs_all_results.tsv"            cmp_fcs/rer_fcs_all.tsv
cp "${BASE}/scoring/fade/fcs_results/fcs_all_results.tsv"           cmp_fcs/fade_fcs_all.tsv
cp "${BASE}/scoring/accum/fcs_results/fcs_all_results.tsv"          cmp_fcs/accum_fcs_all.tsv
cp "${BASE}/string/string_results/"*_enrichment.tsv cmp_string/ 2>/dev/null || true

cat > render_compare.R <<RSCRIPT
rmarkdown::render(
    '12.Comparison_report.Rmd',
    params = list(
        fcs_dir    = 'cmp_fcs',
        string_dir = 'cmp_string',
        fdr_thr    = 0.15,
        pperm_thr  = 0.025,
        top_n      = 20,
        traitname  = '${TRAIT}'
    ),
    output_file = '12.Comparison_report_${TRAIT}.html'
)
RSCRIPT

Rscript render_compare.R || fail "Comparison report failed"
[ -f "12.Comparison_report_${TRAIT}.html" ] || fail "Comparison report produced no html"

mkdir -p "${BASE}/compare"
cp "12.Comparison_report_${TRAIT}.html" "${BASE}/compare/"
cp "12.Comparison_report_${TRAIT}.html" "${BASE}/html_reports/"
echo "[7/7] Comparison report done"
else
    echo "[7/7] RUN_COMPARE=false — skipping"
fi

echo ""
echo "======================================================"
echo " ALL DONE: $TRAIT"
echo "======================================================"
