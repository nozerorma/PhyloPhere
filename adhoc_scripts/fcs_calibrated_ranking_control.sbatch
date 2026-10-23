#!/bin/bash
# WP9: WP9_CONTROL=true uses wp9_gene_stat_control.R to compare every requested variant and every pseudo-observed cycle.
# One job per trait/mode/ranking, no array required. STORE, UNIVERSE and WP9_OUT may override derived paths.
# The WP9 path archives HEAD to a new stable snapshot. It deliberately uses whole-memory validation matrices.
# At 10k, request memory separately (e.g. --mem=64G) and explicitly set MAX_CELLS=200000000 per ranking.
# Whole-memory validation control: explicitly reconstructs a compatibility null from canonical blocks when available.
# fcs_calibrated_ranking_control.sbatch — Does the CAAS FCS separate the observed run from null cycles when the gene ranking is calibrated per gene?
#
# Productive CAAS FCS ranks -log10 of the per-gene permulation p. The base statistic
# (max, F(max), F(max)^n) is selected by gene_stat. The legacy per-cycle control below
# calibrates the canonical statistic with the same inclusive-tie and leave-one-out rules.
# WP9_CONTROL compares the base-statistic variants without changing the productive default.
# It adds disjoint reference/evaluation folds because a marginal LOO p-rate alone is
# largely determined by ranks and cannot validate biological calibration.
#
# Two modes, as in fcs_null_as_observed.sbatch:
#   OBS=0 (default): a null cycle is the observed (its ranking against the other 999); one array task per picked cycle, picks stratified by level.
#   OBS=1: the observed run itself (its gene scores against the 1000 cycles); a one-task array.
# Submit per trait and mode:
#   ssh correfoc "mkdir -p ~/scratch/0.Phylophere/.tmp/fcs_calib_ctrl && bash -lc 'sbatch --export=ALL,TRAIT=neoplasia,MODE=strict,OBS=1 --array=1-1'" < adhoc_scripts/fcs_calibrated_ranking_control.sbatch
#   ssh correfoc "bash -lc 'sbatch --export=ALL,TRAIT=neoplasia,MODE=strict --array=1-20%8'" < adhoc_scripts/fcs_calibrated_ranking_control.sbatch
# TRAIT: malignant | neoplasia. MODE: strict | complete. Optional: N_PER (default 5), SEED (default 1), NUM_G and MAX_G (defaults 10 and 0 (no upper bound)).
# Legacy outputs are under .tmp/fcs_calib_ctrl; WP9 outputs/snapshots default to .tmp/wp9/.
#
#SBATCH --job-name=fcs_calib_ctrl
#SBATCH --partition=std-cpu
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH -o /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_calib_ctrl/slurm-%A_%a.out
#SBATCH -e /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_calib_ctrl/slurm-%A_%a.err

set -euo pipefail
: "${TRAIT:?set TRAIT (malignant|neoplasia)}" "${MODE:?set MODE (strict|complete)}"
export TRAIT MODE
export OBS="${OBS:-0}" N_PER="${N_PER:-5}" SEED="${SEED:-1}" NUM_G="${NUM_G:-10}" MAX_G="${MAX_G:-0}"
if [ "${WP9_CONTROL:-false}" = true ]; then
    case "${SLURM_ARRAY_TASK_ID:-1}" in 1) ;; *) echo 'WP9_CONTROL computes every cycle in one job; do not use an array' >&2; exit 2;; esac
    wp9_root="${PHYLOPHERE_ROOT:-$HOME/scratch/0.Phylophere}"
    wp9_results="$HOME/scratch/2.Primates/2.Primates_results/CAAS_RESULTS"
    wp9_complete="$wp9_results/${TRAIT}_def/${TRAIT}_prevalence_complete"
    wp9_mode="$wp9_results/${TRAIT}_def/${TRAIT}_prevalence_${MODE}"
    wp9_store="$wp9_complete/caas_permulation"
    if [ "$MODE" = strict ]; then wp9_store="$wp9_complete/strict/caas_permulation"; fi
    wp9_rankings="${RANKINGS:-global}"
    wp9_out="${WP9_OUT:-$wp9_root/.tmp/wp9/${TRAIT}_${MODE}_${wp9_rankings//,/_}_${SLURM_JOB_ID:?}}"
    [ ! -e "$wp9_out" ] && [ ! -e "${wp9_out}_code" ] || { echo 'WP9 output/snapshot already exists' >&2; exit 2; }
    mkdir -p "${wp9_out}_code" "${wp9_out}_tmp"
    git -C "$wp9_root" archive HEAD | tar -x -C "${wp9_out}_code"
    git -C "$wp9_root" rev-parse HEAD > "${wp9_out}_code/commit.txt"
    export TMPDIR="${wp9_out}_tmp"
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
    wp9_args=( --store "${STORE:-$wp9_store}" --universe "${UNIVERSE:-$wp9_mode/universes/cleaned_background_main_${TRAIT}_prevalence.txt}"
       --out "$wp9_out" --rankings "$wp9_rankings" --variants "${VARIANTS:-max,fmax,beta025,beta05,beta075,size_adj_max}"
       --min-genes "$NUM_G" --max-genes "$MAX_G" --max-cells "${MAX_CELLS:-50000000}" )
    if [ -n "${EXPECTED_CYCLES:-}" ]; then wp9_args+=( --expected-cycles "$EXPECTED_CYCLES" ); fi
    if [ -n "${DATABASES:-}" ]; then wp9_args+=( --databases "$DATABASES" ); fi
    /usr/bin/time -v "${RSCRIPT:-$HOME/micromamba/envs/phylophere/bin/Rscript}" \
      "${wp9_out}_code/adhoc_scripts/wp9_gene_stat_control.R" "${wp9_args[@]}"
    exit 0
fi

export TMPDIR="$HOME/scratch/0.Phylophere/.tmp/fcs_calib_ctrl/tmp"
mkdir -p "$TMPDIR"
# Env R called directly: conda is not on the compute-node PATH
/homes/users/mramon/micromamba/envs/phylophere/bin/Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({ library(RERconverge); library(Matrix) })
env <- function(k) Sys.getenv(k)
trait <- env("TRAIT"); mode <- env("MODE"); obs_mode <- env("OBS") == "1"; n_per <- as.integer(env("N_PER")); seed <- as.integer(env("SEED"))
num_g <- as.integer(env("NUM_G")); max_g <- as.integer(env("MAX_G")); k <- as.integer(env("SLURM_ARRAY_TASK_ID"))
home <- Sys.getenv("HOME"); root <- file.path(home, "scratch/0.Phylophere")
res  <- file.path(home, "scratch/2.Primates/2.Primates_results/CAAS_RESULTS")
run_complete <- file.path(res, sprintf("%s_def/%s_prevalence_complete", trait, trait))
run_mode     <- file.path(res, sprintf("%s_def/%s_prevalence_%s", trait, trait, mode))
perms_file   <- if (mode == "strict") file.path(run_complete, "strict/caas_permulation/caas_perms.rds") else file.path(run_complete, "caas_permulation/caas_perms.rds")
universe_file <- file.path(run_mode, "universes", sprintf("cleaned_background_main_%s_prevalence.txt", trait))
out_dir <- file.path(root, ".tmp/fcs_calib_ctrl", paste0(trait, "_", mode)); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); setwd(out_dir)

src <- file.path(root, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "evidence_flags.R")); source(file.path(src, "fcs_enrich.R"))
gmts <- fcs_load_gmts(file.path(root, "subworkflows/ENRICHMENT/dat"))
universe <- unique(trimws(readLines(universe_file))); universe <- universe[nzchar(universe)]
source(file.path(root, "subworkflows/SCORING/local/src/caas_block_contract.R"))
x <- caas_control_null(dirname(perms_file)); M <- x$caas_corStat_byrank; rks <- c("global", "top", "bottom"); N <- ncol(M$global)
gs <- read.delim(file.path(run_mode, "scoring/gene_scores.tsv"), stringsAsFactors = FALSE)
obs_col <- c(global = "gene_caas_score", top = "gene_caas_score_top", bottom = "gene_caas_score_bottom")
obs_level <- sum(gs$gene_caas_score > 0, na.rm = TRUE); level <- colSums(M$global > 0, na.rm = TRUE)

# Use the productive inclusive-tie, pseudocount and leave-one-out rules.
source(file.path(src, "gene_p.R"))
loo_log10p <- function(X) gene_p(X, rep(0, nrow(X)))$nlog10_null
obs_log10p <- function(s, X) gene_p(X, s)$nlog10_obs


if (obs_mode) { tag <- "observed"; stratum <- "observed"; j <- NA_integer_ } else {
  set.seed(seed); q <- quantile(level, c(1 / 3, 2 / 3))
  strata <- list(low = which(level <= q[1]), mid = which(level > q[1] & level <= q[2]),
                 high = which(level > q[2] & level < 0.8 * obs_level), observed_like = which(level >= 0.8 * obs_level & level <= 1.5 * obs_level))
  picks <- unlist(lapply(names(strata), function(s) { i <- strata[[s]]; i <- sample(i, min(n_per, length(i))); setNames(i, rep(s, length(i))) }))
  cat(sprintf("%s %s | cycles %d | observed %d scored genes | %d picks | task %d\n", trait, mode, N, obs_level, length(picks), k))
  if (k > length(picks)) { cat("no pick for this task\n"); quit(status = 0) }
  j <- picks[[k]]; tag <- colnames(M$global)[j]; stratum <- names(picks)[k] }
out <- file.path(out_dir, sprintf("%s.rds", tag)); if (file.exists(out)) { cat("exists:", out, "\n"); quit(status = 0) }

t0 <- Sys.time(); rankings <- list(); nulls <- list()
for (r in rks) {
  g <- rownames(M[[r]]); L <- loo_log10p(M[[r]])
  if (obs_mode) { so <- setNames(gs[[obs_col[[r]]]], gs$Gene)[g]; v <- obs_log10p(so, M[[r]]); nulls[[r]] <- L }
  else { v <- L[, j]; nulls[[r]] <- L[, -j, drop = FALSE] }
  v <- setNames(v, g); rankings[[r]] <- fcs_build_vals(v[v > 0], universe); rm(L); gc(verbose = FALSE)
}
cat(sprintf("gene rankings built in %.1f min; nonzero genes: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")),
            paste(sprintf("%s %d", rks, sapply(rankings, function(v) sum(v > 0))), collapse = ", ")))

dummy <- sprintf("MEM_%s", tag); invisible(file.create(dummy))
MEMOBJ <- list(caas_corStat_byrank = nulls, gene_stat = x$gene_stat)
readRDS <- function(file, ...) if (identical(file, dummy)) MEMOBJ else base::readRDS(file, ...)
r <- fcs_run_all(rankings, gmts, num_g = num_g, max_g = max_g, perms_file = dummy, fdr_thr = 0.15, p_perm_thr = 0.025)
keep <- intersect(c("ranking", "database", "pathway", "stat", "pval", "p.adj", "p.perm", "q.perm", "num.genes", "sig_perm", "redlich_level", "direction"), names(r))
saveRDS(list(trait = trait, mode = mode, cycle = tag, stratum = stratum, level = if (obs_mode) obs_level else level[j], obs_level = obs_level,
             num_g = num_g, max_g = max_g, res = as.data.frame(r)[, keep]), out)
cat(sprintf("%s (%s): %d rows, %.1f min total -> %s\n", tag, stratum, nrow(r), as.numeric(difftime(Sys.time(), t0, units = "mins")), out))
unlink(dummy)
RSCRIPT
