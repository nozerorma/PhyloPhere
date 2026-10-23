#!/bin/bash
# Whole-memory validation control: explicitly reconstructs a compatibility null from canonical blocks when available.
# fcs_rawmax_ranking_control.sbatch — Same control as fcs_calibrated_ranking_control.sbatch, with the raw position score as the gene statistic.
#
# The gene statistic is the maximum raw position score of the gene in that direction (global: best side; top, bottom: that side only), rebuilt
# from perm_pos_cycle_caas.tsv.gz (null cycles) and b0/perm_pos_cycle_caas.tsv.gz (observed), instead of size_adj_max = F(max)^n. F(max)^n is a quantile
# within its own labeling, so it cannot carry the overall strength of a labeling. The ranking is the gene permulation p against the other cycles,
#     p_gk = (1 + #{null cycles l != k with max_gl >= max_gk}) / N   if the gene is scored in k, 1 otherwise,   ranking = -log10(p_gk),
# and the null of the set tests is the same transformation of the other cycles. Cycle picks are those of the calibrated control (same seed, same detection).
# The rest is fcs_run_all of the checkout, unchanged. Modes and submission as in fcs_calibrated_ranking_control.sbatch:
#   ssh correfoc "mkdir -p ~/scratch/0.Phylophere/.tmp/fcs_rawmax_ctrl && bash -lc 'sbatch --export=ALL,TRAIT=neoplasia,MODE=strict,OBS=1 --array=1-1'" < adhoc_scripts/fcs_rawmax_ranking_control.sbatch
# Complete runs read ~28M position rows: submit them with --mem=24G.
# Nothing is written to the cluster besides the outputs under ~/scratch/0.Phylophere/.tmp/fcs_rawmax_ctrl/<trait>_<mode>/.
#
#SBATCH --job-name=fcs_rawmax_ctrl
#SBATCH --partition=std-cpu
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH -o /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_rawmax_ctrl/slurm-%A_%a.out
#SBATCH -e /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_rawmax_ctrl/slurm-%A_%a.err

set -euo pipefail
: "${TRAIT:?set TRAIT (malignant|neoplasia)}" "${MODE:?set MODE (strict|complete)}"
export TRAIT MODE
export OBS="${OBS:-0}" N_PER="${N_PER:-5}" SEED="${SEED:-1}" NUM_G="${NUM_G:-10}" MAX_G="${MAX_G:-1000}"
export TMPDIR="$HOME/scratch/0.Phylophere/.tmp/fcs_rawmax_ctrl/tmp"
mkdir -p "$TMPDIR"
# Env R called directly: conda is not on the compute-node PATH
/homes/users/mramon/micromamba/envs/phylophere/bin/Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({ library(RERconverge); library(Matrix); library(data.table) })
env <- function(k) Sys.getenv(k)
trait <- env("TRAIT"); mode <- env("MODE"); obs_mode <- env("OBS") == "1"; n_per <- as.integer(env("N_PER")); seed <- as.integer(env("SEED"))
num_g <- as.integer(env("NUM_G")); max_g <- as.integer(env("MAX_G")); k <- as.integer(env("SLURM_ARRAY_TASK_ID"))
home <- Sys.getenv("HOME"); root <- file.path(home, "scratch/0.Phylophere")
res  <- file.path(home, "scratch/2.Primates/2.Primates_results/CAAS_RESULTS")
run_complete <- file.path(res, sprintf("%s_def/%s_prevalence_complete", trait, trait))
run_mode     <- file.path(res, sprintf("%s_def/%s_prevalence_%s", trait, trait, mode))
perms_file   <- if (mode == "strict") file.path(run_complete, "strict/caas_permulation/caas_perms.rds") else file.path(run_complete, "caas_permulation/caas_perms.rds")
universe_file <- file.path(run_mode, "universes", sprintf("cleaned_background_main_%s_prevalence.txt", trait))
out_dir <- file.path(root, ".tmp/fcs_rawmax_ctrl", paste0(trait, "_", mode)); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); setwd(out_dir)

src <- file.path(root, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "percentile_flags.R")); source(file.path(src, "fcs_enrich.R"))
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



# Raw maxima are validation-only; use the exact WP5 parts and b_0 block.
store <- dirname(perms_file)
block_files <- Sys.glob(file.path(store, "blocks/blk_*/matrices.rds"))
if (length(block_files)) {
  R <- lapply(M, function(m) { m[] <- 0; m })
  for (f in block_files) {
    b <- readRDS(f)
    if (is.null(b$parts$max)) stop("CAAS blocks lack raw maximum parts")
    for (r in rks) R[[r]][b$genes,b$cycles] <- b$parts$max[[r]]
  }
  b0 <- readRDS(file.path(store,"b0/matrices.rds"))
  OBS <- lapply(rks, function(r) {
    v <- setNames(numeric(nrow(R[[r]])),rownames(R[[r]]))
    common <- intersect(names(v), b0$genes)
    v[common] <- b0$parts$max[[r]][common,"b_0"]
    v
  }); names(OBS) <- rks
  for (r in rks) if (!identical(R[[r]] > 0, M[[r]] > 0)) stop("Raw-max detection differs from canonical stat")
} else {
# Raw maximum position score per gene, direction and cycle, in the layout of the pipeline matrices (same rows, columns and NA pattern)
pos_dir <- dirname(perms_file)
read_pos <- function(f) { d <- fread(f, select = c("Gene", "Position", "side", "cycle", "caas_score"), showProgress = FALSE)
  d <- d[!is.na(caas_score) & caas_score > 0]; d[, caas_score := round(caas_score, 12)]; d }
allp <- rbind(read_pos(file.path(pos_dir, "b0/perm_pos_cycle_caas.tsv.gz")), read_pos(file.path(pos_dir, "perm_pos_cycle_caas.tsv.gz")))
by_dir <- list(global = allp[, .(m = max(caas_score)), by = .(cycle, Gene)], top = allp[side == "top", .(m = max(caas_score)), by = .(cycle, Gene)],
               bottom = allp[side == "bottom", .(m = max(caas_score)), by = .(cycle, Gene)]); rm(allp); gc(verbose = FALSE)
R <- list(); OBS <- list()
for (r in rks) {
  rows <- rownames(M[[r]]); cols <- colnames(M[[r]]); dt <- by_dir[[r]]
  Mr <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols)); i <- match(dt$Gene, rows); j <- match(dt$cycle, cols); ok <- !is.na(i) & !is.na(j)
  Mr[cbind(i[ok], j[ok])] <- dt$m[ok]; Mr[is.na(M[[r]])] <- NA
  mism <- sum((Mr > 0) != (M[[r]] > 0), na.rm = TRUE); lost <- sum(is.na(i) & dt$cycle != "b_0")
  cat(sprintf("%s: raw-max matrix %d x %d | scored cells %d | detection mismatches vs pipeline matrix %d | genes absent from the matrix %d\n", r, nrow(Mr), ncol(Mr), sum(Mr > 0, na.rm = TRUE), mism, lost))
  if (mism > 0 || lost > 0) stop("raw-max matrix does not reproduce the detection of the pipeline matrix")
  o <- dt[cycle == "b_0"]; ov <- setNames(numeric(length(rows)), rows); io <- match(o$Gene, rows); ov[io[!is.na(io)]] <- o$m[!is.na(io)]; OBS[[r]] <- ov
  R[[r]] <- Mr }
rm(by_dir)
}
stopifnot(sum(OBS$global > 0) == obs_level)

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
  g <- rownames(R[[r]]); L <- loo_log10p(R[[r]])
  if (obs_mode) { v <- obs_log10p(OBS[[r]], R[[r]]); nulls[[r]] <- L }
  else { v <- L[, j]; nulls[[r]] <- L[, -j, drop = FALSE] }
  v <- setNames(v, g); rankings[[r]] <- fcs_build_vals(v[v > 0], universe); rm(L); gc(verbose = FALSE)
}
cat(sprintf("gene rankings built in %.1f min; nonzero genes: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")),
            paste(sprintf("%s %d", rks, sapply(rankings, function(v) sum(v > 0))), collapse = ", ")))

dummy <- sprintf("MEM_%s", tag); invisible(file.create(dummy))
MEMOBJ <- list(caas_corStat_byrank = nulls, gene_stat = x$gene_stat)
readRDS <- function(file, ...) if (identical(file, dummy)) MEMOBJ else base::readRDS(file, ...)
r <- fcs_run_all(rankings, gmts, num_g = num_g, max_g = max_g, perms_file = dummy, fdr_thr = 0.15, p_perm_thr = 0.025)
keep <- intersect(c("ranking", "database", "pathway", "stat", "pval", "p.adj", "p.perm", "num.genes", "lach_pval", "lach_p.adj",
                    "lach_p.perm", "sig_wilcoxon", "sig_lachenbruch", "evidence_label"), names(r))
saveRDS(list(trait = trait, mode = mode, cycle = tag, stratum = stratum, level = if (obs_mode) obs_level else level[j], obs_level = obs_level,
             num_g = num_g, max_g = max_g, res = as.data.frame(r)[, keep]), out)
cat(sprintf("%s (%s): %d rows, %.1f min total -> %s\n", tag, stratum, nrow(r), as.numeric(difftime(Sys.time(), t0, units = "mins")), out))
unlink(dummy)
RSCRIPT
