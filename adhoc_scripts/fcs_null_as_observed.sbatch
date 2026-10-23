#!/bin/bash
# Whole-memory validation control: explicitly reconstructs a compatibility null from canonical blocks when available.
# fcs_null_as_observed.sbatch — Is the permulation p of the CAAS FCS calibrated when the observed and the null cycles differ in level?
#
# A null cycle is taken as the observed ranking (global, top, bottom gene scores) and the other 999 cycles are its null, through
# fcs_run_all of the pipeline (the checkout's own fcs_enrich.R, unchanged). If a cycle is exchangeable with the others, the
# permulation p of every set is uniform whatever the width of the null. The picked cycles are stratified by their number of scored
# genes (low, mid, high, and observed_like = 0.8 to 1.5 times the observed run's scored genes), so the p distribution can be read
# against the level.
#
# One array task per picked cycle; the picks are a function of SEED, N_PER and the null only. Submit one array per trait and mode
# (a one-task array first, to check paths and memory):
#   ssh correfoc "mkdir -p ~/scratch/0.Phylophere/.tmp/fcs_null_ctrl && bash -lc 'sbatch --export=ALL,TRAIT=neoplasia,MODE=strict --array=1-1'" < adhoc_scripts/fcs_null_as_observed.sbatch
#   ssh correfoc "bash -lc 'sbatch --export=ALL,TRAIT=neoplasia,MODE=strict --array=1-20%10'" < adhoc_scripts/fcs_null_as_observed.sbatch
# TRAIT: malignant | neoplasia. MODE: strict | complete. Optional: N_PER (cycles per stratum, default 5), SEED (default 1),
# NUM_G and MAX_G (set size limits, defaults 10 and 1000). Nothing is written to the cluster besides the outputs under
# ~/scratch/0.Phylophere/.tmp/fcs_null_ctrl/<trait>_<mode>/ (one <cycle>.rds per task: cycle, stratum, level, res).
#
#SBATCH --job-name=fcs_null_obs
#SBATCH --partition=std-cpu
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -o /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_null_ctrl/slurm-%A_%a.out
#SBATCH -e /homes/users/mramon/scratch/0.Phylophere/.tmp/fcs_null_ctrl/slurm-%A_%a.err

set -euo pipefail
: "${TRAIT:?set TRAIT (malignant|neoplasia)}" "${MODE:?set MODE (strict|complete)}"
export TRAIT MODE
export N_PER="${N_PER:-5}" SEED="${SEED:-1}" NUM_G="${NUM_G:-10}" MAX_G="${MAX_G:-1000}"
export TMPDIR="$HOME/scratch/0.Phylophere/.tmp/fcs_null_ctrl/tmp"
mkdir -p "$TMPDIR"
# Env R called directly: conda is not on the compute-node PATH
/homes/users/mramon/micromamba/envs/phylophere/bin/Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({ library(RERconverge); library(Matrix) })
env <- function(k) Sys.getenv(k)
trait <- env("TRAIT"); mode <- env("MODE"); n_per <- as.integer(env("N_PER")); seed <- as.integer(env("SEED"))
num_g <- as.integer(env("NUM_G")); max_g <- as.integer(env("MAX_G")); k <- as.integer(env("SLURM_ARRAY_TASK_ID"))
home <- Sys.getenv("HOME")
root <- file.path(home, "scratch/0.Phylophere")
res  <- file.path(home, "scratch/2.Primates/2.Primates_results/CAAS_RESULTS")
run_complete <- file.path(res, sprintf("%s_def/%s_prevalence_complete", trait, trait))
run_mode     <- file.path(res, sprintf("%s_def/%s_prevalence_%s", trait, trait, mode))
perms_file   <- if (mode == "strict") file.path(run_complete, "strict/caas_permulation/caas_perms.rds") else file.path(run_complete, "caas_permulation/caas_perms.rds")
universe_file <- file.path(run_mode, "universes", sprintf("cleaned_background_main_%s_prevalence.txt", trait))
out_dir <- file.path(root, ".tmp/fcs_null_ctrl", paste0(trait, "_", mode))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); setwd(out_dir)

src <- file.path(root, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "evidence_flags.R")); source(file.path(src, "fcs_enrich.R"))
gmts <- fcs_load_gmts(file.path(root, "subworkflows/ENRICHMENT/dat"))
universe <- unique(trimws(readLines(universe_file))); universe <- universe[nzchar(universe)]

source(file.path(root, "subworkflows/SCORING/local/src/caas_block_contract.R"))
x <- caas_control_null(dirname(perms_file)); M <- x$caas_corStat_byrank
level <- colSums(M$global > 0, na.rm = TRUE)                       # scored genes per cycle (global ranking)
gs <- read.delim(file.path(run_mode, "scoring/gene_scores.tsv"), stringsAsFactors = FALSE)
obs_level <- sum(gs$gene_caas_score > 0, na.rm = TRUE)             # scored genes of the observed run (global ranking)
set.seed(seed)
q <- quantile(level, c(1 / 3, 2 / 3))
strata <- list(low = which(level <= q[1]), mid = which(level > q[1] & level <= q[2]),
               high = which(level > q[2] & level < 0.8 * obs_level),
               observed_like = which(level >= 0.8 * obs_level & level <= 1.5 * obs_level))
picks <- unlist(lapply(names(strata), function(s) { i <- strata[[s]]; i <- sample(i, min(n_per, length(i))); setNames(i, rep(s, length(i))) }))
cat(sprintf("%s %s | null %s | cycles %d | level q33 %.0f q67 %.0f | observed %d | %d picks | task %d\n", trait, mode, basename(perms_file),
            ncol(M$global), q[1], q[2], obs_level, length(picks), k))
if (k > length(picks)) { cat("no pick for this task\n"); quit(status = 0) }

j <- picks[[k]]; cyc <- colnames(M$global)[j]; out <- file.path(out_dir, sprintf("%s.rds", cyc))
if (file.exists(out)) { cat("exists:", out, "\n"); quit(status = 0) }

# The null reaches fcs_run_all in memory: readRDS of a dummy file name is shadowed, so no copy of the null is written.
dummy <- sprintf("MEM_%d", k); invisible(file.create(dummy))
MEMOBJ <- list(caas_corStat_byrank = lapply(M, function(m) m[, -j, drop = FALSE]), gene_stat = x$gene_stat)
readRDS <- function(file, ...) if (identical(file, dummy)) MEMOBJ else base::readRDS(file, ...)

t0 <- Sys.time()
rankings <- lapply(setNames(c("global", "top", "bottom"), c("global", "top", "bottom")), function(r) {
  s <- M[[r]][, j]; names(s) <- rownames(M[[r]]); s <- s[!is.na(s) & s > 0]
  fcs_build_vals(s, universe)
})
r <- fcs_run_all(rankings, gmts, num_g = num_g, max_g = max_g, perms_file = dummy, fdr_thr = 0.15, p_perm_thr = 0.025)
keep <- intersect(c("ranking", "database", "pathway", "stat", "pval", "p.adj", "p.perm", "num.genes", "lach_pval", "lach_p.adj",
                    "lach_p.perm", "sig_wilcoxon", "sig_lachenbruch", "evidence_label"), names(r))
saveRDS(list(trait = trait, mode = mode, cycle = cyc, stratum = names(picks)[k], level = level[j], obs_level = obs_level,
             num_g = num_g, max_g = max_g, res = as.data.frame(r)[, keep]), out)
cat(sprintf("%s (%s, %d scored genes): %d rows, %.1f min -> %s\n", cyc, names(picks)[k], level[j], nrow(r),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out))
unlink(dummy)
RSCRIPT
