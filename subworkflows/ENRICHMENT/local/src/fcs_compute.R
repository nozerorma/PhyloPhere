#!/usr/bin/env Rscript
# =============================================================================
# fcs_compute.R — standalone FCS stats entry point (Nextflow-batchable)
# =============================================================================
# Extracted from 12.FCS_general_report.Rmd's `run` chunk so the expensive part
# — fcs_run_all()'s Wilcoxon-AUC / Lachenbruch / Path-Sum-Permulation tests
# over every GMT database — can run as N independent Nextflow tasks batched by
# GMT file, instead of one monolithic in-Rmd computation. Safe to batch: BH
# correction in fcs_enrich.R is scoped per database throughout (see its own
# header comment + fcs_run_ranking/fcs_run_lachenbruch/fcs_run_permulation),
# so a batch's rows never need reconciling against another batch's — a plain
# row-concat of every batch's output (see FCS_CONCAT in fcs.nf) is exact.
#
# Deliberately stops short of the GMT-description join and evidence_score
# percentile-rank step: `description` comes from the FULL (unbatched) GMT
# directory in 12.FCS_general_report.Rmd itself (cheap, and desc_map is reused
# a second time later in that Rmd — see its own comments), and evidence_score
# ranks each term's statistic as a percentile among ALL terms tested in that
# ranking across every database, so it can only be computed once, after every
# batch's rows are merged back together — also left to the Rmd.
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

# ── minimal flag parser (see subworkflows/SCORING/local/src/scoring_caas_perms.R) ──
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

stats_file    <- get_arg("--stats-file")
universe_file <- get_arg("--universe-file", "NO_FILE")
gmt_dir       <- get_arg("--gmt-dir")
perms_file    <- get_arg("--perms-file", "NO_FILE")
num_g         <- as.numeric(get_arg("--num-g", "10"))
max_g         <- as.numeric(get_arg("--max-g", "0"))
fdr_thr       <- as.numeric(get_arg("--fdr-thr", "0.15"))
fdr_wilcoxon  <- as.numeric(get_arg("--fdr-wilcoxon", as.character(fdr_thr)))
fdr_lachenbruch <- as.numeric(get_arg("--fdr-lachenbruch", as.character(fdr_thr)))
fdr_permsum   <- as.numeric(get_arg("--fdr-permsum", as.character(fdr_thr)))
pperm_thr     <- as.numeric(get_arg("--pperm-thr", "0.025"))
n_perms_sum   <- as.numeric(get_arg("--n-perms-sum", "10000"))
seed_val      <- as.integer(get_arg("--seed", "1998"))
out_file      <- get_arg("--output", "fcs_enrich_partial.tsv")

stopifnot(!is.null(stats_file), file.exists(stats_file))
stopifnot(!is.null(gmt_dir), dir.exists(gmt_dir))
set.seed(seed_val)

this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
src_candidates <- c(file.path(this_dir, "fcs_enrich.R"), "src/fcs_enrich.R", "fcs_enrich.R")
src <- src_candidates[file.exists(src_candidates)][1]
if (is.na(src)) stop("fcs_enrich.R not found next to fcs_compute.R or under src/")
suppressMessages(invisible(capture.output(source(src))))

# ── stats / universe / rankings — identical logic to 12.FCS_general_report.Rmd's `load`/`run` chunks ──
stats <- readr::read_tsv(stats_file, show_col_types = FALSE)
if (!"gene" %in% names(stats)) {
  gcol <- intersect(c("Gene", "GENE"), names(stats))[1]
  if (!is.na(gcol)) stats <- dplyr::rename(stats, gene = !!gcol)
}
stopifnot("gene" %in% names(stats))

universe <- character(0)
if (universe_file != "NO_FILE" && file.exists(universe_file) && !grepl("^NO_", basename(universe_file))) {
  universe <- unique(trimws(readLines(universe_file)))
  universe <- universe[nzchar(universe)]
}
if (length(universe) < 2) universe <- unique(stats$gene)

score_cols <- grep("^score_", names(stats), value = TRUE)
if (length(score_cols) == 0) stop("stats_file has no score_<ranking> columns")

rankings <- list()
for (sc in score_cols) {
  rk <- sub("^score_", "", sc)
  rankings[[rk]] <- fcs_build_vals(setNames(stats[[sc]], stats$gene), universe)
}

gmts <- fcs_load_gmts(gmt_dir)

enrich <- fcs_run_all(rankings, gmts, num_g = num_g, max_g = max_g,
                      perms_file = perms_file, fdr_thr = fdr_thr,
                      fdr_wilcoxon = fdr_wilcoxon, fdr_lachenbruch = fdr_lachenbruch,
                      fdr_permsum = fdr_permsum, p_perm_thr = pperm_thr,
                      n_perms_sum = n_perms_sum)

readr::write_tsv(enrich, out_file)
