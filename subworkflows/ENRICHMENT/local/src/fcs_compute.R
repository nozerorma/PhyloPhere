#!/usr/bin/env Rscript
# fcs_compute.R — FCS statistics of the gene sets of one batch of GMT files.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Called by:  FCS_COMPUTE_BATCHED Nextflow process (fcs.nf → Rscript fcs_compute.R ...)
#
# Runs fcs_run_all() (Wilcoxon-AUC and Lachenbruch tests) over
# every GMT database of --gmt-dir, so the expensive part can run as independent tasks
# batched by GMT file. The BH correction of fcs_enrich.R is scoped per database, so the
# rows of one batch never need reconciling with another and a row-concat of all batches
# (FCS_CONCAT in fcs.nf) is exact.
#
# Left to 12.FCS_general_report.Rmd, after the batches are merged:
#   - the GMT description join, which uses the full (unbatched) GMT directory;
#   - evidence_score, which ranks each statistic as a percentile among all the terms
#     tested in that ranking across every database.
#
# Args (named flags, from task.script):
#   --stats-file        gene scores, one score_<ranking> column per ranking
#   --universe-file     gene universe, or NO_FILE (then the genes of --stats-file)
#   --gmt-dir           directory of the *.gmt files of this batch
#   --perms-file        null permutations for the permulation p of both tests, or NO_FILE
#   --num-g, --max-g    minimum and maximum gene-set size (max-g 0 = no limit)
#   --fdr-thr           default FDR; --fdr-wilcoxon and --fdr-lachenbruch override it per test
#   --pperm-thr         permutation p-value threshold
#   --output            output TSV
# =============================================================================

# ── Dependencies ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────

# Minimal named-flag parser: returns the value after `flag`, or `default` when it is absent or last.
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
pperm_thr     <- as.numeric(get_arg("--pperm-thr", "0.025"))
out_file      <- get_arg("--output", "fcs_enrich_partial.tsv")

stopifnot(!is.null(stats_file), file.exists(stats_file))
stopifnot(!is.null(gmt_dir), dir.exists(gmt_dir))


# ── Load fcs_enrich.R ─────────────────────────────────────────────────────────

# fcs_enrich.R sits next to this script in the task directory, or under src/ when run from the module root.
this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
src_candidates <- c(file.path(this_dir, "fcs_enrich.R"), "src/fcs_enrich.R", "fcs_enrich.R")
src <- src_candidates[file.exists(src_candidates)][1]
if (is.na(src)) stop("fcs_enrich.R not found next to fcs_compute.R or under src/")
suppressMessages(invisible(capture.output(source(src))))

# ── Stats, universe and rankings ──────────────────────────────────────────────

# Same logic as the `load` and `run` chunks of 12.FCS_general_report.Rmd.
stats <- readr::read_tsv(stats_file, show_col_types = FALSE)
if (!"gene" %in% names(stats)) {
  gcol <- intersect(c("Gene", "GENE"), names(stats))[1]
  if (!is.na(gcol)) stats <- dplyr::rename(stats, gene = !!gcol)
}
stopifnot("gene" %in% names(stats))

# Without a universe file (or with fewer than 2 genes in it) the universe is the set of scored genes.
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


# ── Enrichment ────────────────────────────────────────────────────────────────

gmts <- fcs_load_gmts(gmt_dir)

enrich <- fcs_run_all(rankings, gmts, num_g = num_g, max_g = max_g,
                      perms_file = perms_file, fdr_thr = fdr_thr,
                      fdr_wilcoxon = fdr_wilcoxon, fdr_lachenbruch = fdr_lachenbruch,
                      p_perm_thr = pperm_thr)


# ── Output ────────────────────────────────────────────────────────────────────

readr::write_tsv(enrich, out_file)
