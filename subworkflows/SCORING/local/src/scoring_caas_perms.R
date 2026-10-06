#!/usr/bin/env Rscript
# scoring_caas_perms.R — Gene×cycle score table of the permulation null → genes×N matrices in caas_perms.rds.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  CAAS_CORE_MERGE Nextflow process (caas_permulation.nf → Rscript scoring_caas_perms.R ...)
#
# Turns the long table gene_cycle_scores.tsv (one row per gene and permulation cycle) into
# one genes×N matrix per direction and axis, the null that fcs_enrich.R (FCS p.perm),
# scoring_compute.R (gene permulation p) and 11.Scoring_report.Rmd read.
#
# Args (named flags, from task.script):
#   --gene-cycle-scores  TSV with Gene, cycle, global_asr, top_asr, bottom_asr,
#                        global_caas, top_caas, bottom_caas
#   --universe           gene universe, one gene per line (a "Gene" header is ignored),
#                        or NO_FILE (default; then the genes of the table)
#   --cycles             cycles replayed, one tag per line; N counts them even when a
#                        cycle left no row in the table (optional)
#   --output             output RDS (default caas_perms.rds)
#
# Output RDS: list(corStat_byrank = list(global_asr, top_asr, bottom_asr),
#                  caas_corStat_byrank = list(global, top, bottom),
#                  gene_stat, asr_stat); each matrix has genes in rows and cycles in
#                  columns, and 0 where a gene has no signal in a cycle.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
})

# ── Parse arguments ───────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
gcs_file      <- get_arg("--gene-cycle-scores")
universe_file <- get_arg("--universe", "NO_FILE")
out_file      <- get_arg("--output", "caas_perms.rds")
# cycles replayed (one tag per line): N counts them even when a cycle left no row
cycles_file   <- get_arg("--cycles", NULL)

stopifnot(!is.null(gcs_file), file.exists(gcs_file))

# ── Load inputs ───────────────────────────────────────────────────────────────

gcs <- read_tsv(gcs_file, show_col_types = FALSE)

# Statistic stamps: they pin the null to the observed-side formula it must match.
#   gene_stat : the gene-level CAAS aggregator (scoring_compute.R `size_adj_max`, mirrored
#               by gene_wrapper.py `_size_adj_max_null`). fcs_enrich.R (p.perm left NA)
#               and scoring_compute.R (stops) check it before using `caas_corStat_byrank`.
#   asr_stat  : the aggregator of the ASR path scores behind `corStat_byrank` (`*_asr`),
#               the 90th percentile per gene and direction (gene_wrapper.py `_q90`).
GENE_STAT <- "size_adj_max"
ASR_STAT  <- "q90"

if (nrow(gcs) == 0) {
  message("[caas_perms] no scored rows; writing empty RDS")
  saveRDS(list(corStat_byrank = list(), caas_corStat_byrank = list(),
               gene_stat = GENE_STAT, asr_stat = ASR_STAT), out_file)
  quit(status = 0)
}

# ── Cycle roster and gene universe ────────────────────────────────────────────

# Columns: the roster when given (a table cycle absent from it is an error), else the cycles of the table.
cycle_levels <- sort(unique(gcs$cycle))
if (!is.null(cycles_file)) {
  roster <- trimws(readLines(cycles_file)); roster <- roster[nzchar(roster)]
  stray <- setdiff(cycle_levels, roster)
  if (length(stray)) stop(sprintf("gene-cycle scores of %d cycle(s) not in %s: %s", length(stray), cycles_file,
                                  paste(head(stray, 5), collapse = ", ")))
  cycle_levels <- sort(unique(roster))
}
n_perms <- length(cycle_levels)

# Rows: the genes of the table plus those of the universe file when given (genes with no row stay at 0).
perm_genes <- sort(unique(gcs$Gene))
universe <- perm_genes
if (!is.null(universe_file) && universe_file != "NO_FILE" && file.exists(universe_file)) {
  u <- tryCatch(readLines(universe_file), error = function(e) character(0))
  u <- trimws(u); u <- u[nzchar(u) & u != "Gene"]
  if (length(u)) universe <- sort(unique(c(u, perm_genes)))
}

# ── Build the genes×N matrices ────────────────────────────────────────────────

# All six value columns share one lookup of (Gene, cycle) -> (row, column), computed once.
# Each matrix is then filled by a single vectorized linear-index write, so no wide
# intermediate of the genome-wide genes × cycles table is built per column. A missing
# gene/cycle pair stays 0 (no signal) and an NA value is written as 0.
gene_idx  <- match(gcs$Gene, universe)
cycle_idx <- match(gcs$cycle, cycle_levels)
valid     <- !is.na(gene_idx) & !is.na(cycle_idx)
gene_idx  <- gene_idx[valid]
cycle_idx <- cycle_idx[valid]
lin_idx   <- gene_idx + (cycle_idx - 1L) * length(universe)

value_cols <- c("global_asr", "top_asr", "bottom_asr",
                 "global_caas", "top_caas", "bottom_caas")
stopifnot(all(value_cols %in% names(gcs)))

# Matrix of one value column: genes in rows, cycles in columns.
build_matrix <- function(col) {
  mat <- matrix(0.0, nrow = length(universe), ncol = n_perms,
                dimnames = list(universe, cycle_levels))
  v <- gcs[[col]][valid]
  v[is.na(v)] <- 0.0
  mat[lin_idx] <- v
  mat
}

# ── Assemble and save ─────────────────────────────────────────────────────────

# ASR axis (path scores) and CAAS axis (size-adjusted gene scores), one matrix per direction.
corStat_byrank <- list(
  global_asr = build_matrix("global_asr"),
  top_asr    = build_matrix("top_asr"),
  bottom_asr = build_matrix("bottom_asr")
)
caas_corStat_byrank <- list(
  global = build_matrix("global_caas"),
  top    = build_matrix("top_caas"),
  bottom = build_matrix("bottom_caas")
)

saveRDS(list(corStat_byrank = corStat_byrank,
             caas_corStat_byrank = caas_corStat_byrank,
             gene_stat = GENE_STAT, asr_stat = ASR_STAT), out_file)
cat(sprintf("[caas_perms] wrote %s — %d genes × %d cycles, asr+caas nulls (%s)\n",
            out_file, length(universe), n_perms,
            paste(names(corStat_byrank), collapse = ", ")))
