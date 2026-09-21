#!/usr/bin/env Rscript
# =============================================================================
# scoring_caas_perms.R — CAAS permulation-excess null → genes×N matrices
# =============================================================================
suppressPackageStartupMessages({
  library(readr)
})

# ── minimal flag parser ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
gcs_file      <- get_arg("--gene-cycle-scores")
universe_file <- get_arg("--universe", "NO_FILE")
out_file      <- get_arg("--output", "caas_perms.rds")

stopifnot(!is.null(gcs_file), file.exists(gcs_file))

# Read the tiny gene-cycle scores file
gcs <- read_tsv(gcs_file, show_col_types = FALSE)

# Statistic stamps — pin the null to the observed-side formulas it must match.
#   gene_stat : the gene-level CAAS aggregator (scoring_compute.R `size_adj_max`
#               == gene_wrapper.py `_size_adj_max_null`); fcs_enrich.R and
#               scoring_compute.R Tier 1A both gate on this before consuming the
#               `caas_corStat_byrank` matrices (invariant break-point #11).
#   asr_stat  : the ASR path-score aggregator behind `corStat_byrank` (`*_asr`).
GENE_STAT <- "size_adj_max"
ASR_STAT  <- "q90"

if (nrow(gcs) == 0) {
  message("[caas_perms] no scored rows; writing empty RDS")
  saveRDS(list(corStat_byrank = list(), caas_corStat_byrank = list(),
               gene_stat = GENE_STAT, asr_stat = ASR_STAT), out_file)
  quit(status = 0)
}

cycle_levels <- sort(unique(gcs$cycle))
n_perms <- length(cycle_levels)

# Gene universe: cleaned_background if given, else genes present in the table
perm_genes <- sort(unique(gcs$Gene))
universe <- perm_genes
if (!is.null(universe_file) && universe_file != "NO_FILE" && file.exists(universe_file)) {
  u <- tryCatch(readLines(universe_file), error = function(e) character(0))
  u <- trimws(u); u <- u[nzchar(u) & u != "Gene"]
  if (length(u)) universe <- sort(unique(c(u, perm_genes)))
}

# ── Build genes×N matrices for all six direction columns in one pass over the
# long table (absent gene/cycle → 0 = no signal). Each build_matrix() call
# used to be its own select()+tidyr::pivot_wider() reshape of the full
# genome-wide x n_cycles table -- six full wide intermediates for what is,
# for each column, just "look up (Gene, cycle) -> (row, col) and write a
# value". gene_idx/cycle_idx/valid are computed once and reused across all
# six columns; the value assignment itself is a single vectorized linear-index
# write per matrix, so the long table (gcs) is read once instead of six times
# and no wide intermediate is ever materialized.
gene_idx  <- match(gcs$Gene, universe)
cycle_idx <- match(gcs$cycle, cycle_levels)
valid     <- !is.na(gene_idx) & !is.na(cycle_idx)
gene_idx  <- gene_idx[valid]
cycle_idx <- cycle_idx[valid]
lin_idx   <- gene_idx + (cycle_idx - 1L) * length(universe)

value_cols <- c("global_asr", "top_asr", "bottom_asr",
                 "global_caas", "top_caas", "bottom_caas")
stopifnot(all(value_cols %in% names(gcs)))

build_matrix <- function(col) {
  mat <- matrix(0.0, nrow = length(universe), ncol = n_perms,
                dimnames = list(universe, cycle_levels))
  v <- gcs[[col]][valid]
  v[is.na(v)] <- 0.0
  mat[lin_idx] <- v
  mat
}

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
