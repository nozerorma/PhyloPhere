#!/usr/bin/env Rscript
# =============================================================================
# scoring_caas_perms.R — CAAS permulation-excess null → genes×N matrices
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
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

# ── Build a genes×N matrix for one direction (absent gene/cycle → 0 = no signal)
build_matrix <- function(df, col) {
  m <- df %>%
    select(Gene, cycle, !!col := dplyr::all_of(col)) %>%
    tidyr::pivot_wider(names_from = cycle, values_from = !!col)
  
  mat <- matrix(0.0, nrow = length(universe), ncol = n_perms,
                dimnames = list(universe, cycle_levels))
  if (nrow(m) > 0) {
    rn <- m$Gene
    sub <- as.matrix(m[, setdiff(names(m), "Gene"), drop = FALSE])
    sub[is.na(sub)] <- 0.0
    common_cols <- intersect(colnames(sub), cycle_levels)
    common_rows <- intersect(rn, universe)
    if (length(common_rows) && length(common_cols)) {
      ri <- match(common_rows, rn)
      mat[common_rows, common_cols] <- sub[ri, common_cols, drop = FALSE]
    }
  }
  mat
}

corStat_byrank <- list(
  global_asr = build_matrix(gcs, "global_asr"),
  top_asr    = build_matrix(gcs, "top_asr"),
  bottom_asr = build_matrix(gcs, "bottom_asr")
)
caas_corStat_byrank <- list(
  global = build_matrix(gcs, "global_caas"),
  top    = build_matrix(gcs, "top_caas"),
  bottom = build_matrix(gcs, "bottom_caas")
)

saveRDS(list(corStat_byrank = corStat_byrank,
             caas_corStat_byrank = caas_corStat_byrank,
             gene_stat = GENE_STAT, asr_stat = ASR_STAT), out_file)
cat(sprintf("[caas_perms] wrote %s — %d genes × %d cycles, asr+caas nulls (%s)\n",
            out_file, length(universe), n_perms,
            paste(names(corStat_byrank), collapse = ", ")))
