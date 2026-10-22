#!/usr/bin/env Rscript
# scoring_caas_perms.R — The block matrices of the permulation null → the genes×N matrices of caas_perms.rds.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  CAAS_CORE_MERGE_REDUCE Nextflow process (caas_permulation.nf → Rscript scoring_caas_perms.R ...)
#
# Joins the per-block matrices (assemble_block_matrices.R: the genes of one block x its cycles) into one genes×N matrix
# per direction, the null that fcs_enrich.R (FCS p.perm), scoring_compute.R and 11.Scoring_report.Rmd read.
#
# Args (named flags, from task.script):
#   --blocks    the block RDS files (every flag value up to the next flag)
#   --universe  gene universe, one gene per line (a "Gene" header is ignored), or NO_FILE (default; then the genes of
#               the blocks)
#   --cycles    cycles replayed, one tag per line; N counts them even when a cycle left no row (optional)
#   --output    output RDS (default caas_perms.rds)
#
# Output RDS: list(caas_corStat_byrank = list(global, top, bottom), gene_stat); each matrix has genes in rows and cycles
# in columns, and 0 where a gene has no signal in a cycle. gene_stat names the gene statistic of the matrices
# (size_adj_max, fmax or max): fcs_enrich.R compares it with the statistic of the observed gene score.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
flags <- grep("^--", args)
get_args <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) return(NULL)
  j <- flags[flags > i]
  args[(i + 1):(if (length(j)) j[1] - 1 else length(args))]
}
get_arg <- function(flag, default = NULL) {
  v <- get_args(flag)
  if (is.null(v) || !length(v)) default else v[1]
}
block_files   <- get_args("--blocks")
universe_file <- get_arg("--universe", "NO_FILE")
out_file      <- get_arg("--output", "caas_perms.rds")
cycles_file   <- get_arg("--cycles", NULL)
stopifnot(length(block_files) > 0, all(file.exists(block_files)))

blocks <- lapply(block_files, readRDS)
stats <- unique(vapply(blocks, function(b) as.character(b$gene_stat), ""))
if (length(stats) != 1) stop(sprintf("the blocks hold different gene statistics: %s", paste(stats, collapse = ", ")))
GENE_STAT <- stats
dirs <- c("global", "top", "bottom")

if (sum(vapply(blocks, function(b) as.numeric(b$n_units), 0)) == 0) {
  message("[caas_perms] no scored rows; writing empty RDS")
  saveRDS(list(caas_corStat_byrank = list(), gene_stat = GENE_STAT), out_file)
  quit(status = 0)
}

# ── Cycle roster and gene universe ────────────────────────────────────────────

# Columns: the roster when given (a block cycle absent from it is an error), else the cycles of the blocks.
cycle_levels <- sort(unique(unlist(lapply(blocks, function(b) b$cycles))))
if (!is.null(cycles_file)) {
  roster <- trimws(readLines(cycles_file)); roster <- roster[nzchar(roster)]
  stray <- setdiff(cycle_levels, roster)
  if (length(stray)) stop(sprintf("gene-cycle scores of %d cycle(s) not in %s: %s", length(stray), cycles_file,
                                  paste(head(stray, 5), collapse = ", ")))
  cycle_levels <- sort(unique(roster))
}
n_perms <- length(cycle_levels)

# Rows: the genes with a unit in some block plus those of the universe file when given (genes with no row stay at 0).
# A block holds the universe and its own genes, so the union of the blocks is both.
universe <- sort(unique(unlist(lapply(blocks, function(b) b$genes))))
if (!is.null(universe_file) && universe_file != "NO_FILE" && file.exists(universe_file)) {
  u <- tryCatch(readLines(universe_file), error = function(e) character(0))
  u <- trimws(u); u <- u[nzchar(u) & u != "Gene"]
  if (length(u)) universe <- sort(unique(c(u, universe)))
}

# ── Join the blocks ───────────────────────────────────────────────────────────

caas_corStat_byrank <- lapply(dirs, function(d) {
  mat <- matrix(0.0, nrow = length(universe), ncol = n_perms, dimnames = list(universe, cycle_levels))
  for (b in blocks) {
    if (!length(b$genes) || !length(b$cycles)) next
    mat[match(b$genes, universe), match(b$cycles, cycle_levels)] <- b$stat[[d]]
  }
  mat
})
names(caas_corStat_byrank) <- dirs

saveRDS(list(caas_corStat_byrank = caas_corStat_byrank, gene_stat = GENE_STAT), out_file)
cat(sprintf("[caas_perms] wrote %s — %d genes × %d cycles, %d blocks, gene_stat %s\n",
            out_file, length(universe), n_perms, length(blocks), GENE_STAT))
