#!/usr/bin/env Rscript
# scoring_caas_perms.R — The block matrices of the permulation null → the genes×N matrices of caas_perms.rds.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  CAAS_CORE_MERGE_REDUCE (validation; optional export) (caas_permulation.nf → Rscript scoring_caas_perms.R ...)
#
# Joins the per-block matrices (assemble_block_matrices.R: the genes of one block x its cycles) into one genes×N matrix
# per direction for explicit compatibility exports. Productive consumers use blocks and the roster.
#
# Args (named flags, from task.script):
#   --blocks    the block RDS files (every flag value up to the next flag)
#   --universe  gene universe, one gene per line (a "Gene" header is ignored), or NO_FILE (default; then the genes of
#               the blocks)
#   --cycles    complete cycle roster; every cycle must occur in exactly one block (optional for standalone reconstruction)
#   --validate-only  validate blocks without assembling a global matrix
#   --metadata  write the canonical block index JSON
#   --output    output RDS (default caas_perms.rds)
#
# Output RDS: list(caas_corStat_byrank = list(global, top, bottom), gene_stat); each matrix has genes in rows and cycles
# in columns, and 0 where a gene has no signal in a cycle. gene_stat names the gene statistic of the matrices
# (max only): fcs_enrich.R compares it with the statistic of the observed gene score.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
flags <- grep("^--", args)
get_args <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args) || startsWith(args[i+1], "--")) return(NULL)
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

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), "caas_block_contract.R"))
roster <- if (!is.null(cycles_file)) readLines(cycles_file) else NULL
universe <- if (universe_file != "NO_FILE") readLines(universe_file) else character()
universe <- universe[nzchar(universe) & universe != "Gene"]
index <- caas_block_index(block_files, roster, universe)
meta <- get_arg("--metadata")
if (!is.null(meta)) jsonlite::write_json(index, meta, auto_unbox=TRUE, pretty=TRUE)
if (!"--validate-only" %in% args) caas_export_single(block_files, index, out_file)
cat(sprintf("[caas_perms] validated %d blocks, %d genes, %d cycles, gene_stat %s\n",
            length(block_files), length(index$genes), length(index$cycles), index$gene_stat))
