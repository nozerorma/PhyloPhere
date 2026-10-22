#!/usr/bin/env Rscript
# assemble_block_matrices.R — The matrices of one block of cycles of the permulation null, as one RDS.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  CAAS_CORE_MERGE (subworkflows/CT/caas_permulation.nf), after merge_null_block.py
#
# merge_null_block.py writes each genes x cycles matrix as raw little-endian column-major numbers and a manifest.json;
# this script reads them back (base R only) and saves them under their names.
#
# Args (named flags, from task.script):
#   --dir     the matrices/ directory of the block (manifest.json and the raw files)
#   --output  output RDS (default block.matrices.rds)
#
# Output RDS: list(genes, cycles, gene_stat, n_units,
#                  stat  = list(global, top, bottom),                 # the gene statistic chosen with gene_stat
#                  parts = list(max = list(...), fmax = list(...), n = list(...)))
# Every matrix has the genes in rows and the cycles of the block in columns, and 0 where a gene has no value in a cycle
# (it is absent, removed, or has no scored position in the direction). `parts` holds what any gene statistic is built
# from: the largest position score of the gene, F(max) (the ECDF of the cycle's position scores at that score) and the
# number n of scored positions; F(max)^n is the default `stat`.
# =============================================================================

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
dir_in   <- get_arg("--dir")
out_file <- get_arg("--output", "block.matrices.rds")
stopifnot(!is.null(dir_in), file.exists(file.path(dir_in, "manifest.json")))

man    <- fromJSON(file.path(dir_in, "manifest.json"), simplifyVector = FALSE)
genes  <- as.character(unlist(man$genes))
cycles <- as.character(unlist(man$cycles))
n_g <- length(genes); n_c <- length(cycles)

read_matrix <- function(fname) {
  n <- n_g * n_c
  integer_file <- grepl("\\.i4$", fname)
  v <- if (n == 0) numeric(0) else
    readBin(file.path(dir_in, fname), what = if (integer_file) "integer" else "double", n = n,
            size = if (integer_file) 4L else 8L, endian = "little")
  if (length(v) != n) stop(sprintf("%s holds %d values, %d expected", fname, length(v), n))
  matrix(v, nrow = n_g, ncol = n_c, dimnames = list(genes, cycles))
}
dirs <- c("global", "top", "bottom")
named <- function(part) setNames(lapply(dirs, function(d) read_matrix(man$files[[paste0(part, ".", d)]])), dirs)

stat <- setNames(lapply(dirs, function(d) read_matrix(man$files[[man$stat[[d]]]])), dirs)
parts <- list(max = named("max"), fmax = named("fmax"), n = named("n"))

saveRDS(list(genes = genes, cycles = cycles, gene_stat = man$gene_stat, n_units = man$n_units,
             stat = stat, parts = parts), out_file)
cat(sprintf("[block] wrote %s: %d genes x %d cycles, gene_stat %s\n", out_file, n_g, n_c, man$gene_stat))
