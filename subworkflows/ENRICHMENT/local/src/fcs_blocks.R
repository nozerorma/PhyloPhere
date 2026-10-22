#!/usr/bin/env Rscript
# fcs_blocks.R — steps of the CAAS FCS on blocks of the null (WP6). One Rscript call per step; see caas_blocks.R.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
#   transpose   --matrices <blk>/matrices.rds --universe <file> --batch <G> --rank <global|top|bottom> --block <id> --out <dir>
#               one slice per batch of the universe (rows of the batch, cycles of the block)
#   gene_p      --batch <k> --rank <rk> --blocks <id,id,...> --slices <dir> --obs <table.tsv> --obs-col <col> --out <dir>
#               GENE_P of batch k over every block; writes the -log10 p slices and obs_<rk>_<k>.tsv (gene, nlog10_obs)
#   observed    --obs-dir <dir> --universe <file> --rankings <global,top,bottom> --gmt-dir <dir> --num-g --max-g --out <tsv>
#               observed FCS rows of one GMT batch from the obs_<rk>_<k>.tsv files (p.perm not set)
#   null_block  --observed <tsv> --slices <dir> --universe <file> --block <id> --batches <k,k,...> --rank <rk> --alt <greater|two.sided>
#               --gmt-dir <dir> --num-g --max-g --out <tsv>
#               permulation counts of one block for the observed rows of one GMT batch
#   close       --observed <tsv> --counts <tsv,tsv,...> --fdr-thr --pperm-thr --analytic-levels --out <tsv>
#               sums the counts, sets p.perm and q.perm, classifies: the same table fcs_run_all returns
# The steps read and write only the files named above.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i[length(i)] == length(args)) default else args[i[length(i)] + 1]
}
step <- args[1]
src_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(readr); library(dplyr) })
if (!exists("fcs_percentile_flags", mode = "function")) source(file.path(src_dir, "percentile_flags.R"))
source(file.path(src_dir, "fcs_enrich.R"))
source(file.path(src_dir, "gene_p.R"))
source(file.path(src_dir, "caas_blocks.R"))

# Counts of one block: the types are fixed so that an empty or all-NA file still reads as numbers.
read_counts <- function(path) {
  read_tsv(path, col_types = cols(ranking = col_character(), database = col_character(), pathway = col_character(),
                                  count = col_double(), valid = col_double()), show_col_types = FALSE)
}
# Empty FCS table with the schema of fcs_run_all (a batch with no set writes this, as the unblocked code does).
empty_fcs_table <- function() {
  tibble::tibble(
    ranking = character(), database = character(), pathway = character(),
    stat = numeric(), pval = numeric(), p.adj = numeric(), p.perm = numeric(), q.perm = numeric(),
    num.genes = numeric(), gene.vals = character(),
    sig_perm = logical(), direction = character(), redlich_level = numeric())
}
read_universe <- function(path) {
  u <- unique(trimws(readLines(path))); u[nzchar(u)]
}
obs_col_name <- function(rk) c(global = "gene_caas_score", top = "gene_caas_score_top", bottom = "gene_caas_score_bottom")[[rk]]

switch(step,
  transpose = {
    m <- readRDS(get_arg("--matrices"))
    rk <- get_arg("--rank", "global")
    universe <- read_universe(get_arg("--universe"))
    batches <- caas_batches(universe, as.integer(get_arg("--batch")))
    caas_transpose_block(m$stat[[rk]], batches, rk, get_arg("--block"), get_arg("--out"))
  },
  gene_p = {
    k <- as.integer(get_arg("--batch")); rk <- get_arg("--rank", "global")
    blocks <- strsplit(get_arg("--blocks"), ",")[[1]]
    obs_tab <- read_tsv(get_arg("--obs"), show_col_types = FALSE)
    col <- get_arg("--obs-col", obs_col_name(rk))
    obs <- setNames(as.numeric(obs_tab[[col]]), obs_tab$gene)
    out <- get_arg("--out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
    slices <- get_arg("--slices")                         # the -log10 p slices are written beside the input slices
    obs_k <- caas_gene_p_batch(k, rk, blocks, slices, obs)
    write_tsv(tibble::tibble(gene = names(obs_k), nlog10_obs = as.numeric(obs_k)),
              file.path(out, sprintf("obs_%s_%04d.tsv", rk, k)))
  },
  observed = {
    dir <- get_arg("--obs-dir"); universe <- read_universe(get_arg("--universe"))
    rks <- strsplit(get_arg("--rankings", "global"), ",")[[1]]
    rankings <- list()
    for (rk in rks) {
      files <- list.files(dir, pattern = sprintf("^obs_%s_[0-9]{4}\\.tsv$", rk), full.names = TRUE)
      if (length(files) == 0) next
      tab <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE))
      fcs_check_gene_ids(tab$gene, nrow(tab), "observed gene-p batches")
      vals <- setNames(rep(NA_real_, length(universe)), universe)
      vals[tab$gene] <- tab$nlog10_obs
      rankings[[rk]] <- fcs_build_vals(vals, universe)
    }
    gmts <- fcs_load_gmts(get_arg("--gmt-dir"))
    obs <- fcs_observed(rankings, gmts, num_g = as.integer(get_arg("--num-g", 10)), max_g = as.integer(get_arg("--max-g", 0)))
    out <- if (is.null(obs$enrich)) empty_fcs_table() else obs$enrich
    write_tsv(out, get_arg("--out"))
  },
  null_block = {
    obs <- read_tsv(get_arg("--observed"), show_col_types = FALSE)
    rk <- get_arg("--rank", "global")
    obs_rk <- obs[obs$ranking == rk, , drop = FALSE]
    universe <- read_universe(get_arg("--universe"))
    batches <- as.integer(strsplit(get_arg("--batches"), ",")[[1]])
    mats <- list(caas_assemble_block(batches, rk, get_arg("--block"), get_arg("--slices"), universe))
    gmts <- fcs_load_gmts(get_arg("--gmt-dir"))
    counts <- if (nrow(obs_rk) == 0) NULL else
      fcs_null_counts_by_matrices(rk, obs_rk, mats, gmts, get_arg("--alt", "greater"),
                                  num_g = as.integer(get_arg("--num-g", 10)), max_g = as.integer(get_arg("--max-g", 0)))
    write_tsv(if (is.null(counts)) tibble::tibble(database = character(), pathway = character(), count = numeric(),
                                                  valid = numeric(), ranking = character()) else counts, get_arg("--out"))
  },
  close = {
    obs <- read_tsv(get_arg("--observed"), show_col_types = FALSE)
    files <- strsplit(get_arg("--counts", ""), ",")[[1]]
    files <- files[nzchar(files)]
    # the counts of the blocks add up: one row per (ranking, database, pathway) before the close
    counts <- if (length(files) == 0) NULL else
      bind_rows(lapply(files, read_counts)) %>%
      group_by(ranking, database, pathway) %>%
      summarise(count = sum(count), valid = sum(valid), .groups = "drop")
    lv <- as.numeric(strsplit(get_arg("--analytic-levels", "0.05,0.1,0.15"), ",")[[1]])
    if (nrow(obs) == 0) { write_tsv(empty_fcs_table(), get_arg("--out")); quit(save = "no", status = 0) }
    enrich <- fcs_close_counts(as.data.frame(obs), counts)
    enrich <- fcs_classify_significance(enrich, fdr_thr = as.numeric(get_arg("--fdr-thr", 0.15)),
                                        p_perm_thr = as.numeric(get_arg("--pperm-thr", 0.025)), analytic_levels = lv)
    enrich <- dplyr::relocate(enrich, ranking, database, pathway, stat, pval, p.adj, p.perm, q.perm,
                              num.genes, gene.vals, sig_perm, direction, redlich_level)
    write_tsv(enrich, get_arg("--out"))
  },
  stop("unknown step: ", step))
