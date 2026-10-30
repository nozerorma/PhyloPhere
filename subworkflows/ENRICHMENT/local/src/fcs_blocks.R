#!/usr/bin/env Rscript
# fcs_blocks.R — steps of the CAAS FCS on the calibration of SCORING (-log10 p slices of SCORING_EVIDENCE). One Rscript call per step; see caas_blocks.R.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
#   observed    --obs-dir <dir> --universe <file> [--metadata <calibration metadata.json>] --rankings <global,top,bottom> --gmt-dir <dir> --num-g --max-g --out <tsv>
#               observed FCS rows of one GMT batch from the obs_<rk>_<k>.tsv files (p.perm not set)
#   null_block  --observed <tsv> --slices <dir> --universe <file> --block <id> --batch-size <G> --rank <rk> --alt <greater|two.sided>
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
if (!identical(get_arg("--gene-stat", "max"), "max")) stop("Only CAAS gene_stat=max is supported")
step <- args[1]
src_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
suppressPackageStartupMessages({ library(readr); library(dplyr) })
source(file.path(src_dir, "fcs_enrich.R"))
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
read_observed <- function(path) {
  if (grepl('[.]rds$', path)) readRDS(path) else read_tsv(path, show_col_types=FALSE)
}
write_observed <- function(tab, path) {
  if (grepl('[.]rds$', path)) saveRDS(tab, path) else write_tsv(tab, path)
}
read_universe <- function(path) {
  u <- unique(trimws(readLines(path))); u[nzchar(u)]
}

rer_rankings <- function() {
  tab <- read_tsv(get_arg('--obs'), show_col_types = FALSE)
  fcs_check_gene_ids(tab$gene, nrow(tab), 'RER observed')
  u <- read_universe(get_arg('--universe'))
  cols <- grep('^score_', names(tab), value = TRUE)
  setNames(lapply(cols, function(c) fcs_build_vals(setNames(tab[[c]], tab$gene), u, score_source='rer')),
           sub('^score_', '', cols))
}
empty_counts <- function() tibble::tibble(database=character(), pathway=character(), count=numeric(), valid=numeric(), ranking=character())
switch(step,
  rer_observed = {
    obs <- fcs_observed(rer_rankings(), fcs_load_gmts(get_arg('--gmt-dir')),
                        num_g=as.integer(get_arg('--num-g',10)), max_g=as.integer(get_arg('--max-g',0)), score_source='rer')
    write_observed(if(is.null(obs$enrich)) empty_fcs_table() else obs$enrich, get_arg('--out'))
  },
  rer_null = {
    m <- readRDS(get_arg('--matrices'))$corStat
    counts <- empty_counts()
    if (!is.null(m)) {
      m <- as.matrix(m)
      fcs_check_rer_gene_ids(rownames(m), nrow(m), 'RER block')
      obs <- read_observed(get_arg('--observed'))
      rankings <- rer_rankings()
      gmts <- fcs_load_gmts(get_arg('--gmt-dir'))
      counts <- bind_rows(lapply(names(rankings), function(rk) {
        # RER excludes unavailable scores separately in each null column.
        # Genes unavailable in the observed ranking can be evaluable in the null.
        aligned <- m
        rows <- obs[obs$ranking == rk,,drop=FALSE]
        if (!nrow(rows)) return(NULL)
        fcs_null_counts_by_matrices(rk, rows, list(aligned), gmts, 'two.sided',
                                    num_g=as.integer(get_arg('--num-g',10)), max_g=as.integer(get_arg('--max-g',0)))
      }))
      if (!nrow(counts)) counts <- empty_counts()
    }
    write_tsv(counts, get_arg('--out'))
  },
  observed = {
    dir <- get_arg("--obs-dir"); universe <- read_universe(get_arg("--universe"))
    # The calibration was made on one universe: its slices cover exactly those genes in that order. A different universe
    # (another file, order or content) would mix statistics of two gene families, so it is rejected here.
    meta_path <- get_arg("--metadata", "")
    if (nzchar(meta_path)) {
      meta <- jsonlite::read_json(meta_path, simplifyVector = TRUE)
      if (is.null(meta$universe_md5) || !identical(as.character(meta$universe_md5), unname(as.character(tools::md5sum(get_arg("--universe"))))))
        stop("FCS universe differs from the universe of the calibration (universe_md5 of metadata.json): recalibrate on this universe")
      if (!identical(as.integer(meta$genes), length(universe)))
        stop("Calibration gene count (", meta$genes, ") differs from the genes of the FCS universe (", length(universe), ")")
    }
    rks <- strsplit(get_arg("--rankings", "global"), ",")[[1]]
    rankings <- list()
    for (rk in rks) {
      files <- list.files(dir, pattern = sprintf("^obs_%s_[0-9]{4}\\.tsv$", rk), full.names = TRUE)
      if (length(files) == 0) next
      tab <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE))
      fcs_check_gene_ids(tab$gene, nrow(tab), "observed gene-p batches")
      if (!setequal(tab$gene, universe)) stop("Observed gene-p batches do not cover the FCS universe exactly (ranking ", rk, ")")
      vals <- setNames(rep(NA_real_, length(universe)), universe)
      vals[tab$gene] <- tab$nlog10_obs
      rankings[[rk]] <- fcs_build_vals(vals, universe)
    }
    gmts <- fcs_load_gmts(get_arg("--gmt-dir"))
    obs <- fcs_observed(rankings, gmts, num_g = as.integer(get_arg("--num-g", 10)), max_g = as.integer(get_arg("--max-g", 0)))
    out <- if (is.null(obs$enrich)) empty_fcs_table() else obs$enrich
    write_observed(out, get_arg("--out"))
  },
  null_block = {
    obs <- read_observed(get_arg("--observed"))
    rk <- get_arg("--rank", "global")
    obs_rk <- obs[obs$ranking == rk, , drop = FALSE]
    universe <- read_universe(get_arg("--universe"))
    # The gene batches are those of the universe as read_universe() returns it, as in the calibration.
    batch_size <- as.integer(get_arg("--batch-size"))
    if (is.na(batch_size)) stop("null_block needs --batch-size")
    batches <- seq_along(caas_batches(universe, batch_size))
    mats <- list(caas_assemble_block(batches, rk, get_arg("--block"), get_arg("--slices"), universe))
    gmts <- fcs_load_gmts(get_arg("--gmt-dir"))
    counts <- if (nrow(obs_rk) == 0) NULL else
      fcs_null_counts_by_matrices(rk, obs_rk, mats, gmts, get_arg("--alt", "greater"),
                                  num_g = as.integer(get_arg("--num-g", 10)), max_g = as.integer(get_arg("--max-g", 0)))
    write_tsv(if (is.null(counts)) tibble::tibble(database = character(), pathway = character(), count = numeric(),
                                                  valid = numeric(), ranking = character()) else counts, get_arg("--out"))
  },
  close = {
    obs <- read_observed(get_arg("--observed"))
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
