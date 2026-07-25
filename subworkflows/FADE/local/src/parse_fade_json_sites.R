#!/usr/bin/env Rscript
# =============================================================================
# parse_fade_json_sites.R — per-site FADE Bayes Factors, direct from raw JSON
# =============================================================================
# HyPhy FADE reports one Bayes Factor per (site, target amino acid); the
# gene-level FADE_REPORT (6.FADE_report.Rmd) collapses this to a single
# per-gene row and no longer keeps the underlying site-level table (dropped as
# a memory optimisation — nothing downstream used to consume it). This script
# regenerates the piece posenrich needs: one row per site whose max BF (across
# the 20 target AAs) clears the classic significance bar, giving FADE a
# position-level ((Gene,Position)-keyed) evidence layer analogous to UCR
# core/flank and FUBAR positive/purifying.
#
# Validated this session against a real 16,130-file production run: 1,195
# significant (gene,site) rows for the `top` direction, matching the pipeline's
# own gene-level fade_summary_top.tsv aggregate (sum(n_sig_site_aa)) exactly.
#
# Usage:
#   Rscript parse_fade_json_sites.R \
#     --json_dir   <dir with *.FADE.json> \
#     --direction  top|bottom \
#     --bf_thr     100 \
#     --n_cores    8 \
#     --out        fade_sites_top.csv
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(parallel)
})

parse_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}

json_dir  <- parse_arg("--json_dir")
direction <- parse_arg("--direction")
bf_thr    <- as.numeric(parse_arg("--bf_thr", "100"))
n_cores   <- as.integer(parse_arg("--n_cores", "4"))
out_file  <- parse_arg("--out", sprintf("fade_sites_%s.csv", direction))

stopifnot(!is.null(json_dir), !is.null(direction))

AA_LETTERS <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

json_files <- list.files(json_dir, pattern = "\\.FADE\\.json$", full.names = TRUE)
cat(sprintf("[parse_fade_json_sites] direction=%s | %d JSON files | bf_thr=%.1f\n",
            direction, length(json_files), bf_thr))

if (length(json_files) == 0) {
  # Write header-only CSV so the Nextflow output declaration is satisfied and
  # downstream (build_position_gmt.py) sees an honestly-empty FADE layer
  # rather than a missing file.
  writeLines("gene,position,max_bf,target_aa", out_file)
  cat("[parse_fade_json_sites] no JSON files found — wrote empty output\n")
  quit(save = "no", status = 0)
}

parse_one <- function(path) {
  gene_id <- sub("\\.(top|bottom)\\.FADE\\.json$", "", basename(path))
  tryCatch({
    js <- fromJSON(path, simplifyVector = FALSE)
    if (is.null(js[["MLE"]]) || is.null(js[["MLE"]][["content"]])) return(NULL)
    hdrs <- js[["MLE"]][["headers"]]
    bf_idx <- which(sapply(hdrs, function(h) grepl("BayesFactor", h[[1]], ignore.case = TRUE)))
    if (length(bf_idx) == 0) bf_idx <- 4L
    content <- js[["MLE"]][["content"]]

    bf_matrix <- NULL
    for (aa in AA_LETTERS) {
      aa_data <- content[[aa]]
      if (is.null(aa_data) || length(aa_data) == 0) next
      # content[[aa]] is keyed by partition (in practice always one key, "0");
      # the per-gene site index is the POSITION within that partition's list,
      # not the key name.
      site_key <- names(aa_data)[1]
      all_site_rows <- aa_data[[site_key]]
      if (is.null(all_site_rows) || length(all_site_rows) == 0) next
      bf_vals <- vapply(all_site_rows, function(rd) {
        rn <- as.numeric(unlist(rd))
        if (length(rn) >= bf_idx) rn[[bf_idx]] else NA_real_
      }, numeric(1))
      if (is.null(bf_matrix)) {
        bf_matrix <- matrix(NA_real_, nrow = length(bf_vals), ncol = length(AA_LETTERS),
                            dimnames = list(NULL, AA_LETTERS))
      }
      bf_matrix[, aa] <- bf_vals
    }
    if (is.null(bf_matrix)) return(NULL)

    max_bf <- apply(bf_matrix, 1, max, na.rm = TRUE)
    sig_idx <- which(is.finite(max_bf) & max_bf >= bf_thr)
    if (length(sig_idx) == 0) return(NULL)
    top_aa_idx <- apply(bf_matrix[sig_idx, , drop = FALSE], 1, which.max)
    data.frame(
      gene = gene_id, position = sig_idx - 1L,   # 1-based site -> 0-based position,
      max_bf = max_bf[sig_idx],                  # matching CAAS/CT's coordinate system
      target_aa = AA_LETTERS[top_aa_idx],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("[parse_fade_json_sites] failed to parse %s: %s", basename(path), conditionMessage(e)))
    NULL
  })
}

t0 <- Sys.time()
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  results <- mclapply(json_files, parse_one, mc.cores = n_cores)
} else {
  results <- lapply(json_files, parse_one)
}
cat(sprintf("[parse_fade_json_sites] parsed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

results <- results[!vapply(results, is.null, logical(1))]
out_df <- if (length(results)) do.call(rbind, results) else
  data.frame(gene = character(), position = integer(), max_bf = numeric(), target_aa = character())

write.csv(out_df, out_file, row.names = FALSE)
cat(sprintf("[parse_fade_json_sites] wrote %s: %d significant (gene,position) rows across %d genes\n",
            out_file, nrow(out_df), length(unique(out_df$gene))))
