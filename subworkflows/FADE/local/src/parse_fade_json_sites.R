#!/usr/bin/env Rscript
# parse_fade_json_sites.R — Site-level FADE Bayes factors from the raw *.FADE.json files.
# PhyloPhere | subworkflows/FADE/local/src/
# =============================================================================
# Called by:  FADE_JSON_TO_CSV Nextflow process (fade_json_to_csv.nf → Rscript parse_fade_json_sites.R ...)
#
# HyPhy FADE reports one Bayes factor per (site, target amino acid). The gene-level
# FADE report keeps one row per gene, so this script reads the JSON files again and
# writes one row per site whose maximum BF over the 20 amino acids reaches --bf_thr:
# the position-keyed ((Gene, Position)) FADE evidence layer of POSENRICH.
#
# Output columns: gene, position (0-based, the coordinate system of the CAAS Position
# column), max_bf, target_aa (the amino acid with the highest BF). A header-only
# CSV is written when there are no JSON files or no site reaches the threshold.
#
# Args (named flags, from task.script):
#   --json_dir   directory with the <gene>.<direction>.FADE.json files
#   --direction  top | bottom
#   --bf_thr     minimum max BF to keep a site (default 100)
#   --n_cores    parallel workers (default 4; fork-based, serial on non-unix)
#   --out        output CSV (default fade_sites_<direction>.csv)
# =============================================================================

# ── Dependencies ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(jsonlite)
  library(parallel)
})


# ── Arguments ─────────────────────────────────────────────────────────────────

# Value after `flag` on the command line, or `default` when the flag is absent.
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

# Order of the amino acids (keys of the MLE content in the JSON)
AA_LETTERS <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]


# ── Input files ───────────────────────────────────────────────────────────────

json_files <- list.files(json_dir, pattern = "\\.FADE\\.json$", full.names = TRUE)
cat(sprintf("[parse_fade_json_sites] direction=%s | %d JSON files | bf_thr=%.1f\n",
            direction, length(json_files), bf_thr))

if (length(json_files) == 0) {
  # Header-only CSV: the declared output exists and build_position_gmt.py reads an
  # empty FADE layer.
  writeLines("gene,position,max_bf,target_aa", out_file)
  cat("[parse_fade_json_sites] no JSON files found — wrote empty output\n")
  quit(save = "no", status = 0)
}


# ── Per-gene parsing ──────────────────────────────────────────────────────────

# Significant sites of one JSON file as a data frame, NULL when the file has no MLE
# content, no site reaches --bf_thr, or the file cannot be parsed (reported as a message).
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
      # content[[aa]] is keyed by partition (a single key in these runs); the site
      # index is the position within that partition's list, not the key name.
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
      max_bf = max_bf[sig_idx],                  # as in the CAAS Position column
      target_aa = AA_LETTERS[top_aa_idx],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("[parse_fade_json_sites] failed to parse %s: %s", basename(path), conditionMessage(e)))
    NULL
  })
}


# ── Parse all files and write ─────────────────────────────────────────────────

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
