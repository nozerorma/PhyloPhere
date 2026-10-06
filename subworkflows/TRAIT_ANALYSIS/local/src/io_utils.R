# io_utils.R — Directory and table-reading helpers for the trait-analysis reports.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Defines createDir(), read_csv_to_df() and read_tsv_to_df(). Each logs through
# debug_log() so the report log records every path touched.
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)

# Fallback logger, used only when commons.R has not defined debug_log().
if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# ── Directories ───────────────────────────────────────────────────────────────

# Create `directory` (and missing parents); an existing one is left untouched.
createDir <- function(directory) {
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    debug_log("createDir: created %s", directory)
  } else {
    debug_log("createDir: exists %s", directory)
  }
}

# ── Table readers ─────────────────────────────────────────────────────────────

# Comma-separated file with header, read with base R (column types guessed by read.csv).
read_csv_to_df <- function(file) {
  debug_log("read_csv_to_df: %s", file)
  df <- read.csv(file, sep = ",")
  debug_log("read_csv_to_df: rows = %d, cols = %d", nrow(df), ncol(df))
  return(df)
}

# Tab-separated file with header, read with readr (returns a tibble).
read_tsv_to_df <- function(file) {
  debug_log("read_tsv_to_df: %s", file)
  df <- read_tsv(file, col_names = TRUE)
  debug_log("read_tsv_to_df: rows = %d, cols = %d", nrow(df), ncol(df))
  return(df)
}
