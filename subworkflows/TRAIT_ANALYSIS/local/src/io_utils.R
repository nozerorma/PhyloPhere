# io_utils.R — Directory and table-reading helpers for the trait-analysis reports.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Defines createDir(), read_trait_table() and write_trait_archive(). Paths log through
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

# ── Trait tables ───────────────────────────────────────────────────────────────

# Keep names verbatim: configured columns may contain spaces or punctuation.
read_trait_table <- function(path) {
  sep <- if (grepl("\\.(tsv|txt|tab)$", tolower(path))) "\t" else ","
  df <- read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                   stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE, colClasses = "character",
                   na.strings = c("NA", ""))
  if (any(!nzchar(names(df))) || anyDuplicated(names(df))) {
    stop("Trait table has empty or duplicate column names: ", path)
  }
  df[] <- lapply(df, function(col) if (is.character(col)) trimws(col) else col)
  df
}

write_trait_archive <- function(input, output) {
  write.table(read_trait_table(input), output, sep = "\t", row.names = FALSE,
              quote = TRUE, na = "NA")
}
