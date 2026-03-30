# ----------------------------------------
# Input/Output Utilities
# ----------------------------------------

# Load required library
library(readr)
library(dplyr)
library(tidyr)

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Create a directory if it doesn't exist
createDir <- function(directory) {
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    debug_log("createDir: created %s", directory)
  } else {
    debug_log("createDir: exists %s", directory)
  }
}

# Read CSV file into a data frame
read_csv_to_df <- function(file) {
  debug_log("read_csv_to_df: %s", file)
  df <- read.csv(file, sep = ",")
  debug_log("read_csv_to_df: rows = %d, cols = %d", nrow(df), ncol(df))
  return(df)
}

# Read TSV file into a data frame
read_tsv_to_df <- function(file) {
  debug_log("read_tsv_to_df: %s", file)
  df <- read_tsv(file, col_names = TRUE)
  debug_log("read_tsv_to_df: rows = %d, cols = %d", nrow(df), ncol(df))
  return(df)
}
