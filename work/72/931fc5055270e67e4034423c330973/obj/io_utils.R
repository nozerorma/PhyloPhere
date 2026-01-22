# ----------------------------------------
# Input/Output Utilities
# ----------------------------------------

# Load required library
library(readr)
library(dplyr)
library(tidyr)


# Create a directory if it doesn't exist
createDir <- function(directory) {
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

# Read CSV file into a data frame
read_csv_to_df <- function(file) {
  df <- read.csv(file, sep = ",")
  return(df)
}

# Read TSV file into a data frame
read_tsv_to_df <- function(file) {
  df <- read_tsv(file, col_names = TRUE)
  return(df)
}