# Common utilities
## Robust scaling
robust_scale <- function(x) {
  (x - median(x)) / IQR(x)
}

## Directory creation
createDir <- function(directory) {
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

read_tsv_to_df <- function(file) {
  df <- readr::read_tsv(file, col_names = TRUE, show_col_types = FALSE)
  return(df)
}

# Function to read CSV files and return data frames
read_csv_to_df <- function(file) {
  df <- read.csv(file, sep = ",")
  return(df)
}

# Find outliers using the IQR rule
findoutlier <- function(x) {
  q25 <- quantile(x, 0.25, na.rm = TRUE)
  q75 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- IQR(x, na.rm = TRUE)
  x < (q25 - 1.5 * iqr) | x > (q75 + 1.5 * iqr)
}

# Check if two confidence intervals overlap
ci_overlap <- function(lb1, ub1, lb2, ub2) {
  (lb1 <= ub2) & (lb2 <= ub1)
}

# Utility: check if a trait is numeric
is_numeric_trait <- function(x) {
  suppressWarnings(all(!is.na(as.numeric(x)) | is.na(x)))
}
