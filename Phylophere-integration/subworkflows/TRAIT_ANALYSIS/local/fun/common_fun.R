# Common libs
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
  name <- read_tsv(file, col_names = FALSE)
  return(df)
}

# Function to read TSV files and return data frames
read_csv_to_df <- function(file) {
  df <- read.csv(file, sep = "\t")
  return(df)
}
