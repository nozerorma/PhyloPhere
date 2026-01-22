# Common objects for the trait analysis

## Set up variable to control command line arguments
args <- commandArgs(TRUE)

get_arg <- function(args, idx, default = NULL) {
  if (length(args) >= idx && nzchar(args[idx])) {
    return(args[idx])
  }
  default
}


# Set seed for reproducibility
seed_val <- get_arg(args, 4, "")

# ----------------------------------------
# R Markdown Setup
# ----------------------------------------

setup_rmd <- function() {
  knitr::opts_chunk$set(warning = FALSE, message = FALSE)
  knitr::opts_knit$set(root.dir =getwd()) # Set working directory to project root
  set.seed(as.integer(seed_val)) 
}
# Define working and results directories
workingDir <- getwd()
objDir <- file.path(workingDir, "obj")
resultsDir <- get_arg(args, 3, getwd())

# Load palettes, statistical functions, I/O utilities, and general utilities
source(file.path(objDir, "directories.R"))
source(file.path(objDir, "palettes.R"))
source(file.path(objDir, "stats.R"))
source(file.path(objDir, "io_utils.R"))

# Trait objects
## Trait data needs to have standardized colnames (species, family) and be comma-separated (csv format).
trait_path <- get_arg(args, 1, "")
if (!nzchar(trait_path)) {
  stop("Trait file path not provided.")
}

print(paste0("Loading trait data from: ", trait_path))

# Check if the trait file is a csv or tsv
if (endsWith(trait_path, ".csv")) {
  sep_char <- ","
} else if (endsWith(trait_path, ".tsv")) {
  sep_char <- "\t"
} else {
  sep_char <- ","
}

trait_df <- read.csv(trait_path, sep = sep_char, stringsAsFactors = FALSE)
trait_df[] <- lapply(trait_df, function(col) {
  if (is.character(col)) trimws(col) else col
})

if (!"species" %in% names(trait_df)) {
  stop("Trait file must include a 'species' column.")
}

# Phylo objects
tree_path <- get_arg(args, 2, "")
if (!nzchar(tree_path)) {
  stop("Tree file path not provided.")
}
tree <- phytools::read.newick(tree_path) %>% ape::as.phylo() # Tree object
tree_species <- tree$tip.label # Tree species

### OPTS
# Optional parameters for analysis
# Analysis metadata (trait-naive)
clade_name <- get_arg(args, 5, "clade") # Clade of interest (used for titles, file names, etc)
taxon_of_interest <- get_arg(args, 6, "family") # Taxonomic level of interest (column name in the trait file)
trait <- get_arg(args, 7, "trait") # Trait of interest (column name in the trait file)

# Check if sample size (N) column is provided
source(file.path(workingDir, "./sample_size.R"))

# Load tax_id.R to handle TAX_ID mapping if needed
source(file.path(workingDir, "./tax_id.R"))
