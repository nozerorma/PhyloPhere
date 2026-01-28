# Common objects for the trait analysis

if (is.null(getOption("phylo_debug"))) {
  options(phylo_debug = TRUE)
}
if (!exists("phylo_debug_log", envir = .GlobalEnv)) {
  phylo_debug_log <- character()
}
debug_log <- function(...) {
  if (isTRUE(getOption("phylo_debug", FALSE))) {
    msg <- sprintf(...)
    if (exists("phylo_debug_log", envir = .GlobalEnv)) {
      phylo_debug_log <<- c(phylo_debug_log, msg)
    }
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

## Set up variable to control command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "--args") {
  args <- args[-1]
}
debug_log("args = %s", ifelse(length(args) > 0, paste(args, collapse = " | "), "<none>"))

get_arg <- function(args, idx, default = NULL) {
  if (length(args) >= idx && nzchar(args[idx])) {
    return(args[idx])
  }
  default
}


# Set seed for reproducibility
seed_val <- get_arg(args, 4, "")
debug_log("seed_val = %s", ifelse(nzchar(seed_val), seed_val, "<empty>"))

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
debug_log("workingDir = %s", workingDir)
debug_log("objDir = %s", objDir)
debug_log("resultsDir = %s", resultsDir)

# Load palettes, statistical functions, I/O utilities, and general utilities
source(file.path(objDir, "io_utils.R"))
source(file.path(objDir, "directories.R"))
source(file.path(objDir, "palettes.R"))
source(file.path(objDir, "stats.R"))
source(file.path(objDir, "plotting_fun.R"))

# ----------------------------------------

# Trait objects
## Trait data needs to have standardized colnames (species, family) and be comma-separated (csv format).
trait_path <- get_arg(args, 1, "")
if (!nzchar(trait_path)) {
  stop("Trait file path not provided.")
}

print(paste0("Loading trait data from: ", trait_path))
debug_log("trait_path = %s", trait_path)

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
debug_log("trait_df rows = %d, cols = %d", nrow(trait_df), ncol(trait_df))
debug_log("trait_df columns: %s", paste(names(trait_df), collapse = ", "))

if (!"species" %in% names(trait_df)) {
  stop("Trait file must include a 'species' column.")
}
debug_log("trait_df species unique = %d", length(unique(trait_df$species)))

# Phylo objects
tree_path <- get_arg(args, 2, "")
if (!nzchar(tree_path)) {
  stop("Tree file path not provided.")
}
debug_log("tree_path = %s", tree_path)
tree <- phytools::read.newick(tree_path) %>% ape::as.phylo() # Tree object
tree_species <- tree$tip.label # Tree species
debug_log("tree tips = %d, nodes = %d", length(tree$tip.label), tree$Nnode)

### OPTS
# Optional parameters for analysis
# Analysis metadata (trait-naive)
clade_name <- get_arg(args, 5, "clade") # Clade of interest (used for titles, file names, etc)
taxon_of_interest <- get_arg(args, 6, "family") # Taxonomic level of interest (column name in the trait file)
trait <- get_arg(args, 7, "trait") # Trait of interest (column name in the trait file)
debug_log("clade_name = %s", clade_name)
debug_log("taxon_of_interest = %s", taxon_of_interest)
debug_log("trait = %s", trait)

# Check if sample size (N) column is provided
source(file.path(objDir, "sample_size.R"))

# Load phylo.R to handle phylogenetic tree processing and tax_id mapping if needed
source(file.path(objDir, "phylo.R"))

# Is there a second trait of interest? (Optional, for fan tree visualization. Related to co-visualization)
secondary_trait <- get_arg(args, 11, "")
debug_log("secondary_trait = %s", ifelse(nzchar(secondary_trait), secondary_trait, "<none>"))
has.secondary <- FALSE
if (nzchar(secondary_trait) && secondary_trait %in% names(trait_df)) {
  has.secondary <- TRUE
  debug_log("has.secondary = TRUE, missing = %d", sum(is.na(trait_df[[secondary_trait]])))
} else {
  message("No valid secondary trait provided; proceeding without it.")
  debug_log("has.secondary = FALSE")
}

# Is there any other trait you'd like to include in the fan tree visualization? (Optional, for branch coloring)
branch_trait <- get_arg(args, 12, "")
debug_log("branch_trait = %s", ifelse(nzchar(branch_trait), branch_trait, "<none>"))
has.branch <- FALSE
if (nzchar(branch_trait) && branch_trait %in% names(trait_df)) {
  has.branch <- TRUE
  debug_log("has.branch = TRUE, missing = %d", sum(is.na(trait_df[[branch_trait]])))
} else {
  message("No valid branch trait provided; proceeding without it.")
  debug_log("has.branch = FALSE")
}
