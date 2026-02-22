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

# ----------------------------------------
# Parameter Access via YAML params
# ----------------------------------------
# NOTE: get_arg() function is DEPRECATED - use params$ directly in Rmd files
# This function is kept only for backward compatibility but should NOT be used in new code
#
# Modern approach (use this):
#   trait_path <- params$trait_file
#   seed_val <- params$seed
#
# Old approach (deprecated):
#   trait_path <- get_arg(args, 1, "")
#
get_arg <- function(args, idx, default = NULL) {
  warning("get_arg() is deprecated. Use params$ from YAML header instead.")
  if (length(args) >= idx && nzchar(args[idx])) {
    return(args[idx])
  }
  default
}

# Access parameters from YAML header (passed via rmarkdown::render params list)
# These variables are set from params in the calling Rmd file's setup chunk
trait_path <- if (exists("params") && !is.null(params$trait_file)) params$trait_file else stop("trait_file parameter required")
tree_path <- if (exists("params") && !is.null(params$tree_file)) params$tree_file else stop("tree_file parameter required")
resultsDir <- if (exists("params")) params$output_dir else getwd()
seed_val <- if (exists("params")) params$seed else ""

debug_log("trait_path = %s", trait_path)
debug_log("tree_path = %s", tree_path)
debug_log("resultsDir = %s", resultsDir)
debug_log("seed_val = %s", ifelse(nzchar(seed_val), seed_val, "<empty>"))

# ----------------------------------------
# R Markdown Setup
# ----------------------------------------

setup_rmd <- function() {
  knitr::opts_chunk$set(warning = FALSE, message = FALSE)
  knitr::opts_knit$set(root.dir =getwd()) # Set working directory to project root
  if (nzchar(seed_val)) {
    set.seed(as.integer(seed_val))
  }
}

# Define working and results directories
workingDir <- getwd()
objDir <- file.path(workingDir, "obj")
debug_log("workingDir = %s", workingDir)
debug_log("objDir = %s", objDir)

# Load palettes, statistical functions, I/O utilities, and general utilities
source(file.path(objDir, "io_utils.R"))
source(file.path(objDir, "directories.R"))
source(file.path(objDir, "palettes.R"))
source(file.path(objDir, "plotting_fun.R"))

# ----------------------------------------
# Trait objects
# ----------------------------------------
## Trait data needs to have standardized colnames (species, family) and be comma-separated (csv format).

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
debug_log("trait_df rows = %d, cols = %d", nrow(trait_df), ncol(trait_df))
debug_log("trait_df columns: %s", paste(names(trait_df), collapse = ", "))

if (!"species" %in% names(trait_df)) {
  stop("Trait file must include a 'species' column.")
}
debug_log("trait_df species unique = %d", length(unique(trait_df$species)))

# ----------------------------------------
# Phylo objects
# ----------------------------------------

tree <- phytools::read.newick(tree_path) %>% ape::as.phylo() # Tree object
tree_species <- tree$tip.label # Tree species
debug_log("tree tips = %d, nodes = %d", length(tree$tip.label), tree$Nnode)

# ----------------------------------------
# Optional parameters (from YAML params)
# ----------------------------------------

clade_name <- if (exists("params")) params$clade_name else "clade"
taxon_of_interest <- if (exists("params")) params$taxon_of_interest else "family"
trait <- if (exists("params")) params$traitname else "trait"

debug_log("clade_name = %s", clade_name)
debug_log("taxon_of_interest = %s", taxon_of_interest)
debug_log("trait = %s", trait)

# Validate that required columns exist in trait_df
if (!taxon_of_interest %in% names(trait_df)) {
  stop(sprintf("Column '%s' (taxon_of_interest) not found in trait file. Available columns: %s", 
               taxon_of_interest, paste(names(trait_df), collapse=", ")))
}
if (!trait %in% names(trait_df)) {
  stop(sprintf("Column '%s' (trait) not found in trait file. Available columns: %s", 
               trait, paste(names(trait_df), collapse=", ")))
}

# Check if sample size (N) column is provided
source(file.path(objDir, "sample_size.R"))

# Load phylo.R to handle phylogenetic tree processing and tax_id mapping if needed
source(file.path(objDir, "phylo.R"))

# ----------------------------------------
# Secondary traits (optional)
# ----------------------------------------

secondary_trait <- if (exists("params")) params$secondary_trait else ""
debug_log("secondary_trait = %s", ifelse(nzchar(secondary_trait), secondary_trait, "<none>"))
has.secondary <- FALSE
if (nzchar(secondary_trait) && secondary_trait %in% names(trait_df)) {
  has.secondary <- TRUE
  debug_log("has.secondary = TRUE, missing = %d", sum(is.na(trait_df[[secondary_trait]])))
} else {
  message("No valid secondary trait provided; proceeding without it.")
  debug_log("has.secondary = FALSE")
}

branch_trait <- if (exists("params")) params$branch_trait else ""
debug_log("branch_trait = %s", ifelse(nzchar(branch_trait), branch_trait, "<none>"))
has.branch <- FALSE
if (nzchar(branch_trait) && branch_trait %in% names(trait_df)) {
  has.branch <- TRUE
  debug_log("has.branch = TRUE, missing = %d", sum(is.na(trait_df[[branch_trait]])))
} else {
  message("No valid branch trait provided; proceeding without it.")
  debug_log("has.branch = FALSE")
}

# ----------------------------------------
# Categorization / discretization method
# ----------------------------------------
# Used by stats.f() for global_label and by the CI-composition Rmd.
# Options: quartile, quintile, decile, median_sd, parameterized

discrete_method <- if (exists("params") && !is.null(params$discrete_method) && nzchar(params$discrete_method)) params$discrete_method else "decile"
top_quantile    <- if (exists("params") && !is.null(params$top_quantile)    && nzchar(params$top_quantile))    as.numeric(params$top_quantile)    else 0.90
bottom_quantile <- if (exists("params") && !is.null(params$bottom_quantile) && nzchar(params$bottom_quantile)) as.numeric(params$bottom_quantile) else 0.10
debug_log("discrete_method = %s, top_quantile = %.2f, bottom_quantile = %.2f",
          discrete_method, top_quantile, bottom_quantile)

# Maximum iterations for the contrast selection algorithm
contrast_max_iter <- if (exists("params") && !is.null(params$contrast_max_iter) && nzchar(params$contrast_max_iter)) as.integer(params$contrast_max_iter) else 4L
debug_log("contrast_max_iter = %d", contrast_max_iter)

# Lets try sourcing stats.R after all the parameters and data are loaded, since it relies on some of these variables being defined
source(file.path(objDir, "stats.R"))