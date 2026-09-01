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
    cat("[DEBUG] ", msg, "\n", sep = "", file = stderr())
    flush.console()
  }
}

debug_stage <- function(label, expr) {
  debug_log("[STAGE START] %s @ %s", label, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  debug_log("[STAGE END] %s elapsed=%.3fs", label, proc.time()[["elapsed"]] - t0)
  value
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
  knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo = FALSE)
  knitr::opts_knit$set(root.dir =getwd()) # Set working directory to project root
  if (nzchar(seed_val)) {
    set.seed(as.integer(seed_val))
  }
}

# Define working and results directories
workingDir <- getwd()
objDir <- file.path(workingDir, "src")
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

debug_log("tree_path exists = %s", file.exists(tree_path))
tree_preview <- tryCatch(
  readLines(tree_path, n = 2, warn = FALSE),
  error = function(e) paste0("<readLines error: ", conditionMessage(e), ">")
)
debug_log("tree preview: %s", paste(tree_preview, collapse = " | "))

tree <- debug_stage(
  "read tree",
  ape::read.tree(file = tree_path)
)
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

# Resolve a requested trait name against the trait_df columns, tolerating
# case differences between GUI/params defaults and curated file headers
# (e.g. params `branch_trait = "LQ"` vs a column literally named "lq").
resolve_trait_column <- function(trait_name, df_names) {
  if (!nzchar(trait_name)) return(NA_character_)
  if (trait_name %in% df_names) return(trait_name)
  match_idx <- match(tolower(trait_name), tolower(df_names))
  if (!is.na(match_idx)) return(df_names[match_idx])
  NA_character_
}

secondary_trait_requested <- if (exists("params")) params$secondary_trait else ""
secondary_trait_resolved <- resolve_trait_column(secondary_trait_requested, names(trait_df))
debug_log("secondary_trait = %s", ifelse(nzchar(secondary_trait_requested), secondary_trait_requested, "<none>"))
has.secondary <- FALSE
if (!is.na(secondary_trait_resolved)) {
  secondary_trait <- secondary_trait_resolved
  has.secondary <- TRUE
  debug_log("has.secondary = TRUE, resolved column = %s, missing = %d", secondary_trait, sum(is.na(trait_df[[secondary_trait]])))
} else {
  secondary_trait <- secondary_trait_requested
  message("No valid secondary trait provided; proceeding without it.")
  debug_log("has.secondary = FALSE")
}

branch_trait_requested <- if (exists("params")) params$branch_trait else ""
branch_trait_resolved <- resolve_trait_column(branch_trait_requested, names(trait_df))
debug_log("branch_trait = %s", ifelse(nzchar(branch_trait_requested), branch_trait_requested, "<none>"))
has.branch <- FALSE
if (!is.na(branch_trait_resolved)) {
  branch_trait <- branch_trait_resolved
  has.branch <- TRUE
  debug_log("has.branch = TRUE, resolved column = %s, missing = %d", branch_trait, sum(is.na(trait_df[[branch_trait]])))
} else {
  branch_trait <- branch_trait_requested
  message("No valid branch trait provided; proceeding without it.")
  debug_log("has.branch = FALSE")
}

# ----------------------------------------
# Trait type and PSS parameters
# ----------------------------------------
trait_type <- if (exists("params") && !is.null(params$trait_type) && nzchar(params$trait_type)) tolower(params$trait_type) else "auto"
pss_top_pct <- if (exists("params") && !is.null(params$pss_top_pct) && nzchar(as.character(params$pss_top_pct))) as.numeric(params$pss_top_pct) else 0.01
perm_strategy <- if (exists("params") && !is.null(params$perm_strategy) && nzchar(params$perm_strategy)) params$perm_strategy else "best_model"

debug_log("trait_type = %s, pss_top_pct = %.4f, perm_strategy = %s", trait_type, pss_top_pct, perm_strategy)

# Maximum contrasts for the contrast selection algorithm (0 or Inf = dynamic discovery)
max_contrasts <- if (exists("params") && !is.null(params$max_contrasts) && nzchar(as.character(params$max_contrasts)) && as.integer(params$max_contrasts) > 0L) {
  as.integer(params$max_contrasts)
} else if (exists("params") && !is.null(params$contrast_max_iter) && nzchar(as.character(params$contrast_max_iter)) && as.integer(params$contrast_max_iter) > 0L) {
  as.integer(params$contrast_max_iter)
} else {
  Inf
}
debug_log("max_contrasts = %s", ifelse(is.finite(max_contrasts), as.character(max_contrasts), "<dynamic>"))

# Lets try sourcing stats.R after all the parameters and data are loaded, since it relies on some of these variables being defined
source(file.path(objDir, "stats.R"))
