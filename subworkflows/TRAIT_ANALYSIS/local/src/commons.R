# commons.R — Shared setup of the trait-analysis reports: parameters, trait table, tree and helpers.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: 0.Data_pruning.Rmd, 1.Dataset_exploration.Rmd, 2.Phenotype_exploration.Rmd,
#             3.CI-composition.Rmd, 4.Independent_contrasts.Rmd (setup chunk, then setup_rmd())
#
# Reads the rmarkdown `params` of the calling report (trait_file, tree_file, output_dir,
# seed, clade_name, taxon_of_interest, sp_colname, traitname, secondary_trait,
# branch_trait, trait_type, pss_top_pct, perm_strategy, max_contrasts), loads the trait
# table (trait_df, with a `species` column) and the tree (tree), resolves the optional
# secondary and branch traits (has.secondary, has.branch), and sources the other files
# of src/. It stops when a required parameter or column is missing. The report must
# be rendered with the working directory that holds src/.
# =============================================================================


# ── Debug logging ──────────────────────────────────────────────────────────────

# debug_log() prints "[DEBUG] ..." to stderr and appends the message to
# phylo_debug_log, which the reports can print into the HTML.

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

# Evaluates `expr`, logging its start time and elapsed seconds, and returns its value.
debug_stage <- function(label, expr) {
  debug_log("[STAGE START] %s @ %s", label, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  debug_log("[STAGE END] %s elapsed=%.3fs", label, proc.time()[["elapsed"]] - t0)
  value
}

# ── Parameters ────────────────────────────────────────────────────────────────

# get_arg() reads positional command-line arguments and warns on every call. No report
# uses it; parameters come from the YAML `params` (params$trait_file, params$seed).
get_arg <- function(args, idx, default = NULL) {
  warning("get_arg() is deprecated. Use params$ from YAML header instead.")
  if (length(args) >= idx && nzchar(args[idx])) {
    return(args[idx])
  }
  default
}

# The calling report passes `params` through rmarkdown::render(); the required ones stop the run when absent.
trait_path <- if (exists("params") && !is.null(params$trait_file)) params$trait_file else stop("trait_file parameter required")
tree_path <- if (exists("params") && !is.null(params$tree_file)) params$tree_file else stop("tree_file parameter required")
resultsDir <- if (exists("params")) params$output_dir else getwd()
seed_val <- if (exists("params")) params$seed else ""

debug_log("trait_path = %s", trait_path)
debug_log("tree_path = %s", tree_path)
debug_log("resultsDir = %s", resultsDir)
debug_log("seed_val = %s", ifelse(nzchar(seed_val), seed_val, "<empty>"))

# ── R Markdown setup ──────────────────────────────────────────────────────────

# Sets the chunk defaults (silent chunks, no echo), fixes the knitr root directory to the
# working directory and seeds the RNG when a seed was given.
setup_rmd <- function() {
  knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo = FALSE)
  knitr::opts_knit$set(root.dir =getwd()) # the working directory that holds src/
  if (nzchar(seed_val)) {
    set.seed(as.integer(seed_val))
  }
}

# Working directory and the src/ directory that holds the other helper files.
workingDir <- getwd()
objDir <- file.path(workingDir, "src")
debug_log("workingDir = %s", workingDir)
debug_log("objDir = %s", objDir)

# Helper files: I/O utilities, results directories, palettes and plotting functions.
source(file.path(objDir, "io_utils.R"))
source(file.path(objDir, "directories.R"))
source(file.path(objDir, "palettes.R"))
source(file.path(objDir, "plotting_fun.R"))

# ── Trait table ───────────────────────────────────────────────────────────────

# The trait file is read as tab-separated when it ends in .tsv and as comma-separated otherwise.
# It needs a species column (named by sp_colname) and the taxon_of_interest and trait columns.

print(paste0("Loading trait data from: ", trait_path))

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

sp_colname <- if (exists("params") && !is.null(params$sp_colname) && nzchar(params$sp_colname)) params$sp_colname else "species"
debug_log("sp_colname = %s", sp_colname)

if (!sp_colname %in% names(trait_df)) {
  stop(sprintf("Trait file must include a '%s' column (sp_colname). Available columns: %s",
               sp_colname, paste(names(trait_df), collapse = ", ")))
}
if (sp_colname != "species") {
  trait_df$species <- trait_df[[sp_colname]]
}
debug_log("trait_df species unique = %d", length(unique(trait_df$species)))

# ── Species tree ──────────────────────────────────────────────────────────────
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

# ── Ultrametric check ─────────────────────────────────────────────────────────

# Contrast independence (the modified Dunn index in selection_algorithm.R) and the
# OU/BM Phylogenetic Shift Score (pss_core.R) both assume a time tree. On a
# phylogram, a fast-evolving lineage has a long terminal branch that reflects rate,
# not time: it inflates the diameter of its contrast pair and its patristic distance
# to its sister, so a valid independent contrast can get a Dunn index below 1 and be
# dropped. PhyloPhere expects a dated (ultrametric) species tree and only warns here:
# rooting and penalized-likelihood dating of an arbitrary phylogram is not robust
# (midpoint rooting fails under the same rate variation that motivates it).
if (!ape::is.ultrametric(tree, tol = 1e-6)) {
  warning("commons.R: input tree is NOT ultrametric (phylogram). Contrast ",
          "independence (Dunn) and the OU/BM PSS assume a time tree; ",
          "rate-variation branch-length artifacts may distort pair selection. ",
          "Supply a dated species tree.")
}

tree_species <- tree$tip.label
debug_log("tree tips = %d, nodes = %d, ultrametric = %s",
          length(tree$tip.label), tree$Nnode, ape::is.ultrametric(tree, tol = 1e-6))

# ── Clade, taxon and trait ────────────────────────────────────────────────────

# Defaults apply only when the report has no `params`; the taxon and trait columns must exist.

clade_name <- if (exists("params")) params$clade_name else "clade"
taxon_of_interest <- if (exists("params")) params$taxon_of_interest else "family"
trait <- if (exists("params")) params$traitname else "trait"

debug_log("clade_name = %s", clade_name)
debug_log("taxon_of_interest = %s", taxon_of_interest)
debug_log("trait = %s", trait)

if (!taxon_of_interest %in% names(trait_df)) {
  stop(sprintf("Column '%s' (taxon_of_interest) not found in trait file. Available columns: %s", 
               taxon_of_interest, paste(names(trait_df), collapse=", ")))
}
if (!trait %in% names(trait_df)) {
  stop(sprintf("Column '%s' (trait) not found in trait file. Available columns: %s", 
               trait, paste(names(trait_df), collapse=", ")))
}

# Detects the optional count columns (n_trait, c_trait) and the sample-size column.
source(file.path(objDir, "sample_size.R"))

# Match the trait table to the tree tips.
source(file.path(objDir, "phylo.R"))

# ── Secondary and branch traits ───────────────────────────────────────────────

# Resolves a requested trait name against the trait_df columns, tolerating case
# differences between the parameter and the file header (params `branch_trait = "LQ"`
# against a column named "lq"). Returns the column name, or NA when there is no match.
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

# ── Trait type and PSS parameters ─────────────────────────────────────────────

# trait_type: "auto" or the forced type; pss_top_pct: fraction of the hi>lo pairs kept by the
# continuous-trait gate; perm_strategy: evolutionary model strategy (best_model, bm or ou).
trait_type <- if (exists("params") && !is.null(params$trait_type) && nzchar(params$trait_type)) tolower(params$trait_type) else "auto"
pss_top_pct <- if (exists("params") && !is.null(params$pss_top_pct) && nzchar(as.character(params$pss_top_pct))) as.numeric(params$pss_top_pct) else 0.01
perm_strategy <- if (exists("params") && !is.null(params$perm_strategy) && nzchar(params$perm_strategy)) params$perm_strategy else "best_model"

debug_log("trait_type = %s, pss_top_pct = %.4f, perm_strategy = %s", trait_type, pss_top_pct, perm_strategy)

# Maximum number of contrasts for the selection algorithm (0 or unset gives Inf: as many as the data support).
max_contrasts <- if (exists("params") && !is.null(params$max_contrasts) && nzchar(as.character(params$max_contrasts)) && as.integer(params$max_contrasts) > 0L) {
  as.integer(params$max_contrasts)
} else {
  Inf
}
debug_log("max_contrasts = %s", ifelse(is.finite(max_contrasts), as.character(max_contrasts), "<dynamic>"))

# stats.R relies on the variables defined above, so it is sourced last.
source(file.path(objDir, "stats.R"))
