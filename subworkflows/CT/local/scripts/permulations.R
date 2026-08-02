#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/permulations.R
# Branch: perms_lambda
# =============================================================================
# Generates null trait permulations supporting BM, FGBG, and Pagel's Lambda modes.
# Implements 2-Tier pool harvesting (Tier 1: N_obs, Tier 2: N_obs - 1) with overall Dunn index validation.
# Auto-detects observed independent pair count (N_pairs_obs) directly from config.file (V3 cluster IDs).
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
  library(ape)
  library(geiger)
  library(dplyr)
})

# Functions used by CT Resample
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm = ratematrix(treewithbranchlengths, namedvec)
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  simulatedvec
}

simpermvec <- function(namedvec, treewithbranchlengths) {
  vec = simulatevec(namedvec, treewithbranchlengths)
  simsorted = sort(vec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}

# Parse positional CLI arguments
args = commandArgs(trailingOnly=TRUE)

tree               <- args[1]
config.file        <- args[2]
number.of.cycles   <- as.integer(args[3])
selection.strategy <- tolower(as.character(args[4])) # "bm", "fgbg", or "lambda"
phenotypes         <- args[5]
outdir             <- args[6]
chunk.size         <- ifelse(length(args) >= 7, as.integer(args[7]), 500)
include.b0         <- TRUE
if (length(args) >= 8) {
  include.b0 <- tolower(as.character(args[8])) %in% c("1", "true", "t", "yes", "y")
}
target_pairs    <- ifelse(length(args) >= 9, as.integer(args[9]), 0L)
discrete_method <- ifelse(length(args) >= 10, as.character(args[10]), "quintile")
max_tries       <- ifelse(length(args) >= 11, as.integer(args[11]), 1000000L)
pheno_col_name  <- ifelse(length(args) >= 12, as.character(args[12]), "")

# Create output directory
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  write(paste("[INFO]", Sys.time(), "Created output directory:", outdir), stdout())
}

tree.o <- read.tree(tree)

# Read config file (species | binary label | cluster)
cfg <- read.table(config.file, sep ="\t", header = FALSE, stringsAsFactors = FALSE)
foreground.df <- subset(cfg, cfg$V2 == "1")
background.df <- subset(cfg, cfg$V2 == "0")

foreground.species <- foreground.df$V1
background.species <- background.df$V1

# Auto-detect observed pair count (N_pairs_obs) from config V3 cluster IDs if target_pairs <= 0
if (target_pairs <= 0) {
  if (ncol(cfg) >= 3) {
    auto_obs_pairs <- max(suppressWarnings(as.integer(cfg$V3)), na.rm = TRUE)
    if (is.finite(auto_obs_pairs) && auto_obs_pairs > 0) {
      target_pairs <- auto_obs_pairs
      write(paste("[INFO]", Sys.time(), sprintf("Auto-detected observed independent pair count from config (V3 cluster IDs): N_pairs_obs = %d", target_pairs)), stdout())
    }
  }
}

# Read phenotype file (supports 2-col TSV or multi-col TSV with header)
pheno_raw <- read.delim(
  phenotypes,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN", "nan", "NULL", "null")
)

# Determine species and trait columns
if ("species" %in% names(pheno_raw)) {
  sp_col <- "species"
} else if ("V1" %in% names(pheno_raw)) {
  sp_col <- "V1"
} else {
  sp_col <- names(pheno_raw)[1]
}

if (nzchar(pheno_col_name) && pheno_col_name %in% names(pheno_raw)) {
  val_col <- pheno_col_name
} else if ("V2" %in% names(pheno_raw)) {
  val_col <- "V2"
} else {
  num_cols <- names(pheno_raw)[sapply(pheno_raw, function(x) canCoerce(x, "numeric") || is.numeric(x))]
  num_cols <- setdiff(num_cols, sp_col)
  val_col <- if (length(num_cols) > 0) num_cols[1] else names(pheno_raw)[2]
}

phenotype.df <- data.frame(
  V1 = as.character(pheno_raw[[sp_col]]),
  V2 = suppressWarnings(as.numeric(pheno_raw[[val_col]])),
  stringsAsFactors = FALSE
)

phenotype.df <- subset(phenotype.df, !is.na(V1) & !is.na(V2) & is.finite(V2))

if (nrow(phenotype.df) == 0) {
  stop("No valid phenotype rows remain after NA filtering")
}

available.species <- phenotype.df$V1
foreground.species <- intersect(foreground.species, available.species)
background.species <- intersect(background.species, available.species)

foreground.size <- length(foreground.species)
background.size <- length(background.species)

if (foreground.size == 0 || background.size == 0) {
  stop("Foreground or background group is empty after missing-trait filtering")
}

foreground.values <- subset(phenotype.df, phenotype.df$V1 %in% foreground.species)$V2
background.values <- subset(phenotype.df, phenotype.df$V1 %in% background.species)$V2

starting.values <- phenotype.df$V2
all.species     <- phenotype.df$V1
names(starting.values) <- all.species

pruned.tree.o <- drop.tip(tree.o, setdiff(tree.o$tip.label, all.species))
starting.values <- starting.values[names(starting.values) %in% pruned.tree.o$tip.label]

# Ensure tree is strictly bifurcating with positive branch lengths
pruned.tree.o <- multi2di(pruned.tree.o)
pruned.tree.o$edge.length[pruned.tree.o$edge.length <= 0] <- 1e-8

# Tree Rescaling for Pagel's Lambda Mode
simulation_tree <- pruned.tree.o
if (selection.strategy == "lambda") {
  write(paste("[INFO]", Sys.time(), "Fitting ML Pagel's lambda on empirical trait vector..."), stdout())
  fitl <- fitContinuous(pruned.tree.o, starting.values, model = "lambda")
  lambda_hat <- as.numeric(fitl$opt$lambda)
  write(paste("[INFO]", Sys.time(), sprintf("Estimated ML Pagel's lambda: %.4f", lambda_hat)), stdout())
  
  simulation_tree <- rescale(pruned.tree.o, "lambda", lambda_hat)
  simulation_tree$edge.length[simulation_tree$edge.length <= 0] <- 1e-8
}

# Source lean contrast selector if pair filtering is requested or lambda strategy is active
use_lean_filter <- (target_pairs > 0)
if (use_lean_filter) {
  possible_paths <- c(
    file.path(dirname(normalizePath(sys.frames()[[1]]$ofile %||% ".")), "lean_contrast_selector.R"),
    "subworkflows/CT/local/scripts/lean_contrast_selector.R",
    "/homes/users/mramon/scratch/PhyloPhere/subworkflows/CT/local/scripts/lean_contrast_selector.R"
  )
  lean_script <- possible_paths[file.exists(possible_paths)][1]
  
  if (!is.na(lean_script) && file.exists(lean_script)) {
    source(lean_script)
    write(paste("[INFO]", Sys.time(), sprintf("Lean contrast filtering ENABLED from %s (target_pairs = %d, discrete_method = %s)", lean_script, target_pairs, discrete_method)), stdout())
  } else {
    warning("lean_contrast_selector.R not found; proceeding without pair count filtering")
    use_lean_filter <- FALSE
  }
}

# Write b_0 if requested
if (include.b0) {
  write(paste("[INFO]", Sys.time(), "Writing b_0 (original trait configuration) to resample_000.tab"), stdout())
  b0.df <- data.frame(matrix(ncol = 3, nrow = 1))
  b0.outline <- c("b_0", paste(foreground.species, collapse=","), paste(background.species, collapse=","))
  b0.df[1,] <- b0.outline
  b0.filepath <- file.path(outdir, "resample_000.tab")
  write.table(b0.df, sep="\t", col.names=FALSE, row.names=FALSE, file=b0.filepath, quote=FALSE)
}

# HARVESTING LOOP WITH TIER 1 & TIER 2 POOLING
start.time <- Sys.time()
write(paste("[START]", start.time, sprintf("Beginning permulation generation (strategy: %s, max_tries: %d)...", selection.strategy, max_tries)), stdout())

tier1_list <- list() # N_pairs_perm == target_pairs & mod_dunn >= 1.0
tier2_list <- list() # N_pairs_perm == target_pairs - 1 & mod_dunn >= 1.0

total_draws <- 0L

while (length(tier1_list) < number.of.cycles && total_draws < max_tries) {
  total_draws <- total_draws + 1L
  
  if (selection.strategy == "fgbg") {
    permulated_phenotype <- sample(starting.values)
  } else {
    permulated_phenotype <- simpermvec(starting.values, simulation_tree)
  }
  
  if (use_lean_filter) {
    eval_res <- evaluate_lean_contrast_selection(
      trait_vec = permulated_phenotype,
      tree = pruned.tree.o,
      target_pairs = target_pairs,
      discrete_method = discrete_method
    )
    
    if (eval_res$tier == 1L) {
      tier1_list[[length(tier1_list) + 1L]] <- permulated_phenotype
    } else if (eval_res$tier == 2L && length(tier2_list) < number.of.cycles) {
      tier2_list[[length(tier2_list) + 1L]] <- permulated_phenotype
    }
  } else {
    tier1_list[[length(tier1_list) + 1L]] <- permulated_phenotype
  }
  
  if (total_draws %% 500 == 0 || length(tier1_list) == number.of.cycles) {
    elapsed <- as.numeric(difftime(Sys.time(), start.time, units="secs"))
    write(sprintf("[%s] Draws: %d | Tier 1: %d/%d | Tier 2: %d | Elapsed: %.1f min",
                  format(Sys.time(), "%H:%M:%S"),
                  total_draws,
                  length(tier1_list),
                  number.of.cycles,
                  length(tier2_list),
                  elapsed / 60), stdout())
  }
}

# Assemble final accepted pool of permulations (Priority: Tier 1 -> Tier 2)
n_tier1 <- length(tier1_list)
needed  <- number.of.cycles - n_tier1

if (needed <= 0) {
  final_pool <- tier1_list[1:number.of.cycles]
  write(paste("[INFO]", Sys.time(), sprintf("Target pool achieved using 100%% Tier 1 permulations (%d/%d)", number.of.cycles, number.of.cycles)), stdout())
} else {
  n_tier2_use <- min(needed, length(tier2_list))
  final_pool <- c(tier1_list, tier2_list[seq_len(n_tier2_use)])
  write(paste("[WARN]", Sys.time(), sprintf("Harvested %d Tier 1 permulations + %d Tier 2 permulations (Total: %d/%d)",
                                            n_tier1, n_tier2_use, length(final_pool), number.of.cycles)), stdout())
}

if (length(final_pool) == 0) {
  stop("Failed to harvest any valid permulations within max_tries limit")
}

# Write resample chunks to disk
simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))
chunk.counter <- 1
file.counter  <- 1

for (b in seq_along(final_pool)) {
  pvec <- final_pool[[b]]
  cycle.tag <- paste("b", as.character(b), sep = "_")
  x <- enframe(pvec)
  
  potential.fg.df <- subset(x, value %in% foreground.values)
  potential.bg.df <- subset(x, value %in% background.values)
  
  fg.species <- sample(potential.fg.df$name, foreground.size)
  bg.species <- sample(potential.bg.df$name, background.size)
  
  outline <- c(cycle.tag, paste(fg.species, collapse=","), paste(bg.species, collapse=","))
  simulated.traits.df[nrow(simulated.traits.df) + 1, ] <- outline
  
  chunk.counter <- chunk.counter + 1
  
  if (chunk.counter > chunk.size || b == length(final_pool)) {
    filename <- sprintf("resample_%03d.tab", file.counter)
    filepath <- file.path(outdir, filename)
    
    write.table(simulated.traits.df, sep="\t", col.names=FALSE, row.names=FALSE, file=filepath, quote=FALSE)
    
    simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))
    chunk.counter <- 1
    file.counter  <- file.counter + 1
  }
}

end.time <- Sys.time()
total.elapsed <- as.numeric(difftime(end.time, start.time, units="mins"))
write(paste("[COMPLETE]", end.time, "| Generated", length(final_pool), "permulation records from", total_draws, "draws in", round(total.elapsed, 2), "minutes"), stdout())
