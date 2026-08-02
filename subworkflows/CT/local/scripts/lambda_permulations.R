#!/usr/bin/env Rscript
# =============================================================================
# subworkflows/TRAIT_ANALYSIS/local/src/lambda_permulations.R
# Branch: perms_lambda
# =============================================================================
# Implements lambda-rescaled permulations for phenotype lists and benchmarks
# the computational cost of running the selection algorithm on each permulation.
#
# Logic for Lambda:
#   Taken from Correlations_cancer/src/permulate_response.R:
#   1. Fit ML Pagel's lambda using geiger::fitContinuous(tree, trait, model = "lambda")
#   2. Rescale tree via geiger::rescale(tree, "lambda", lambda_hat)
#   3. Draw permulated vectors via simpermvec(trait, phy_lambda)
#   4. Execute pair selection algorithm on each permulated trait vector
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Load selection algorithm modules
script_dir <- dirname(suppressWarnings(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)))
if (!is.character(script_dir) || length(script_dir) == 0 || !nzchar(script_dir)) {
  script_dir <- "subworkflows/TRAIT_ANALYSIS/local/src"
}
source(file.path(script_dir, "selection_algorithm.R"))

# ---- RERconverge functions (verbatim from PhyloPhere / permulate_response.R) ----
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm = geiger::ratematrix(treewithbranchlengths, namedvec)
  sims = geiger::sim.char(treewithbranchlengths, rm, nsim = 1)
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

# ---- Function to evaluate selection algorithm on a single trait vector ---------
run_selection_on_vector <- function(trait_vec, tree, has_n = FALSE, n_vec = NULL) {
  # Build long traits_df format expected by pair_sel.f
  df_trait <- data.frame(
    species = names(trait_vec),
    trait = "test_trait",
    value = as.numeric(trait_vec),
    stringsAsFactors = FALSE
  )
  if (has_n && !is.null(n_vec)) {
    df_trait$n_trait <- as.numeric(n_vec[names(trait_vec)])
  }
  
  dist_matrix <- calculate_patristic_distances(tree, df_trait)
  
  # Compute pairwise differences
  species_names <- df_trait$species
  n_sp <- length(species_names)
  overlap_list <- list()
  idx <- 1
  for (i in 1:(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      sp1 <- species_names[i]
      sp2 <- species_names[j]
      val1 <- trait_vec[sp1]
      val2 <- trait_vec[sp2]
      diff_val <- val1 - val2
      overlap_list[[idx]] <- data.frame(
        species1 = sp1, species2 = sp2, trait_diff = diff_val, stringsAsFactors = FALSE
      )
      idx <- idx + 1
      overlap_list[[idx]] <- data.frame(
        species1 = sp2, species2 = sp1, trait_diff = -diff_val, stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  overlap_df <- do.call(rbind, overlap_list)
  
  # Environment flag for pair_sel.f
  assign("has.n", has_n, envir = .GlobalEnv)
  
  # Run pair selection
  res <- try(pair_sel.f(dist_matrix, overlap_df, df_trait, "test_trait"), silent = TRUE)
  if (inherits(res, "try-error")) {
    return(list(n_pairs = 0, dunn_results = data.frame(), ok = FALSE))
  }
  
  n_pairs <- if (!is.null(res$dunn_results)) nrow(res$dunn_results) + 1 else 1
  return(list(n_pairs = n_pairs, dunn_results = res$dunn_results, ok = TRUE))
}

# ---- Benchmark function across a phenotype list -------------------------------
benchmark_phenotypes <- function(tree_file, trait_file, phenotype_list, n_perm = 50, seed = 42) {
  set.seed(seed)
  tree <- read.tree(tree_file)
  df_raw <- read.delim(trait_file, stringsAsFactors = FALSE)
  if (!("species" %in% colnames(df_raw))) stop("species column required in trait_file")
  
  cat("=========================================================================\n")
  cat(sprintf("BENCHMARKING LAMBDA-RESCALED PERMULATIONS & SELECTION ALGORITHM\n"))
  cat(sprintf("Phenotypes: %d | Permulations per phenotype: %d\n", length(phenotype_list), n_perm))
  cat("=========================================================================\n\n")
  
  summary_results <- list()
  
  for (pheno in phenotype_list) {
    cat(sprintf("---> Processing phenotype: %s <----\n", pheno))
    if (!(pheno %in% colnames(df_raw))) {
      cat(sprintf("  [SKIP] Column '%s' not found in trait file.\n\n", pheno))
      next
    }
    
    # Clean trait values
    sub_df <- df_raw[!is.na(df_raw[[pheno]]), c("species", pheno)]
    sub_df[[pheno]] <- as.numeric(sub_df[[pheno]])
    sub_df <- sub_df[is.finite(sub_df[[pheno]]), ]
    
    # Intersect with tree
    common_sp <- intersect(sub_df$species, tree$tip.label)
    if (length(common_sp) < 4) {
      cat(sprintf("  [SKIP] Fewer than 4 species shared with tree (%d found).\n\n", length(common_sp)))
      next
    }
    
    sub_df <- sub_df[sub_df$species %in% common_sp, ]
    sub_df <- sub_df[!duplicated(sub_df$species), ]
    phy <- keep.tip(tree, sub_df$species)
    phy <- multi2di(phy)
    phy$edge.length[phy$edge.length <= 0] <- 1e-8
    
    vec <- setNames(sub_df[[pheno]], sub_df$species)[phy$tip.label]
    
    # 1. Fit ML Lambda (Correlations Workflow logic)
    t0_fit <- proc.time()[["elapsed"]]
    fitl <- try(fitContinuous(phy, vec, model = "lambda"), silent = TRUE)
    t_fit <- proc.time()[["elapsed"]] - t0_fit
    
    if (inherits(fitl, "try-error")) {
      lambda_hat <- 1.0
      phy_lambda <- phy
      cat(sprintf("  [WARNING] ML Lambda fit failed. Falling back to BM (lambda = 1.0).\n"))
    } else {
      lambda_hat <- as.numeric(fitl$opt$lambda)
      phy_lambda <- rescale(phy, "lambda", lambda_hat)
      phy_lambda$edge.length[phy_lambda$edge.length <= 0] <- 1e-8
      cat(sprintf("  [ML FIT] Lambda-hat = %.4f (fit time: %.3f s)\n", lambda_hat, t_fit))
    }
    
    # 2. Observed selection baseline
    obs_sel <- run_selection_on_vector(vec, phy)
    n_obs_pairs <- obs_sel$n_pairs
    cat(sprintf("  [OBSERVED] Selection algorithm extracted %d pairs for real trait.\n", n_obs_pairs))
    
    # 3. Run Permulations & Benchmark
    t0_perm <- proc.time()[["elapsed"]]
    pair_counts <- numeric(n_perm)
    success_count <- 0
    
    for (b in seq_len(n_perm)) {
      pvec <- simpermvec(vec, phy_lambda)
      sel_res <- run_selection_on_vector(pvec, phy)
      pair_counts[b] <- sel_res$n_pairs
      if (sel_res$n_pairs >= n_obs_pairs) {
        success_count <- success_count + 1
      }
    }
    t_perm_total <- proc.time()[["elapsed"]] - t0_perm
    t_per_perm <- t_perm_total / n_perm
    
    pass_rate <- (success_count / n_perm) * 100
    avg_pairs <- mean(pair_counts)
    
    # Estimates for larger scale runs
    est_1k   <- (t_per_perm * 1000) / 60
    est_10k  <- (t_per_perm * 10000) / 60
    est_100k <- (t_per_perm * 100000) / 3600
    
    cat(sprintf("  [BENCHMARK] Total perm time (%d perms): %.2f s (%.4f s / perm)\n", n_perm, t_perm_total, t_per_perm))
    cat(sprintf("  [YIELD] Permulations with pairs >= observed (%d): %d / %d (%.2f%%)\n", n_obs_pairs, success_count, n_perm, pass_rate))
    cat(sprintf("  [METRICS] Avg pairs extracted per perm: %.2f (range %d .. %d)\n", avg_pairs, min(pair_counts), max(pair_counts)))
    cat(sprintf("  [COST PROJECTION] Estimated time for:\n"))
    cat(sprintf("    - 1,000 perms  : %.2f mins\n", est_1k))
    cat(sprintf("    - 10,000 perms : %.2f mins\n", est_10k))
    cat(sprintf("    - 100,000 perms: %.2f hours\n\n", est_100k))
    
    summary_results[[pheno]] <- data.frame(
      phenotype = pheno,
      n_species = length(vec),
      lambda_hat = round(lambda_hat, 4),
      n_obs_pairs = n_obs_pairs,
      sec_per_perm = round(t_per_perm, 4),
      pass_rate_pct = round(pass_rate, 2),
      avg_perm_pairs = round(avg_pairs, 2),
      est_mins_1k = round(est_1k, 2),
      est_mins_10k = round(est_10k, 2),
      stringsAsFactors = FALSE
    )
  }
  
  df_summary <- do.call(rbind, summary_results)
  cat("=========================================================================\n")
  cat("FINAL PHENOTYPE BENCHMARK SUMMARY\n")
  cat("=========================================================================\n")
  print(df_summary)
  return(df_summary)
}

# Run directly if executed as standalone script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  tree_f  <- if (length(args) >= 1) args[1] else "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_true2/malignant_prevalence/data_exploration/2.CT/3.Tree/pruned_tree_file.nwk"
  trait_f <- if (length(args) >= 2) args[2] else "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_true2/malignant_prevalence/data_exploration/0.Data-pruning/pruned_trait_file.tsv"
  n_p     <- if (length(args) >= 3) as.integer(args[3]) else 50
  
  pheno_list <- c(
    "malignant_prevalence",
    "neoplasia_prevalence",
    "benign_prevalence",
    "body_mass",
    "LQ",
    "mature_age"
  )
  
  benchmark_phenotypes(tree_f, trait_f, pheno_list, n_perm = n_p)
}
