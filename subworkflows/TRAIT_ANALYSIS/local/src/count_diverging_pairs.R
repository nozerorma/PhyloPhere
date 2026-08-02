#!/usr/bin/env Rscript
# =============================================================================
# subworkflows/TRAIT_ANALYSIS/local/src/count_diverging_pairs.R
# Branch: perms_lambda
# =============================================================================
# Counts diverging pairs under Quintile Discretization (Top Q80 vs Bottom Q20)
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(dplyr)
  library(tidyr)
})

tree_f  <- "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_true2/malignant_prevalence/data_exploration/2.CT/3.Tree/pruned_tree_file.nwk"
trait_f <- "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final_true2/malignant_prevalence/data_exploration/0.Data-pruning/pruned_trait_file.tsv"

tree <- read.tree(tree_f)
df_raw <- read.delim(trait_f, stringsAsFactors = FALSE)

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

mod_dunn <- function(distance, clusters, selected_cluster) {
  c1 <- which(clusters == selected_cluster)
  intraClust <- if (length(c1) > 1) max(distance[c1, c1]) else 0
  if (intraClust == 0) return(Inf)
  nc <- max(clusters)
  interClust <- rep(NA, nc)
  for (j in seq_len(nc)) {
    if (j == selected_cluster) next
    c2 <- which(clusters == j)
    interClust[j] <- min(distance[c1, c2])
  }
  min(interClust, na.rm = TRUE) / intraClust
}

# Count pairs enforcing Quintile discretization (Top >= Q80 vs Bottom <= Q20)
run_pair_count_quintiles <- function(trait_vec, tree) {
  D <- cophenetic(tree)
  sp_names <- names(trait_vec)
  
  # Quintile thresholds
  q20 <- unname(quantile(trait_vec, 0.20, na.rm = TRUE))
  q80 <- unname(quantile(trait_vec, 0.80, na.rm = TRUE))
  
  top_sp <- sp_names[trait_vec >= q80]
  bot_sp <- sp_names[trait_vec <= q20]
  
  pairs_list <- list()
  idx <- 1
  for (t_s in top_sp) {
    for (b_s in bot_sp) {
      if (t_s == b_s) next
      abs_diff <- abs(trait_vec[t_s] - trait_vec[b_s])
      if (abs_diff > 0) {
        pairs_list[[idx]] <- data.frame(sp1 = t_s, sp2 = b_s, dist = as.numeric(D[t_s, b_s]), abs_diff = as.numeric(abs_diff), stringsAsFactors = FALSE)
        idx <- idx + 1
      }
    }
  }
  pairs_df <- do.call(rbind, pairs_list)
  if (is.null(pairs_df) || nrow(pairs_df) == 0) return(0)
  pairs_df <- pairs_df[order(pairs_df$dist, -pairs_df$abs_diff), ]
  
  selected_sp <- data.frame(sp = c(pairs_df$sp1[1], pairs_df$sp2[1]), cluster = c(1, 1), stringsAsFactors = FALSE)
  n_clusters <- 1
  
  while (n_clusters < 10) {
    cand_found <- FALSE
    for (r in seq_len(nrow(pairs_df))) {
      q1 <- pairs_df$sp1[r]; q2 <- pairs_df$sp2[r]
      if (q1 %in% selected_sp$sp || q2 %in% selected_sp$sp) next
      
      test_sp <- rbind(selected_sp, data.frame(sp = c(q1, q2), cluster = c(n_clusters + 1, n_clusters + 1)))
      dsub <- D[test_sp$sp, test_sp$sp]
      cl_vec <- setNames(test_sp$cluster, test_sp$sp)
      
      md <- min(sapply(1:(n_clusters + 1), function(k) mod_dunn(dsub, cl_vec, k)))
      if (md >= 1.0) {
        selected_sp <- test_sp
        n_clusters <- n_clusters + 1
        cand_found <- TRUE
        break
      }
    }
    if (!cand_found) break
  }
  return(n_clusters)
}

vals <- suppressWarnings(as.numeric(df_raw$malignant_prevalence))
good <- !is.na(vals) & is.finite(vals)
sub_sp <- df_raw$species[good]
sub_val <- vals[good]
common_sp <- intersect(sub_sp, tree$tip.label)

phy <- keep.tip(tree, common_sp)
phy <- multi2di(phy); phy$edge.length[phy$edge.length <= 0] <- 1e-8
vec <- setNames(sub_val, sub_sp)[phy$tip.label]

fitl <- fitContinuous(phy, vec, model = "lambda")
lambda_hat <- as.numeric(fitl$opt$lambda)
phy_lam <- rescale(phy, "lambda", lambda_hat)
phy_lam$edge.length[phy_lam$edge.length <= 0] <- 1e-8

n_perm <- 200
set.seed(42)
counts <- numeric(n_perm)

for (b in seq_len(n_perm)) {
  pvec <- simpermvec(vec, phy_lam)
  counts[b] <- run_pair_count_quintiles(pvec, phy)
}

cat("=========================================================================\n")
cat(sprintf("QUINTILE-RESTRICTED PAIR COUNTS DISTRIBUTION (%d PERMUTATIONS)\n", n_perm))
cat(sprintf("Phenotype: malignant_prevalence | ML lambda = %.4f\n", lambda_hat))
cat("=========================================================================\n")
tbl <- table(counts)
print(tbl)
cat("\n")
cat(sprintf("Total Permulations                      : %d\n", n_perm))
cat(sprintf("Permulations with EXACTLY 3 pairs       : %d (%.1f%%)\n", sum(counts == 3), 100*mean(counts == 3)))
cat(sprintf("Permulations with MORE than 3 pairs (>3): %d (%.1f%%)\n", sum(counts > 3), 100*mean(counts > 3)))
cat(sprintf("Permulations with LESS than 3 pairs (<3): %d (%.1f%%)\n", sum(counts < 3), 100*mean(counts < 3)))
cat(sprintf("Average Quintile Pairs Extracted        : %.2f\n", mean(counts)))
