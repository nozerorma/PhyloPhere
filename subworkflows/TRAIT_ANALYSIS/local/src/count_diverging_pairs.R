#!/usr/bin/env Rscript
# =============================================================================
# subworkflows/TRAIT_ANALYSIS/local/src/count_diverging_pairs.R
# Branch: perms_lambda
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

run_pair_count <- function(trait_vec, tree) {
  D <- cophenetic(tree)
  sp_names <- names(trait_vec)
  n_sp <- length(sp_names)
  pairs_list <- list()
  idx <- 1
  for (i in 1:(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      s1 <- sp_names[i]; s2 <- sp_names[j]
      abs_diff <- abs(trait_vec[s1] - trait_vec[s2])
      if (abs_diff > 0) {
        pairs_list[[idx]] <- data.frame(sp1 = s1, sp2 = s2, dist = as.numeric(D[s1, s2]), abs_diff = as.numeric(abs_diff), stringsAsFactors = FALSE)
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

n_perm <- 100
set.seed(42)
counts <- numeric(n_perm)

for (b in seq_len(n_perm)) {
  pvec <- simpermvec(vec, phy_lam)
  counts[b] <- run_pair_count(pvec, phy)
}

cat("=========================================================================\n")
cat(sprintf("PAIR COUNTS DISTRIBUTION FOR MALIGNANT_PREVALENCE (%d PERMUTATIONS)\n", n_perm))
cat("=========================================================================\n")
tbl <- table(counts)
print(tbl)
cat("\n")
cat(sprintf("Total Permulations            : %d\n", n_perm))
cat(sprintf("Permulations with EXACTLY 3   : %d (%.1f%%)\n", sum(counts == 3), 100*mean(counts == 3)))
cat(sprintf("Permulations with > 3 pairs   : %d (%.1f%%)\n", sum(counts > 3), 100*mean(counts > 3)))
cat(sprintf("Permulations with >= 4 pairs  : %d (%.1f%%)\n", sum(counts >= 4), 100*mean(counts >= 4)))
cat(sprintf("Permulations with >= 5 pairs  : %d (%.1f%%)\n", sum(counts >= 5), 100*mean(counts >= 5)))
