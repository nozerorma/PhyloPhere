#!/usr/bin/env Rscript
# =============================================================================
# subworkflows/CT/local/scripts/lean_contrast_selector.R
# Branch: perms_lambda
# =============================================================================
# Lean contrast selection & Dunn validation module for CAASTOOLS (CT).
# Evaluates whether a permulated trait vector yields N_pairs_obs independent
# contrasts (Tier 1) or N_pairs_obs - 1 contrasts (Tier 2) with overall mod_dunn >= 1.0.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(dplyr)
  library(tidyr)
})

# Modified Dunn index for phylogenetic cluster validation
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

# Calculate overall minimum Dunn index across all clusters
calculate_overall_dunn <- function(D_sub, cluster_vec) {
  nc <- max(cluster_vec)
  if (nc <= 1) return(Inf)
  dunn_vals <- sapply(seq_len(nc), function(k) mod_dunn(D_sub, cluster_vec, k))
  min(dunn_vals, na.rm = TRUE)
}

#' Run Lean Contrast Selection and Tiered Dunn Validation
#' 
#' @param trait_vec Named numeric vector of species trait values.
#' @param tree phylo object.
#' @param target_pairs Integer: target number of independent pairs to extract (N_pairs_obs).
#' @param discrete_method Discretization method ("quintile", "quartile", "decile", "median_sd", "parameterized").
#' @param top_quantile Numeric quantile threshold for top extreme (default: 0.80).
#' @param bottom_quantile Numeric quantile threshold for bottom extreme (default: 0.20).
#' @return List with integer 'tier' (1, 2, or 0), logical 'valid', integer 'n_pairs', and overall 'dunn_min'.
evaluate_lean_contrast_selection <- function(trait_vec,
                                            tree,
                                            target_pairs,
                                            discrete_method = "quintile",
                                            top_quantile = 0.80,
                                            bottom_quantile = 0.20) {
  
  if (target_pairs <= 0) {
    return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
  }
  
  sp_names <- names(trait_vec)
  if (is.null(sp_names) || length(sp_names) < 4) {
    return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
  }
  
  common_sp <- intersect(sp_names, tree$tip.label)
  if (length(common_sp) < 4) {
    return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
  }
  
  D <- cophenetic(tree)
  trait_vec <- trait_vec[common_sp]
  
  # Quantile thresholds
  g_med <- median(trait_vec, na.rm = TRUE)
  g_sd  <- sd(trait_vec, na.rm = TRUE)
  
  lower_thresh <- switch(discrete_method,
    "quartile"      = unname(quantile(trait_vec, 0.25, na.rm = TRUE)),
    "quintile"      = unname(quantile(trait_vec, 0.20, na.rm = TRUE)),
    "decile"        = unname(quantile(trait_vec, 0.10, na.rm = TRUE)),
    "median_sd"     = g_med - g_sd,
    "parameterized" = unname(quantile(trait_vec, bottom_quantile, na.rm = TRUE)),
    unname(quantile(trait_vec, 0.20, na.rm = TRUE))
  )
  
  upper_thresh <- switch(discrete_method,
    "quartile"      = unname(quantile(trait_vec, 0.75, na.rm = TRUE)),
    "quintile"      = unname(quantile(trait_vec, 0.80, na.rm = TRUE)),
    "decile"        = unname(quantile(trait_vec, 0.90, na.rm = TRUE)),
    "median_sd"     = g_med + g_sd,
    "parameterized" = unname(quantile(trait_vec, top_quantile, na.rm = TRUE)),
    unname(quantile(trait_vec, 0.80, na.rm = TRUE))
  )
  
  # Partition tips into high_extreme vs low_extreme
  high_extreme_sp <- names(trait_vec)[trait_vec >= upper_thresh & trait_vec > g_med]
  low_extreme_sp  <- names(trait_vec)[trait_vec <= lower_thresh & trait_vec < g_med]
  
  if (length(high_extreme_sp) == 0 || length(low_extreme_sp) == 0) {
    return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
  }
  
  # Build candidate pairs
  pairs_list <- list()
  idx <- 1
  for (h_s in high_extreme_sp) {
    for (l_s in low_extreme_sp) {
      if (h_s == l_s) next
      abs_diff <- abs(trait_vec[h_s] - trait_vec[l_s])
      if (abs_diff > 0) {
        pairs_list[[idx]] <- data.frame(
          sp1 = h_s,
          sp2 = l_s,
          dist = as.numeric(D[h_s, l_s]),
          abs_diff = as.numeric(abs_diff),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }
  
  pairs_df <- do.call(rbind, pairs_list)
  if (is.null(pairs_df) || nrow(pairs_df) == 0) {
    return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
  }
  
  # Sort by patristic distance (asc) then trait difference (desc)
  pairs_df <- pairs_df[order(pairs_df$dist, -pairs_df$abs_diff), ]
  
  # Function to try selecting k non-overlapping pairs and evaluating overall Dunn
  extract_k_pairs <- function(k_target) {
    selected_sp <- data.frame(sp = character(0), cluster = integer(0), stringsAsFactors = FALSE)
    n_cl <- 0L
    for (r in seq_len(nrow(pairs_df))) {
      q1 <- pairs_df$sp1[r]; q2 <- pairs_df$sp2[r]
      if (q1 %in% selected_sp$sp || q2 %in% selected_sp$sp) next
      
      selected_sp <- rbind(
        selected_sp,
        data.frame(sp = c(q1, q2), cluster = c(n_cl + 1L, n_cl + 1L), stringsAsFactors = FALSE)
      )
      n_cl <- n_cl + 1L
      if (n_cl == k_target) break
    }
    
    if (n_cl < k_target) return(list(valid = FALSE, dunn_min = 0, pairs = NULL))
    
    dsub <- D[selected_sp$sp, selected_sp$sp]
    cl_vec <- setNames(selected_sp$cluster, selected_sp$sp)
    dunn_min <- calculate_overall_dunn(dsub, cl_vec)
    
    list(valid = (dunn_min >= 1.0), dunn_min = dunn_min, pairs = selected_sp)
  }
  
  # Check Tier 1: Exactly target_pairs with overall mod_dunn >= 1.0
  res_t1 <- extract_k_pairs(target_pairs)
  if (res_t1$valid) {
    return(list(
      tier = 1L,
      valid = TRUE,
      n_pairs = target_pairs,
      dunn_min = res_t1$dunn_min,
      pairs = res_t1$pairs
    ))
  }
  
  # Check Tier 2: Fallback target_pairs - 1 (if target_pairs > 1)
  if (target_pairs > 1) {
    res_t2 <- extract_k_pairs(target_pairs - 1L)
    if (res_t2$valid) {
      return(list(
        tier = 2L,
        valid = TRUE,
        n_pairs = target_pairs - 1L,
        dunn_min = res_t2$dunn_min,
        pairs = res_t2$pairs
      ))
    }
  }
  
  # Invalid
  return(list(tier = 0L, valid = FALSE, n_pairs = 0L, dunn_min = 0, pairs = NULL))
}
