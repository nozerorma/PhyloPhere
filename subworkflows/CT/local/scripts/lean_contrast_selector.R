#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/lean_contrast_selector.R
# =============================================================================
# Lean port of TRAIT_ANALYSIS/local/src/selection_algorithm.R::pair_sel.f, for
# use inside the permulation harvesting loop.
#
# The production selector picks contrast pairs from CI-non-overlap candidates;
# this lean version substitutes percentile discretisation (there are no
# confidence intervals for a simulated trait vector) but keeps the selection
# ALGORITHM identical: seed with the closest/most-divergent pair, then greedily
# add the candidate with the HIGHEST modified Dunn index, gated at Dunn >= 1.
#
# Selecting by highest Dunn (rather than by smallest patristic distance) is what
# makes an accepted permulation comparable to the observed design: both were
# produced by the same rule.
#
# Everything is base R over a precomputed distance matrix — no data frames in
# the inner loop — because this runs once per permulation draw (up to max_tries).
# =============================================================================

# Modified Dunn index for one cluster, verbatim semantics of
# selection_algorithm.R::mod_dunn: (min distance to any other cluster) divided
# by (this cluster's own diameter).
#
# @param D        Full patristic distance matrix (dimnames = species).
# @param members  List of character vectors, one per cluster.
# @param k        Index of the cluster to score.
mod_dunn_lean <- function(D, members, k) {
  c1 <- members[[k]]
  intra <- if (length(c1) > 1) max(D[c1, c1]) else 0
  if (intra == 0) return(Inf)
  inter <- Inf
  for (j in seq_along(members)) {
    if (j == k) next
    inter <- min(inter, min(D[c1, members[[j]]]))
  }
  if (!is.finite(inter)) return(Inf)
  inter / intra
}

# Overall Dunn across every cluster = min_k mod_dunn(k). This is the validation
# condition: all clusters mutually independent iff overall >= 1.
#
# Recomputed on the FINAL set rather than trusted from the incremental gate,
# because adding cluster k can only shrink an earlier cluster's inter-cluster
# minimum — so a set built from individually-passing pairs can still end below 1.
overall_dunn_lean <- function(D, members) {
  if (length(members) <= 1) return(Inf)
  min(vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1)))
}

# Resolve the low/high trait cut-offs for a discretisation method.
# Mirrors the discrete_method / top_quantile / bottom_quantile semantics used by
# the contrast selection step (conf/common.config).
lean_thresholds <- function(trait_vec, discrete_method, bottom_quantile, top_quantile) {
  g_med <- median(trait_vec, na.rm = TRUE)
  g_sd  <- stats::sd(trait_vec, na.rm = TRUE)
  q <- function(p) unname(stats::quantile(trait_vec, p, na.rm = TRUE))
  switch(discrete_method,
    "quartile"      = list(lower = q(0.25), upper = q(0.75), med = g_med),
    "quintile"      = list(lower = q(0.20), upper = q(0.80), med = g_med),
    "decile"        = list(lower = q(0.10), upper = q(0.90), med = g_med),
    "median_sd"     = list(lower = g_med - g_sd, upper = g_med + g_sd, med = g_med),
    "parameterized" = list(lower = q(bottom_quantile), upper = q(top_quantile), med = g_med),
    list(lower = q(0.20), upper = q(0.80), med = g_med)   # default: quintile
  )
}

#' Lean contrast selection + tiered Dunn validation for one permulated vector.
#'
#' @param trait_vec       Named numeric vector of permulated trait values.
#' @param D               Precomputed patristic distance matrix on the REAL
#'                        (unrescaled) tree. Hoisted by the caller — never
#'                        recompute cophenetic() inside the draw loop.
#' @param target_pairs    N_pairs_obs, the observed independent pair count.
#' @param discrete_method "quintile" (default), "quartile", "decile",
#'                        "median_sd", or "parameterized" — as in common.config.
#' @param bottom_quantile Lower cut-off when method = "parameterized".
#' @param top_quantile    Upper cut-off when method = "parameterized".
#' @return list(tier, n_pairs, dunn_min, n_below, fg, bg, reason)
#'   Both accepted tiers carry EXACTLY target_pairs pairs — the permulated design
#'   must contain the same number of species as the observed one. The tiers differ
#'   only in how strictly independence is enforced:
#'   tier 1 = every pair has mod_dunn >= 1 (fully independent, like the observed)
#'   tier 2 = exactly one pair falls below mod_dunn 1 (relaxed fallback)
#'   tier 0 = rejected: could not form target_pairs pairs, or two or more pairs
#'            fall below 1. `reason` says which.
#'   fg/bg  = high-extreme and low-extreme member of each pair, in matched order
evaluate_lean_contrast_selection <- function(trait_vec,
                                             D,
                                             target_pairs,
                                             discrete_method = "quintile",
                                             bottom_quantile = 0.10,
                                             top_quantile = 0.90) {

  reject <- function(reason, n_pairs = 0L, dunn = 0, n_below = NA_integer_) {
    list(tier = 0L, n_pairs = n_pairs, dunn_min = dunn, n_below = n_below,
         fg = NULL, bg = NULL, reason = reason)
  }

  if (target_pairs <= 0L) return(reject("target_pairs <= 0"))

  sp <- intersect(names(trait_vec), rownames(D))
  if (length(sp) < 2L * target_pairs) return(reject("too few species with distances"))
  trait_vec <- trait_vec[sp]

  th <- lean_thresholds(trait_vec, discrete_method, bottom_quantile, top_quantile)
  high_sp <- names(trait_vec)[trait_vec >= th$upper & trait_vec > th$med]
  low_sp  <- names(trait_vec)[trait_vec <= th$lower & trait_vec < th$med]
  if (length(high_sp) < target_pairs || length(low_sp) < target_pairs) {
    return(reject("not enough extreme species to form target_pairs"))
  }

  # ---- candidate pairs: every high x low combination -------------------------
  n_cand <- length(high_sp) * length(low_sp)
  c_hi   <- character(n_cand); c_lo <- character(n_cand)
  c_dist <- numeric(n_cand);   c_dif <- numeric(n_cand)
  i <- 0L
  for (h in high_sp) for (l in low_sp) {
    d <- abs(trait_vec[[h]] - trait_vec[[l]])
    if (d <= 0) next
    i <- i + 1L
    c_hi[i] <- h; c_lo[i] <- l; c_dist[i] <- D[h, l]; c_dif[i] <- d
  }
  if (i == 0L) return(reject("no candidate pair with non-zero trait difference"))
  keep   <- seq_len(i)
  c_hi   <- c_hi[keep]; c_lo <- c_lo[keep]; c_dist <- c_dist[keep]; c_dif <- c_dif[keep]

  # Production order: smallest patristic distance first, largest trait diff to break ties.
  ord    <- order(c_dist, -c_dif)
  c_hi   <- c_hi[ord]; c_lo <- c_lo[ord]; c_dist <- c_dist[ord]; c_dif <- c_dif[ord]

  # ---- seed with the top-ranked pair (pair_sel.f lines 235-249) --------------
  members  <- list(c(c_hi[1], c_lo[1]))
  sel_fg   <- c_hi[1]; sel_bg <- c_lo[1]
  used     <- c(c_hi[1], c_lo[1])

  # ---- greedily add the HIGHEST-Dunn candidate until target_pairs is reached --
  # No Dunn gate during construction: the pair COUNT is the hard constraint (a
  # permulated design must contain as many species as the observed one), and how
  # independent the assembled set turned out is scored afterwards as the tier.
  while (length(members) < target_pairs) {
    best_d <- -Inf; best_i <- NA_integer_
    for (j in seq_along(c_hi)) {
      if (c_hi[j] %in% used || c_lo[j] %in% used) next
      cand <- c(members, list(c(c_hi[j], c_lo[j])))
      dv   <- mod_dunn_lean(D, cand, length(cand))
      # ties broken by larger trait difference, matching the production arrange()
      if (dv > best_d || (dv == best_d && !is.na(best_i) && c_dif[j] > c_dif[best_i])) {
        best_d <- dv; best_i <- j
      }
    }
    if (is.na(best_i)) break                   # no non-overlapping candidate left
    members <- c(members, list(c(c_hi[best_i], c_lo[best_i])))
    sel_fg  <- c(sel_fg, c_hi[best_i]); sel_bg <- c(sel_bg, c_lo[best_i])
    used    <- c(used, c_hi[best_i], c_lo[best_i])
  }

  n_pairs <- length(members)
  if (n_pairs < target_pairs) {
    return(reject("could not form target_pairs non-overlapping pairs", n_pairs))
  }

  # Score independence on the FINAL set: adding a cluster can only shrink an
  # earlier cluster's inter-cluster minimum, so per-cluster Dunn must be
  # recomputed once the set is complete rather than trusted from build time.
  dunn_vec <- vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1))
  n_below  <- sum(dunn_vec < 1)
  dunn_min <- min(dunn_vec)

  if (n_below >= 2L) {
    return(reject("two or more pairs below Dunn 1", n_pairs, dunn_min, n_below))
  }

  # The permulated trait values of the chosen species, plus the thresholds they
  # were selected against, are returned so a run can be audited after the fact:
  # every FG value must be >= upper and every BG value <= lower. Without this
  # there is no way to confirm from the output alone that the written species
  # really were the extremes of their own permulated vector.
  list(tier = if (n_below == 0L) 1L else 2L,
       n_pairs = n_pairs, dunn_min = dunn_min, n_below = n_below,
       fg = sel_fg, bg = sel_bg,
       fg_values = unname(trait_vec[sel_fg]), bg_values = unname(trait_vec[sel_bg]),
       q_lower = th$lower, q_upper = th$upper, reason = "accepted")
}
