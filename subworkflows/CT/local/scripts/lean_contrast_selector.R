#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/lean_contrast_selector.R
# =============================================================================
# Lean port of contrast selection for use inside the permulation harvesting loop.
#
# For count data (with confidence bounds), candidates are filtered by CI non-overlap.
# For continuous traits without CIs, candidates are filtered by the Global Top 1%
# Phylogenetic Shift Score (PSS) under the fitted evolutionary model (OU / BM).
#
# Selection algorithm: seed with the highest-ranked candidate pair, then greedily
# add candidate pairs that maximize the modified Dunn index, enforcing independence
# across clusters (mod_dunn >= 1.0).
# =============================================================================

# Modified Dunn index for one cluster, verbatim semantics of
# selection_algorithm.R::mod_dunn: (min distance to any other cluster) divided
# by (this cluster's own diameter).
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

# Overall Dunn across every cluster = min_k mod_dunn(k).
overall_dunn_lean <- function(D, members) {
  if (length(members) <= 1) return(Inf)
  min(vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1)))
}

# Analytical shift score S and composite PSS score
calc_pss_scores <- function(hi, lo, dif, dist, C_mat) {
  pv <- diag(C_mat)[hi] + diag(C_mat)[lo] - 2 * C_mat[cbind(hi, lo)]
  tolerance <- 100 * .Machine$double.eps * max(abs(pv), 1)
  pv[pv < 0 & pv >= -tolerance] <- 0
  z <- dif / sqrt(pmax(pv, 1e-12))
  S <- -expm1(log(2) + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
  max_dif <- max(dif, na.rm = TRUE)
  max_dst <- max(dist, na.rm = TRUE)
  if (max_dif <= 0) max_dif <- 1
  if (max_dst <= 0) max_dst <- 1
  S * (dif / max_dif) / (dist / max_dst)
}

#' Lean contrast selection + tiered Dunn validation for one permulated vector.
#'
#' @param trait_vec       Named numeric vector of permulated trait values.
#' @param D               Precomputed patristic distance matrix on the REAL tree.
#' @param target_pairs    N_pairs_obs, the observed independent pair count.
#' @param ci_lb,ci_ub     Optional named numeric vectors of per-tip confidence bounds.
#' @param cov_matrix      Optional phylogenetic covariance matrix (OU / BM) for PSS calculation.
#' @param top_pct         Top percentile of PSS scores to retain (default: 0.01 for Top 1%).
#' @return list(tier, n_pairs, dunn_min, n_below, fg, bg, reason)
evaluate_lean_contrast_selection <- function(trait_vec,
                                             D,
                                             target_pairs,
                                             ci_lb = NULL,
                                             ci_ub = NULL,
                                             cov_matrix = NULL,
                                             top_pct = 0.01) {

  reject <- function(reason, n_pairs = 0L, dunn = 0, n_below = NA_integer_) {
    list(tier = 0L, n_pairs = n_pairs, dunn_min = dunn, n_below = n_below,
         fg = NULL, bg = NULL, reason = reason)
  }

  if (target_pairs <= 0L) return(reject("target_pairs <= 0"))

  sp <- intersect(names(trait_vec), rownames(D))
  if (length(sp) < 2L * target_pairs) return(reject("too few species with distances"))
  trait_vec <- trait_vec[sp]

  use_ci <- !is.null(ci_lb) && !is.null(ci_ub)

  if (use_ci) {
    lb <- ci_lb[sp]; ub <- ci_ub[sp]
    ok <- is.finite(lb) & is.finite(ub)
    if (sum(ok) < 2L * target_pairs) return(reject("too few species with usable CIs"))
    sp <- sp[ok]; trait_vec <- trait_vec[sp]; lb <- lb[sp]; ub <- ub[sp]
  }

  # Build all ordered pairwise candidates
  c_hi_all <- rep(sp, each = length(sp))
  c_lo_all <- rep(sp, times = length(sp))
  dif_all  <- trait_vec[c_hi_all] - trait_vec[c_lo_all]
  valid    <- (c_hi_all != c_lo_all) & (dif_all > 0)

  if (!any(valid)) return(reject("no candidate pairs with positive trait difference"))

  c_hi   <- c_hi_all[valid]
  c_lo   <- c_lo_all[valid]
  c_dif  <- dif_all[valid]
  c_dist <- D[cbind(c_hi, c_lo)]

  if (use_ci) {
    # Count data: candidate gate is non-overlapping Bayesian CIs
    ci_nonoverlap <- (lb[c_hi] > ub[c_lo])
    if (!any(ci_nonoverlap)) return(reject("no pair with non-overlapping CIs"))
    c_hi   <- c_hi[ci_nonoverlap]
    c_lo   <- c_lo[ci_nonoverlap]
    c_dif  <- c_dif[ci_nonoverlap]
    c_dist <- c_dist[ci_nonoverlap]
    c_score <- -c_dist + (c_dif / max(c_dif)) * 1e-4
  } else if (!is.null(cov_matrix)) {
    # Continuous trait: candidate gate is Global Top 1% PSS under the evolutionary model
    pss_scores <- calc_pss_scores(c_hi, c_lo, c_dif, c_dist, cov_matrix)
    n_keep <- max(8L, ceiling(length(pss_scores) * top_pct))
    thresh <- sort(pss_scores, decreasing = TRUE)[n_keep]
    top_mask <- (pss_scores >= thresh)
    if (!any(top_mask)) return(reject("no candidate pairs in top PSS slice"))
    c_hi   <- c_hi[top_mask]
    c_lo   <- c_lo[top_mask]
    c_dif  <- c_dif[top_mask]
    c_dist <- c_dist[top_mask]
    c_score <- pss_scores[top_mask]
  } else {
    # Fallback when no covariance matrix is provided
    c_score <- -c_dist + (c_dif / max(c_dif)) * 1e-4
  }

  # Sort candidate pairs by score descending
  ord    <- order(-c_score)
  c_hi   <- c_hi[ord]; c_lo <- c_lo[ord]; c_dist <- c_dist[ord]; c_dif <- c_dif[ord]; c_score <- c_score[ord]

  # Seed with top-ranked candidate pair
  members   <- list(c(c_hi[1], c_lo[1]))
  sel_fg    <- c_hi[1]; sel_bg <- c_lo[1]
  used      <- c(c_hi[1], c_lo[1])
  used_flag <- setNames(logical(length(unique(c(c_hi, c_lo)))), unique(c(c_hi, c_lo)))
  used_flag[used] <- TRUE

  # Greedily add the HIGHEST-Dunn candidate until target_pairs is reached
  while (length(members) < target_pairs) {
    avail <- !used_flag[c_hi] & !used_flag[c_lo]
    if (!any(avail)) break
    idx   <- which(avail)
    hi_a  <- c_hi[idx]; lo_a <- c_lo[idx]; dist_a <- c_dist[idx]; dif_a <- c_dif[idx]

    D_hi <- D[hi_a, used, drop = FALSE]
    D_lo <- D[lo_a, used, drop = FALSE]
    inter_hi <- D_hi[, 1]; if (ncol(D_hi) > 1) for (cc in 2:ncol(D_hi)) inter_hi <- pmin(inter_hi, D_hi[, cc])
    inter_lo <- D_lo[, 1]; if (ncol(D_lo) > 1) for (cc in 2:ncol(D_lo)) inter_lo <- pmin(inter_lo, D_lo[, cc])
    inter <- pmin(inter_hi, inter_lo)
    dv    <- ifelse(dist_a == 0, Inf, inter / dist_a)

    best_local <- order(-dv, -dif_a)[1]
    best_i <- idx[best_local]

    members <- c(members, list(c(c_hi[best_i], c_lo[best_i])))
    sel_fg  <- c(sel_fg, c_hi[best_i]); sel_bg <- c(sel_bg, c_lo[best_i])
    used    <- c(used, c_hi[best_i], c_lo[best_i])
    used_flag[c(c_hi[best_i], c_lo[best_i])] <- TRUE
  }

  n_pairs <- length(members)
  if (n_pairs < target_pairs) {
    return(reject("could not form target_pairs non-overlapping pairs", n_pairs))
  }

  # Score independence on the FINAL set
  dunn_vec <- vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1))
  n_below  <- sum(dunn_vec < 1)
  dunn_min <- min(dunn_vec)

  if (n_below >= 2L) {
    return(reject("two or more pairs below Dunn 1", n_pairs, dunn_min, n_below))
  }

  list(tier = if (n_below == 0L) 1L else 2L,
       n_pairs = n_pairs, dunn_min = dunn_min, n_below = n_below,
       fg = sel_fg, bg = sel_bg,
       fg_values = unname(trait_vec[sel_fg]), bg_values = unname(trait_vec[sel_bg]),
       mean_pd = mean(vapply(members, function(x) D[x[1], x[2]], numeric(1))),
       mean_df = mean(vapply(members, function(x) abs(trait_vec[x[1]] - trait_vec[x[2]]), numeric(1))),
       mode = if (use_ci) "ci" else "pss_top1pct",
       reason = "accepted")
}
