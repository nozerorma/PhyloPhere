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
#
# Candidate construction and the greedy add-loop are vectorized rather than
# scalar R `for` loops: candidate count grows as (percentile_frac * n_tips)^2,
# so a scalar loop is fine at ~31 tips (primates) but collapses at ~150+ tips
# (e.g. Carn: 1,024 candidates vs ~42, a 24x jump that measured as a 25x drop
# in draws/s — the dominant cost, not the BM/lambda simulation itself, which
# barely slows down over the same range). The vectorized form computes the
# SAME quantities via matrix indexing instead of per-candidate function calls;
# see the tie-break note by the greedy loop for why this is bit-identical to
# the original scalar implementation, not just equivalent in expectation.
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

# Is this permulated vector an ordinal fg/bg code (not a continuous measure)?
# Mirrors stats.R::is_ordinal_trait so the permulation applies the SAME rule the
# observed selection used. `trait_type`: "ordinal"/"continuous" force it,
# ""/"auto" infers from 2-5 all-integer levels.
lean_is_ordinal <- function(trait_vec, trait_type = "auto") {
  tt <- tolower(as.character(trait_type)[1])
  if (is.na(tt) || !nzchar(tt)) tt <- "auto"
  if (tt == "ordinal") return(TRUE)
  if (tt == "continuous") return(FALSE)
  u <- unique(trait_vec[!is.na(trait_vec)])
  length(u) >= 2 && length(u) <= 5 && all(u == round(u))
}

# Resolve the low/high trait cut-offs for a discretisation method.
# Mirrors the discrete_method / top_quantile / bottom_quantile semantics used by
# the contrast selection step (conf/common.config). For an ordinal fg/bg code
# the cut-offs are the extreme levels themselves (no quantiles).
lean_thresholds <- function(trait_vec, discrete_method, bottom_quantile, top_quantile,
                            trait_type = "auto") {
  g_med <- median(trait_vec, na.rm = TRUE)
  g_sd  <- stats::sd(trait_vec, na.rm = TRUE)
  q <- function(p) unname(stats::quantile(trait_vec, p, na.rm = TRUE))
  if (lean_is_ordinal(trait_vec, trait_type)) {
    lv <- sort(unique(trait_vec[!is.na(trait_vec)]))
    return(list(lower = lv[1], upper = lv[length(lv)], med = g_med))
  }
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
#' @param ci_lb,ci_ub     Optional named numeric vectors of per-tip confidence
#'                        bounds for THIS permulated vector. When supplied the
#'                        selector uses the count-data rule of
#'                        TRAIT_ANALYSIS/3.CI-composition.Rmd — a pair is a valid
#'                        candidate iff its intervals do not overlap — instead of
#'                        the percentile rule. Traits with no n/c columns (the
#'                        categorical path) simply omit these and keep the
#'                        percentile behaviour, which is itself the faithful
#'                        replication of that path's `trait_overlap`.
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
                                             top_quantile = 0.90,
                                             ci_lb = NULL,
                                             ci_ub = NULL,
                                             trait_type = "auto") {

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
    # Count-data path: candidates are ALL ordered pairs whose confidence
    # intervals do not overlap, with no percentile restriction — matching
    # create_pairwise_data()'s trait_overlap = ci_overlap(...) branch.
    lb <- ci_lb[sp]; ub <- ci_ub[sp]
    ok <- is.finite(lb) & is.finite(ub)
    if (sum(ok) < 2L * target_pairs) return(reject("too few species with usable CIs"))
    sp <- sp[ok]; trait_vec <- trait_vec[sp]; lb <- lb[sp]; ub <- ub[sp]
    high_sp <- low_sp <- sp
  } else {
    th <- lean_thresholds(trait_vec, discrete_method, bottom_quantile, top_quantile, trait_type)
    # The `>= th$med` / `<= th$med` guard (was `>` / `<`) only matters when the
    # median coincides with a threshold, which is exactly the binary / bimodal
    # case: a minority-foreground 0/1 vector has median == lower == 0, so the
    # strict `< th$med` made `low_sp` empty and the permulation pool unfillable.
    # The production discrete path (3.CI-composition.Rmd) has no median gate at
    # all — it categorises purely on `trait >= upper_thresh` / `<= lower_thresh` —
    # so relaxing to non-strict here brings the two into line and is a no-op for
    # every continuous trait (upper_thresh >= median >= lower_thresh always).
    high_sp <- names(trait_vec)[trait_vec >= th$upper & trait_vec >= th$med]
    low_sp  <- names(trait_vec)[trait_vec <= th$lower & trait_vec <= th$med]
    if (length(high_sp) < target_pairs || length(low_sp) < target_pairs) {
      return(reject("not enough extreme species to form target_pairs"))
    }
  }

  # ---- candidate pairs (vectorized) ------------------------------------------
  # rep(..., each=)/rep(..., times=) reproduces `for (h in high_sp) for (l in
  # low_sp)` row-for-row (h varies slowest) — kept identical to the original
  # generation order so that if two candidates ever tie exactly on both sort
  # keys below, the tie-break (stable order()) resolves the same way.
  c_hi_all <- rep(high_sp, each = length(low_sp))
  c_lo_all <- rep(low_sp, times = length(high_sp))
  dif_all  <- trait_vec[c_hi_all] - trait_vec[c_lo_all]
  # Orientation matters: the high-extreme member must actually be the higher one.
  keep <- (c_hi_all != c_lo_all) & (dif_all > 0)
  if (use_ci) {
    # ci_overlap(lb1, ub1, lb2, ub2) = (lb1 <= ub2) & (lb2 <= ub1); keep only
    # NON-overlapping pairs, exactly as pair_sel.f's input filter does.
    overlap <- (lb[c_hi_all] <= ub[c_lo_all]) & (lb[c_lo_all] <= ub[c_hi_all])
    keep <- keep & !overlap
  }
  if (!any(keep)) {
    return(reject(if (use_ci) "no pair with non-overlapping CIs"
                  else "no candidate pair with non-zero trait difference"))
  }
  c_hi  <- c_hi_all[keep]; c_lo <- c_lo_all[keep]; c_dif <- dif_all[keep]
  c_dist <- D[cbind(c_hi, c_lo)]   # name-pair matrix indexing, vectorized

  # Production order: smallest patristic distance first, largest trait diff to break ties.
  ord    <- order(c_dist, -c_dif)
  c_hi   <- c_hi[ord]; c_lo <- c_lo[ord]; c_dist <- c_dist[ord]; c_dif <- c_dif[ord]

  # ---- seed with the top-ranked pair (pair_sel.f lines 235-249) --------------
  members  <- list(c(c_hi[1], c_lo[1]))
  sel_fg   <- c_hi[1]; sel_bg <- c_lo[1]
  used     <- c(c_hi[1], c_lo[1])
  used_flag <- setNames(logical(length(unique(c(c_hi, c_lo)))), unique(c(c_hi, c_lo)))
  used_flag[used] <- TRUE

  # ---- greedily add the HIGHEST-Dunn candidate until target_pairs is reached --
  # No Dunn gate during construction: the pair COUNT is the hard constraint (a
  # permulated design must contain as many species as the observed one), and how
  # independent the assembled set turned out is scored afterwards as the tier.
  #
  # Vectorized re-derivation of mod_dunn_lean(D, c(members, candidate), last) for
  # every remaining candidate at once: intra of a 2-member cluster is just its
  # own hi-lo distance (== c_dist, already computed); inter is the closest
  # distance from either member to ANY already-placed species. min() is exactly
  # associative under IEEE754 (a pure comparison, no rounding), so taking the
  # min over the union `used` in one shot equals mod_dunn_lean's min-of-per-
  # cluster-mins — same bits, not merely the same value in expectation. That, plus
  # generating candidates in the original nested-loop order, is what makes
  # order(-dv, -dif) select the identical candidate order(): a stable sort
  # breaks exact ties in favour of the earlier (dist-ascending) entry, same as
  # the scalar loop only ever overwriting `best` on a strict improvement.
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
       # In CI mode there is no percentile cut-off to audit against; the auditable
       # invariant is instead that every emitted pair had non-overlapping intervals.
       q_lower = if (use_ci) NA_real_ else th$lower,
       q_upper = if (use_ci) NA_real_ else th$upper,
       mode = if (use_ci) "ci" else "percentile",
       reason = "accepted")
}
