#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/lean_contrast_selector.R
# =============================================================================
# Lean port of contrast selection for use inside the permulation harvesting loop.
#
# For count data (with confidence bounds), candidates are filtered by CI non-overlap.
# For ordinal fg/bg-coded traits, candidates are max-level vs min-level pairs only.
# For continuous traits without CIs, candidates are filtered by the Global Top 1%
# Phylogenetic Shift Score (PSS) under the fitted evolutionary model (OU / BM).
# These mirror the three branches of 3.CI-composition.Rmd exactly so the
# permulation null applies the same rule the observed selection used.
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

# PSS scoring is provided by the vendored phyloq engine (src/pss_core.R):
# analytical_s() and calculate_pairwise_scores(). This file only needs to be
# sourced alongside it — permulations.R and the report modules stage both.
if (!exists("calculate_pairwise_scores", mode = "function")) {
  .pss_core_paths <- c(
    file.path(getwd(), "src", "pss_core.R"),
    file.path(getwd(), "pss_core.R"),
    "/home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/CT/local/scripts/pss_core.R"
  )
  .hit <- .pss_core_paths[file.exists(.pss_core_paths)]
  if (length(.hit)) source(.hit[1])
}

# Ordinal fg/bg code auto-detection (mirrors stats.R::is_ordinal_trait "auto").
.is_ordinal_vec <- function(v) {
  u <- unique(v[!is.na(v)])
  length(u) >= 2 && length(u) <= 5 && all(u == round(u))
}

# ---------------------------------------------------------------------------
# SHARED CORE — used by both the observed selector (selection_algorithm.R::
# pair_sel.f) and the permulation null (evaluate_lean_contrast_selection below)
# so the two apply an identical candidate-ranking + greedy-Dunn rule.
# ---------------------------------------------------------------------------

#' Unified candidate ranking policy.
#'
#' @param df data.frame of candidate pairs. Required: species1, species2,
#'   distance, abs_diff. Optional: pss_score (divergence signal), pair_n
#'   (combined sample size, count data only).
#' @return df row-reordered, best candidate first.
#'
#' Primary key: PSS score descending when any finite pss_score is present
#' (continuous, count and ordinal all get an OU/BM PSS); patristic distance
#' ascending otherwise (PSS fit failed). Ties: |trait difference| desc, then
#' combined pair sample size desc when available.
rank_candidates <- function(df) {
  if (nrow(df) == 0) return(df)
  has_pss <- "pss_score" %in% names(df) && any(is.finite(df$pss_score))
  keys <- if (has_pss) list(-df$pss_score) else list(df$distance)
  keys <- c(keys, list(-df$abs_diff))
  if ("pair_n" %in% names(df)) keys <- c(keys, list(-df$pair_n))
  keys <- c(keys, list(seq_len(nrow(df))))  # stable: preserve incoming order on full ties
  df[do.call(order, keys), , drop = FALSE]
}

#' Unified greedy Dunn-gated pair assembly.
#'
#' @param ranked candidate pairs already ordered best-first by rank_candidates().
#'   Needs species1, species2, distance, abs_diff (+ pair_n optional).
#' @param D patristic distance matrix (species x species).
#' @param target stop after this many pairs (Inf = run until enforce_dunn stops it).
#' @param enforce_dunn TRUE  -> only accept a pair that keeps every cluster's
#'                              modified Dunn >= 1; stop when none qualifies
#'                              (observed selector: variable pair count).
#'                     FALSE -> always take the highest-Dunn available candidate
#'                              up to `target`, even below 1; caller grades the
#'                              result (permulation null: fixed N, tiered).
#' @return list(selected = data.frame of chosen rows + Dunn_index + cluster,
#'              members  = list of c(species1, species2)).
greedy_dunn_select <- function(ranked, D, target = Inf, enforce_dunn = TRUE) {
  D <- as.matrix(D)
  empty <- ranked[0, , drop = FALSE]
  if (nrow(ranked) == 0 || target < 1) return(list(selected = empty, members = list()))

  seed <- ranked[1, , drop = FALSE]
  seed$Dunn_index <- Inf
  seed$cluster    <- 1L
  selected <- seed
  members  <- list(c(seed$species1, seed$species2))
  used     <- c(seed$species1, seed$species2)

  while (length(members) < target) {
    avail <- !(ranked$species1 %in% used | ranked$species2 %in% used)
    if (!any(avail)) break
    cand <- ranked[avail, , drop = FALSE]

    dunn <- vapply(seq_len(nrow(cand)), function(i) {
      mod_dunn_lean(D, c(members, list(c(cand$species1[i], cand$species2[i]))), length(members) + 1L)
    }, numeric(1))
    cand$Dunn_index <- round(dunn, 4)

    if (enforce_dunn) {
      cand <- cand[cand$Dunn_index >= 1, , drop = FALSE]
      if (nrow(cand) == 0) break
    }

    # Highest Dunn; ties resolved by the incoming rank order (stable sort).
    cand <- cand[order(-cand$Dunn_index, seq_len(nrow(cand))), , drop = FALSE]
    best <- cand[1, , drop = FALSE]

    new_members <- c(members, list(c(best$species1, best$species2)))
    if (enforce_dunn && overall_dunn_lean(D, new_members) < 1) break

    best$cluster <- length(members) + 1L
    selected <- rbind(selected, best)
    members  <- new_members
    used     <- c(used, best$species1, best$species2)
  }

  list(selected = selected, members = members)
}

#' Candidate-pair gate + ranking for one (permulated) trait vector.
#'
#' Shared by evaluate_lean_contrast_selection() (the tiered Dunn null) and
#' lean_fop_harvest() (the per-cycle FOP mirror). Reproduces stages 1-2 of the
#' observed selector: PSS via the vendored phyloq engine on the fixed observed-
#' model covariances, then the trait-type non-overlap gate (CI / ordinal /
#' continuous top_pct), then rank_candidates().
#'
#' @return list(cand_df = ranked data.frame | NULL, mode = "ci"|"ordinal"|"pss",
#'              reason = NULL | character). cand_df carries species1, species2,
#'              distance, abs_diff, pss_score (+ pair_n when n_vec given).
lean_candidate_df <- function(trait_vec, D, target_pairs,
                              tree = NULL, cov_bm = NULL, cov_ou = NULL,
                              selected_model = "OU", ci_lb = NULL, ci_ub = NULL,
                              top_pct = 0.01, ordinal = NULL, n_vec = NULL) {
  D <- as.matrix(D)
  sp <- intersect(names(trait_vec), rownames(D))
  if (length(sp) < 2L * target_pairs) return(list(cand_df = NULL, mode = "na", reason = "too few species with distances"))
  trait_vec <- trait_vec[sp]

  use_ci <- !is.null(ci_lb) && !is.null(ci_ub)
  if (is.null(ordinal)) ordinal <- !use_ci && .is_ordinal_vec(trait_vec)
  if (use_ci) {
    lb <- ci_lb[sp]; ub <- ci_ub[sp]
    ok <- is.finite(lb) & is.finite(ub)
    if (sum(ok) < 2L * target_pairs) return(list(cand_df = NULL, mode = "ci", reason = "too few species with usable CIs"))
    sp <- sp[ok]; trait_vec <- trait_vec[sp]; lb <- lb[sp]; ub <- ub[sp]
  }
  if (is.null(tree) || is.null(cov_bm) || is.null(cov_ou)) {
    return(list(cand_df = NULL, mode = "na", reason = "PSS inputs missing (tree / cov_bm / cov_ou)"))
  }
  tr <- if (length(sp) < length(tree$tip.label)) ape::drop.tip(tree, setdiff(tree$tip.label, sp)) else tree
  tr_sp <- tr$tip.label
  sc <- calculate_pairwise_scores(trait_vec[tr_sp], tr,
                                  cov_bm[tr_sp, tr_sp, drop = FALSE],
                                  cov_ou[tr_sp, tr_sp, drop = FALSE], selected_model)
  hi_is_1 <- sc$TraitValue1 >= sc$TraitValue2
  c_hi  <- ifelse(hi_is_1, sc$Species1, sc$Species2)
  c_lo  <- ifelse(hi_is_1, sc$Species2, sc$Species1)
  c_dif <- abs(sc$TraitValue1 - sc$TraitValue2)
  c_dist <- sc$PatristicDistance
  pss   <- sc$FinalScore
  drop0 <- c_dif > 0
  if (!any(drop0)) return(list(cand_df = NULL, mode = if (use_ci) "ci" else if (isTRUE(ordinal)) "ordinal" else "pss",
                               reason = "no candidate pairs with positive trait difference"))
  c_hi <- c_hi[drop0]; c_lo <- c_lo[drop0]; c_dif <- c_dif[drop0]
  c_dist <- c_dist[drop0]; pss <- pss[drop0]

  if (use_ci) {
    type_mask <- (lb[c_hi] > ub[c_lo])
  } else if (isTRUE(ordinal)) {
    lv_hi <- max(trait_vec, na.rm = TRUE); lv_lo <- min(trait_vec, na.rm = TRUE)
    type_mask <- (trait_vec[c_hi] >= lv_hi) & (trait_vec[c_lo] <= lv_lo)
  } else {
    type_mask <- rep(TRUE, length(c_hi))
  }
  if (!any(type_mask)) return(list(cand_df = NULL, mode = if (use_ci) "ci" else if (isTRUE(ordinal)) "ordinal" else "pss",
                                   reason = "no pair passes the trait-type gate"))
  s1 <- which(type_mask)
  if (!use_ci && !isTRUE(ordinal)) {
    n_keep <- min(max(1L, ceiling(length(s1) * top_pct)), length(s1))
    keep_i <- s1[order(pss[s1], decreasing = TRUE)][seq_len(n_keep)]
  } else {
    keep_i <- s1
  }
  keep <- logical(length(c_hi)); keep[keep_i] <- TRUE

  cand_df <- data.frame(
    species1 = c_hi[keep], species2 = c_lo[keep],
    distance = c_dist[keep], abs_diff = c_dif[keep],
    pss_score = pss[keep], stringsAsFactors = FALSE
  )
  if (!is.null(n_vec)) {
    cand_df$pair_n <- as.numeric(n_vec[cand_df$species1]) + as.numeric(n_vec[cand_df$species2])
  }
  list(cand_df = rank_candidates(cand_df),
       mode = if (use_ci) "ci" else if (isTRUE(ordinal)) "ordinal" else "pss",
       reason = NULL)
}

#' Per-cycle FOP hypothesis harvest — the null mirror of
#' selection_algorithm.R::fop_pair_sel.f steps 2-4.
#'
#' Same candidate gate/rank as the observed harvest and as
#' evaluate_lean_contrast_selection (shared lean_candidate_df). Canonical H1 =
#' greedy_dunn_select(enforce_dunn = TRUE); Voronoi-partition the species by
#' nearest canonical pair; draw one in-domain candidate per domain (exhaustive
#' <= ITER_CAP, else capped random); keep Dunn >= 1; rank by min-PSS and keep
#' the top (max_fop - 1). Returns per-hypothesis fg/bg species and per-domain
#' PSS so the permulation null domain-pools exactly as the observed §2b does.
#'
#' @return list(hypotheses = named list H1..Hn of data.frame(species1, species2,
#'              pss_score, cluster), K = int) or list(hypotheses = list(), K = 0).
#' @param canon_pairs optional data.frame(species1, species2) — the canonical
#'   H1 pairs from the null's own tiered fixed-N selection
#'   (evaluate_lean_contrast_selection). When supplied, H1 is taken verbatim and
#'   the Voronoi partition is seeded from it, so the FOP mirror harvests
#'   alternatives around the SAME canonical contrast the null already accepted
#'   (rather than re-running an enforce_dunn greedy that could diverge from it).
lean_fop_harvest <- function(trait_vec, D, target_pairs,
                             tree = NULL, cov_bm = NULL, cov_ou = NULL,
                             selected_model = "OU", ci_lb = NULL, ci_ub = NULL,
                             top_pct = 0.01, ordinal = NULL, n_vec = NULL,
                             max_fop = 100L, seed = 42L, canon_pairs = NULL) {
  set.seed(seed)
  cc <- lean_candidate_df(trait_vec, D, target_pairs, tree, cov_bm, cov_ou,
                          selected_model, ci_lb, ci_ub, top_pct, ordinal, n_vec)
  cand_df <- cc$cand_df
  if (is.null(cand_df) || nrow(cand_df) < 1) return(list(hypotheses = list(), K = 0L))

  Dm <- as.matrix(D)
  add_cluster <- function(df) { df$cluster <- seq_len(nrow(df)); df }

  # PSS lookup over the candidate pool (either orientation) for canon_pairs.
  .pss_of <- function(s1, s2) {
    hit <- which((cand_df$species1 == s1 & cand_df$species2 == s2) |
                 (cand_df$species1 == s2 & cand_df$species2 == s1))
    if (length(hit)) cand_df$pss_score[hit[1]] else NA_real_
  }

  if (!is.null(canon_pairs) && nrow(canon_pairs) >= 1) {
    cp <- data.frame(
      species1  = as.character(canon_pairs$species1),
      species2  = as.character(canon_pairs$species2),
      pss_score = mapply(.pss_of, canon_pairs$species1, canon_pairs$species2),
      stringsAsFactors = FALSE
    )
    canon <- list(selected = cp,
                  members = lapply(seq_len(nrow(cp)), function(i) c(cp$species1[i], cp$species2[i])))
  } else {
    canon <- greedy_dunn_select(cand_df, Dm, target = target_pairs, enforce_dunn = TRUE)
  }
  K <- length(canon$members)
  if (K < 1) return(list(hypotheses = list(), K = 0L))
  cm <- canon$members

  hyps <- list(H1 = add_cluster(canon$selected[, c("species1", "species2", "pss_score"), drop = FALSE]))
  if (K < target_pairs) return(list(hypotheses = hyps, K = K))  # H1 only, degenerate

  all_sp <- rownames(Dm)
  dom <- vapply(all_sp, function(s)
    which.min(vapply(seq_len(K), function(k) min(Dm[s, cm[[k]][1]], Dm[s, cm[[k]][2]]), numeric(1))),
    integer(1)); names(dom) <- all_sp
  pools <- lapply(seq_len(K), function(k) {
    ds <- names(dom)[dom == k]
    cand_df[cand_df$species1 %in% ds & cand_df$species2 %in% ds, , drop = FALSE]
  })
  psz <- vapply(pools, nrow, integer(1))
  if (any(psz == 0L)) return(list(hypotheses = hyps, K = K))

  ITER_CAP <- as.integer(max_fop) * 20L
  total <- prod(as.numeric(psz))
  idx <- if (total <= ITER_CAP)
    do.call(expand.grid, c(lapply(psz, seq_len), list(KEEP.OUT.ATTRS = FALSE)))
  else
    as.data.frame(lapply(psz, function(n) sample.int(n, ITER_CAP, replace = TRUE)))

  seen <- list(paste(sort(c(hyps$H1$species1, hyps$H1$species2)), collapse = "|"))
  harvested <- list(); minpss <- numeric(0)
  for (it in seq_len(nrow(idx))) {
    ci <- as.integer(unlist(idx[it, , drop = FALSE], use.names = FALSE))
    rows <- do.call(rbind, lapply(seq_len(K), function(k) pools[[k]][ci[k], ]))
    spv <- c(rows$species1, rows$species2)
    if (length(unique(spv)) < 2L * K) next
    sig <- paste(sort(spv), collapse = "|")
    if (sig %in% seen) next
    seen <- c(seen, sig)
    mem <- lapply(seq_len(K), function(i) c(rows$species1[i], rows$species2[i]))
    if (overall_dunn_lean(Dm, mem) >= 1.0) {
      harvested[[length(harvested) + 1L]] <- add_cluster(
        rows[, c("species1", "species2", "pss_score"), drop = FALSE])
      minpss <- c(minpss, suppressWarnings(min(rows$pss_score, na.rm = TRUE)))
    }
  }
  if (length(harvested) > 0L) {
    ord  <- order(-replace(minpss, is.na(minpss), -Inf))
    keep <- head(ord, max(0L, as.integer(max_fop) - 1L))
    harvested <- harvested[keep]
    for (m in seq_along(harvested)) hyps[[paste0("H", m + 1L)]] <- harvested[[m]]
  }
  list(hypotheses = hyps, K = K)
}

#' Lean contrast selection + tiered Dunn validation for one permulated vector.
#'
#' Candidate gate + ranking are identical to the observed selector
#' (selection_algorithm.R::pair_sel.f) via the shared rank_candidates() /
#' greedy_dunn_select() core. The only deliberate difference: the null runs to
#' exactly `target_pairs` and grades independence into tiers, rather than
#' stopping when overall Dunn drops below 1.
#'
#' @param trait_vec       Named numeric vector of permulated trait values.
#' @param D               Patristic distance matrix on the REAL tree.
#' @param target_pairs    N_pairs_obs, the observed independent pair count.
#' @param tree            phylo on the analysis species (for calculate_pairwise_scores).
#' @param cov_bm,cov_ou   BM / OU covariance matrices from the OBSERVED fit
#'                        (pss_core.R::covariances_from_fits), dimnamed by species.
#' @param selected_model  "BM" | "OU" — the observed AIC-selected model.
#' @param ci_lb,ci_ub     Optional per-tip Jeffreys bounds (count data → CI gate).
#' @param top_pct         Top PSS fraction kept as the gate (params.pss_top_pct).
#' @param ordinal         TRUE/FALSE to force the ordinal level gate; NULL = auto.
#' @param n_vec           Optional named per-tip sample sizes → pair_n tiebreak.
#' @return list(tier, n_pairs, dunn_min, n_below, fg, bg, reason)
evaluate_lean_contrast_selection <- function(trait_vec,
                                             D,
                                             target_pairs,
                                             tree = NULL,
                                             cov_bm = NULL,
                                             cov_ou = NULL,
                                             selected_model = "OU",
                                             ci_lb = NULL,
                                             ci_ub = NULL,
                                             top_pct = 0.01,
                                             ordinal = NULL,
                                             n_vec = NULL) {

  reject <- function(reason, n_pairs = 0L, dunn = 0, n_below = NA_integer_) {
    list(tier = 0L, n_pairs = n_pairs, dunn_min = dunn, n_below = n_below,
         fg = NULL, bg = NULL, reason = reason)
  }

  if (target_pairs <= 0L) return(reject("target_pairs <= 0"))
  D <- as.matrix(D)

  sp <- intersect(names(trait_vec), rownames(D))
  if (length(sp) < 2L * target_pairs) return(reject("too few species with distances"))
  trait_vec <- trait_vec[sp]

  use_ci <- !is.null(ci_lb) && !is.null(ci_ub)
  if (is.null(ordinal)) ordinal <- !use_ci && .is_ordinal_vec(trait_vec)

  if (use_ci) {
    lb <- ci_lb[sp]; ub <- ci_ub[sp]
    ok <- is.finite(lb) & is.finite(ub)
    if (sum(ok) < 2L * target_pairs) return(reject("too few species with usable CIs"))
    sp <- sp[ok]; trait_vec <- trait_vec[sp]; lb <- lb[sp]; ub <- ub[sp]
  }

  # --- PSS: vendored phyloq calculate_pairwise_scores (verbatim), on the fixed
  #     observed-model covariances. Returns unordered pairs, FinalScore desc.
  if (is.null(tree) || is.null(cov_bm) || is.null(cov_ou)) {
    return(reject("PSS inputs missing (tree / cov_bm / cov_ou)"))
  }
  tr <- if (length(sp) < length(tree$tip.label)) ape::drop.tip(tree, setdiff(tree$tip.label, sp)) else tree
  tr_sp <- tr$tip.label
  cb <- cov_bm[tr_sp, tr_sp, drop = FALSE]
  co <- cov_ou[tr_sp, tr_sp, drop = FALSE]
  sc <- calculate_pairwise_scores(trait_vec[tr_sp], tr, cb, co, selected_model)

  # Orient each unordered pair: hi = higher trait value (foreground).
  hi_is_1 <- sc$TraitValue1 >= sc$TraitValue2
  c_hi   <- ifelse(hi_is_1, sc$Species1, sc$Species2)
  c_lo   <- ifelse(hi_is_1, sc$Species2, sc$Species1)
  c_dif  <- abs(sc$TraitValue1 - sc$TraitValue2)
  c_dist <- sc$PatristicDistance
  pss    <- sc$FinalScore

  drop0 <- c_dif > 0
  if (!any(drop0)) return(reject("no candidate pairs with positive trait difference"))
  c_hi <- c_hi[drop0]; c_lo <- c_lo[drop0]; c_dif <- c_dif[drop0]
  c_dist <- c_dist[drop0]; pss <- pss[drop0]

  # --- Stage 1: trait-type non-overlap ---
  if (use_ci) {
    type_mask <- (lb[c_hi] > ub[c_lo])
  } else if (isTRUE(ordinal)) {
    lv_hi <- max(trait_vec, na.rm = TRUE); lv_lo <- min(trait_vec, na.rm = TRUE)
    type_mask <- (trait_vec[c_hi] >= lv_hi) & (trait_vec[c_lo] <= lv_lo)
  } else {
    type_mask <- rep(TRUE, length(c_hi))
  }
  if (!any(type_mask)) return(reject("no pair passes the trait-type gate"))
  s1 <- which(type_mask)

  # --- Stage 2 (continuous only): top `top_pct` of stage-1 survivors by PSS.
  #     Count/ordinal keep every stage-1 pair; PSS just ranks them. ---
  if (!use_ci && !isTRUE(ordinal)) {
    n_keep <- min(max(1L, ceiling(length(s1) * top_pct)), length(s1))
    keep_i <- s1[order(pss[s1], decreasing = TRUE)][seq_len(n_keep)]
  } else {
    keep_i <- s1
  }
  keep <- logical(length(c_hi)); keep[keep_i] <- TRUE

  cand_df <- data.frame(
    species1 = c_hi[keep], species2 = c_lo[keep],
    distance = c_dist[keep], abs_diff = c_dif[keep],
    pss_score = pss[keep], stringsAsFactors = FALSE
  )
  if (!is.null(n_vec)) {
    cand_df$pair_n <- as.numeric(n_vec[cand_df$species1]) + as.numeric(n_vec[cand_df$species2])
  }

  ranked <- rank_candidates(cand_df)
  res <- greedy_dunn_select(ranked, D, target = target_pairs, enforce_dunn = FALSE)

  members <- res$members
  n_pairs <- length(members)
  if (n_pairs < target_pairs) {
    return(reject("could not form target_pairs non-overlapping pairs", n_pairs))
  }

  dunn_vec <- vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1))
  n_below  <- sum(dunn_vec < 1)
  dunn_min <- min(dunn_vec)
  if (n_below >= 2L) {
    return(reject("two or more pairs below Dunn 1", n_pairs, dunn_min, n_below))
  }

  sel_fg <- res$selected$species1
  sel_bg <- res$selected$species2
  list(tier = if (n_below == 0L) 1L else 2L,
       n_pairs = n_pairs, dunn_min = dunn_min, n_below = n_below,
       fg = sel_fg, bg = sel_bg,
       fg_values = unname(trait_vec[sel_fg]), bg_values = unname(trait_vec[sel_bg]),
       mean_pd = mean(vapply(members, function(x) D[x[1], x[2]], numeric(1))),
       mean_df = mean(vapply(members, function(x) abs(trait_vec[x[1]] - trait_vec[x[2]]), numeric(1))),
       mode = if (use_ci) "ci" else if (isTRUE(ordinal)) "ordinal" else "pss",
       reason = "accepted")
}
