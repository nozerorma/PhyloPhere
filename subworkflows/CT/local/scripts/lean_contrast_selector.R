#!/usr/bin/env Rscript
# lean_contrast_selector.R — Contrast selection shared by the observed selector and the permulation null.
# PhyloPhere | subworkflows/CT/local/scripts/
# Sourced by: permulations.R, selection_algorithm.R (fop_pair_sel.f, observed selector of TRAIT_ANALYSIS)
# =============================================================================
# Candidates pass the trait-type gate (count: Jeffreys CI non-overlap; ordinal: top
# level vs bottom level; continuous: top `pss_top_pct` by Phylogenetic Shift Score under
# the AIC-selected OU/BM model) and are ranked by PSS. The canonical contrast (H1) is
# assembled greedily under the modified Dunn index; the FOP harvest adds Dunn-independent
# alternative hypotheses (H2..Hn) drawn from the Voronoi domains of the canonical pairs.
# The file defines functions only. The PSS engine is pss_core.R.
# =============================================================================

# ── Dunn index ────────────────────────────────────────────────────────────────

# Modified Dunn index for one cluster: (min distance to any other cluster)
# divided by (this cluster's own diameter).
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

# Integer-indexed specialization of overall_dunn_lean for the hot loop of the FOP
# harvest: `pi1`/`pi2` are length-K integer row/col indices into the bare matrix `Dm`
# (cluster k is the pair pi1[k]-pi2[k]). It returns the same min-over-k modified Dunn as
# overall_dunn_lean(Dm, <same members as names>), without character subscripting,
# per-call closure allocation or the `members` list. A cluster with zero diameter
# contributes Inf (skipped), as in mod_dunn_lean.
overall_dunn_int <- function(Dm, pi1, pi2, K) {
  best <- Inf
  for (k in seq_len(K)) {
    a <- pi1[k]; b <- pi2[k]
    intra <- Dm[a, b]
    if (intra == 0) next
    inter <- Inf
    for (j in seq_len(K)) {
      if (j == k) next
      cc <- pi1[j]; dd <- pi2[j]
      m <- min(Dm[a, cc], Dm[a, dd], Dm[b, cc], Dm[b, dd])
      if (m < inter) inter <- m
    }
    r <- inter / intra
    if (r < best) best <- r
  }
  best
}

# ── PSS engine ────────────────────────────────────────────────────────────────

# PSS scoring comes from the vendored phyloq engine (pss_core.R): analytical_s() and
# calculate_pairwise_scores(). If they are not loaded yet, pss_core.R is looked for in
# src/ and in the working directory, then at a fixed absolute path of the development
# checkout.
if (!exists("calculate_pairwise_scores", mode = "function")) {
  .pss_core_paths <- c(
    file.path(getwd(), "src", "pss_core.R"),
    file.path(getwd(), "pss_core.R"),
    "/home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/CT/local/scripts/pss_core.R"
  )
  .hit <- .pss_core_paths[file.exists(.pss_core_paths)]
  if (length(.hit)) source(.hit[1])
}

# ── Shared selection core ─────────────────────────────────────────────────────

# Ordinal code auto-detection: two to five integer levels (mirrors stats.R::is_ordinal_trait "auto").
.is_ordinal_vec <- function(v) {
  u <- unique(v[!is.na(v)])
  length(u) >= 2 && length(u) <= 5 && all(u == round(u))
}

# The observed selector (selection_algorithm.R::fop_pair_sel.f) and the permulation null
# (permulations.R) select contrasts through these functions only: selection_context() ->
# lean_candidate_df() -> greedy_dunn_select() -> lean_fop_harvest(). Run on the real
# labeling, the null therefore reproduces the observed selection.

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
#' combined pair sample size desc when available, then species names.
#'
#' Numeric keys are rounded (PSS to 4 decimals, distance / |difference| to 6)
#' so floating-point noise never decides a tie, and full ties fall back to the
#' species names (C-locale radix order), never to the incoming row order. The
#' ranking is therefore a function of the candidate set alone. Binary and
#' ordinal traits make PSS ties the rule rather than the exception.
rank_candidates <- function(df) {
  if (nrow(df) == 0) return(df)
  has_pss <- "pss_score" %in% names(df) && any(is.finite(df$pss_score))
  keys <- if (has_pss) list(-round(df$pss_score, 4)) else list(round(df$distance, 6))
  keys <- c(keys, list(-round(df$abs_diff, 6)))
  if ("pair_n" %in% names(df)) keys <- c(keys, list(-df$pair_n))
  keys <- c(keys, list(as.character(df$species1), as.character(df$species2)))
  df[do.call(order, c(keys, list(method = "radix"))), , drop = FALSE]
}

#' Tree, distances and evolutionary model for contrast selection on one trait.
#'
#' The single place both the observed selector and the permulation null derive
#' their selection inputs from: tree pruned to the trait species, dichotomised
#' (zero-length edges nudged to 1e-8), patristic distances, and the AIC-selected
#' BM/OU fit with its covariances (vendored phyloq, pss_core.R).
#'
#' @param trait_vec   named numeric trait values (NA dropped).
#' @param tree        phylo covering the trait species.
#' @param force_model NULL (AIC) | "BM" | "OU".
#' @return list(tree, trait_vec, D, fits, selected_model, cov_bm, cov_ou)
selection_context <- function(trait_vec, tree, force_model = NULL) {
  trait_vec <- trait_vec[is.finite(trait_vec)]
  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, names(trait_vec)))
  tr <- ape::multi2di(tr, random = FALSE)
  tr$edge.length[tr$edge.length <= 0] <- 1e-8
  trait_vec <- trait_vec[tr$tip.label]
  fits  <- fit_models(tr, trait_vec)
  model <- select_model(fits, force_model = force_model)
  cv    <- covariances_from_fits(tr, fits)
  list(tree = tr, trait_vec = trait_vec, D = ape::cophenetic.phylo(tr),
       fits = fits, selected_model = model, cov_bm = cv$BM, cov_ou = cv$OU)
}

#' Evaluate `expr` under set.seed(seed) and restore the caller's RNG state, so
#' a seeded draw inside the selector never resets the stream of the code that
#' calls it (e.g. the permulation simulations in permulations.R).
.with_seed <- function(seed, expr) {
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", old, envir = genv)
    }
  })
  set.seed(seed)
  expr
}

#' Voronoi domain of every species: index of the nearest canonical pair
#' (distance to its closer member). Distances are rounded to 6 decimals so
#' floating-point noise never decides a tie; exact ties go to the lowest pair.
voronoi_domains <- function(Dm, members) {
  all_sp <- rownames(Dm)
  dom <- vapply(all_sp, function(s)
    which.min(round(vapply(seq_along(members), function(k)
      min(Dm[s, members[[k]][1]], Dm[s, members[[k]][2]]), numeric(1)), 6)),
    integer(1))
  names(dom) <- all_sp
  dom
}

#' Unified greedy Dunn-gated pair assembly with maximum phylogenetic tree dispersion.
#'
#' @param ranked candidate pairs already ordered best-first by rank_candidates().
#'   Needs species1, species2, distance, abs_diff (+ pair_n optional).
#' @param D patristic distance matrix (species x species).
#' @param target stop after this many pairs (Inf = run until enforce_dunn stops it).
#' @param enforce_dunn TRUE  -> only accept a pair that keeps every cluster's
#'                              modified Dunn >= 1; stop when none qualifies
#'                              (observed selector: variable pair count).
#'                     FALSE -> take the best candidate up to `target`, even below 1;
#'                              caller grades the result (permulation null: fixed N, tiered).
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

    # Minimum patristic distance from this pair to any already-selected pair
    inter_dists <- vapply(seq_len(nrow(cand)), function(i) {
      s1 <- cand$species1[i]; s2 <- cand$species2[i]
      min(vapply(members, function(m) {
        min(D[s1, m[1]], D[s1, m[2]], D[s2, m[1]], D[s2, m[2]])
      }, numeric(1)))
    }, numeric(1))
    cand$inter_dist <- round(inter_dists, 4)

    if (enforce_dunn) {
      cand <- cand[cand$Dunn_index >= 1, , drop = FALSE]
      if (nrow(cand) == 0) break

      # Order by: 1) maximum inter-cluster distance (broadest tree dispersion / coverage),
      # 2) highest Dunn index, 3) stable incoming rank order (PSS / rank_candidates)
      cand <- cand[order(-cand$inter_dist, -cand$Dunn_index, seq_len(nrow(cand))), , drop = FALSE]

      # Find first candidate whose addition keeps overall Dunn >= 1 across all clusters
      picked <- FALSE
      for (i in seq_len(nrow(cand))) {
        best_cand <- cand[i, , drop = FALSE]
        new_members <- c(members, list(c(best_cand$species1, best_cand$species2)))
        if (overall_dunn_lean(D, new_members) >= 1) {
          best_cand$cluster <- length(members) + 1L
          selected <- rbind(selected, best_cand[, setdiff(names(best_cand), "inter_dist"), drop = FALSE])
          members  <- new_members
          used     <- c(used, best_cand$species1, best_cand$species2)
          picked   <- TRUE
          break
        }
      }
      if (!picked) break
    } else {
      # Without enforce_dunn (tiered null fallback): order by inter-cluster distance then Dunn
      cand <- cand[order(-cand$inter_dist, -cand$Dunn_index, seq_len(nrow(cand))), , drop = FALSE]
      best <- cand[1, , drop = FALSE]
      best$cluster <- length(members) + 1L
      selected <- rbind(selected, best[, setdiff(names(best), "inter_dist"), drop = FALSE])
      members  <- c(members, list(c(best$species1, best$species2)))
      used     <- c(used, best$species1, best$species2)
    }
  }

  list(selected = selected, members = members)
}

#' Candidate-pair gate and ranking for one (permulated) trait vector.
#'
#' Shared by evaluate_lean_contrast_selection() (the tiered Dunn null) and
#' lean_fop_harvest() (the per-cycle FOP mirror). Reproduces stages 1-2 of the
#' observed selector: PSS via the vendored phyloq engine on the fixed observed-
#' model covariances, then the trait-type gate (CI non-overlap / ordinal levels /
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

#' Canonical pairs of a permulated trait whose PSS profile matches the observed one.
#'
#' The observed run assembles its pairs until the Dunn index stops it, so its last pairs are the
#' poorest the trait offers; a null that keeps the K best pairs of a richer pool is systematically
#' closer. Here the null builds its K pairs one by one: for observed pair i (in the order the
#' observed selector chose them) it takes the candidate whose PSS is closest to the observed PSS_i
#' (in log scale) among those that keep every pair Dunn-independent (modified Dunn >= 1).
#' The draw is comparable when the worst pair is within `tol` of its target, i.e.
#' max_i |log(PSS_i / target_i)| <= log(1 + tol).
#'
#' @param ranked      candidate pairs (lean_candidate_df()$cand_df), with species1, species2, pss_score.
#' @param D           patristic distance matrix.
#' @param target_pss  PSS of the observed canonical pairs, in selection order.
#' @param tol         relative tolerance of the PSS of each pair.
#' @param max_probe   nearest candidates tested against the Dunn gate before giving up on a pair.
#' @return list(selected = rows of `ranked` in selection order | NULL, members, mismatch = max |log ratio|,
#'              reason)
match_pss_select <- function(ranked, D, target_pss, tol = 0.25, max_probe = 60L) {
  D <- as.matrix(D)
  fail <- function(reason) list(selected = NULL, members = list(), mismatch = NA_real_, reason = reason)
  if (is.null(ranked) || nrow(ranked) == 0L) return(fail("no candidate pairs"))
  lp <- log(pmax(ranked$pss_score, 1e-12)); lt <- log(pmax(target_pss, 1e-12))
  members <- list(); used <- character(0); picked <- integer(0)
  for (i in seq_along(lt)) {
    avail <- which(!(ranked$species1 %in% used | ranked$species2 %in% used))
    if (!length(avail)) return(fail("could not form target_pairs non-overlapping pairs"))
    ord <- avail[order(abs(lp[avail] - lt[i]))]
    pick <- NA_integer_
    for (j in utils::head(ord, max_probe)) {
      m2 <- c(members, list(c(ranked$species1[j], ranked$species2[j])))
      if (length(members) == 0L || (mod_dunn_lean(D, m2, length(m2)) >= 1 && overall_dunn_lean(D, m2) >= 1)) { pick <- j; break }
    }
    if (is.na(pick)) return(fail("no Dunn-independent candidate for a pair of the PSS profile"))
    members <- c(members, list(c(ranked$species1[pick], ranked$species2[pick])))
    used <- c(used, ranked$species1[pick], ranked$species2[pick]); picked <- c(picked, pick)
  }
  sel <- ranked[picked, , drop = FALSE]
  mismatch <- max(abs(log(pmax(sel$pss_score, 1e-12)) - lt))
  if (mismatch > log1p(tol)) return(list(selected = NULL, members = members, mismatch = mismatch,
                                         reason = "PSS profile not matched within the tolerance"))
  list(selected = sel, members = members, mismatch = mismatch, reason = NULL)
}

#' Pool and capacity of one trait vector under the observed selection rule.
#'
#' The observed selector assembles pairs until the Dunn index stops it, so its K is the capacity of the trait:
#' the most mutually independent pairs its candidate pool can give. This returns the size of the pool with the
#' PSS top_pct gate and without it (the pairs the gate discards), and that capacity on the gated pool.
#' Costly (an exhaustive greedy assembly): meant for the observed trait and a sample of the null draws.
#' @return c(pool_gated, pool_ungated, capacity); NA where the trait has no candidate pair.
lean_draw_capacity <- function(trait_vec, D, tree, cov_bm, cov_ou, selected_model,
                               ci_lb = NULL, ci_ub = NULL, top_pct = 0.01, ordinal = NULL, n_vec = NULL) {
  D <- as.matrix(D)
  gated   <- lean_candidate_df(trait_vec, D, 1L, tree, cov_bm, cov_ou, selected_model, ci_lb, ci_ub, top_pct, ordinal, n_vec)$cand_df
  ungated <- lean_candidate_df(trait_vec, D, 1L, tree, cov_bm, cov_ou, selected_model, ci_lb, ci_ub, 1,       ordinal, n_vec)$cand_df
  if (is.null(gated) || is.null(ungated)) return(c(pool_gated = NA_real_, pool_ungated = NA_real_, capacity = NA_real_))
  cap <- nrow(greedy_dunn_select(gated, D, target = Inf, enforce_dunn = TRUE)$selected)
  c(pool_gated = nrow(gated), pool_ungated = nrow(ungated), capacity = cap)
}

#' FOP multi-hypothesis harvest, shared by the observed selector
#' (selection_algorithm.R::fop_pair_sel.f) and the permulation null.
#'
#' Candidate gate/rank from lean_candidate_df. Canonical H1 = `canon_pairs`
#' when supplied, else greedy_dunn_select(enforce_dunn = TRUE);
#' Voronoi-partition the species by nearest canonical pair (voronoi_domains);
#' draw one in-domain candidate per domain (exhaustive <= ITER_CAP, else capped
#' random draws seeded with `seed`); keep Dunn >= 1; rank by min PSS, mean PSS,
#' overall Dunn (PSS/Dunn rounded to 4 decimals, then the species signature) and
#' keep the top (max_fop - 1). Pools are taken in rank_candidates order, so the
#' draws depend only on the candidate set and the seed.
#'
#' @param seed        RNG seed for the capped random draws (the pipeline seed).
#' @param canon_pairs optional data.frame(species1, species2) — the canonical
#'   H1 pairs (observed: greedy_dunn_select(enforce_dunn = TRUE); null: the
#'   tiered fixed-N selection of evaluate_lean_contrast_selection). H1 is taken
#'   verbatim and seeds the Voronoi partition.
#' @return list(hypotheses = named list H1..Hn of data.frame(species1, species2,
#'              distance, abs_diff, pss_score, cluster), dunn = named numeric
#'              overall Dunn per hypothesis, species_domain = named integer, K)
#'              or list(hypotheses = list(), K = 0).
lean_fop_harvest <- function(trait_vec, D, target_pairs,
                             tree = NULL, cov_bm = NULL, cov_ou = NULL,
                             selected_model = "OU", ci_lb = NULL, ci_ub = NULL,
                             top_pct = 0.01, ordinal = NULL, n_vec = NULL,
                             max_fop = 100L, seed, canon_pairs = NULL) {
  cc <- lean_candidate_df(trait_vec, D, target_pairs, tree, cov_bm, cov_ou,
                          selected_model, ci_lb, ci_ub, top_pct, ordinal, n_vec)
  cand_df <- cc$cand_df
  if (is.null(cand_df) || nrow(cand_df) < 1) return(list(hypotheses = list(), K = 0L))

  Dm <- as.matrix(D)
  add_cluster <- function(df) { df$cluster <- seq_len(nrow(df)); df }

  # Candidate-pool row (either orientation) for a canonical pair.
  .cand_of <- function(s1, s2) {
    hit <- which((cand_df$species1 == s1 & cand_df$species2 == s2) |
                 (cand_df$species1 == s2 & cand_df$species2 == s1))
    if (length(hit)) hit[1] else NA_integer_
  }

  if (!is.null(canon_pairs) && nrow(canon_pairs) >= 1) {
    hit <- mapply(.cand_of, canon_pairs$species1, canon_pairs$species2)
    cp <- data.frame(
      species1  = as.character(canon_pairs$species1),
      species2  = as.character(canon_pairs$species2),
      distance  = if ("distance" %in% names(canon_pairs)) canon_pairs$distance
                  else cand_df$distance[hit],
      abs_diff  = if ("abs_diff" %in% names(canon_pairs)) canon_pairs$abs_diff
                  else cand_df$abs_diff[hit],
      pss_score = if ("pss_score" %in% names(canon_pairs)) canon_pairs$pss_score else cand_df$pss_score[hit],
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

  h_cols <- c("species1", "species2", "distance", "abs_diff", "pss_score")
  hyps <- list(H1 = add_cluster(canon$selected[, h_cols, drop = FALSE]))
  dunn_out <- c(H1 = overall_dunn_lean(Dm, cm))
  dom <- voronoi_domains(Dm, cm)
  if (K < target_pairs) {  # H1 only, degenerate
    return(list(hypotheses = hyps, dunn = dunn_out, species_domain = dom, K = K))
  }

  pools <- lapply(seq_len(K), function(k) {
    ds <- names(dom)[dom == k]
    cand_df[cand_df$species1 %in% ds & cand_df$species2 %in% ds, , drop = FALSE]
  })
  psz <- vapply(pools, nrow, integer(1))
  if (any(psz == 0L)) {
    return(list(hypotheses = hyps, dunn = dunn_out, species_domain = dom, K = K))
  }

  ITER_CAP <- as.integer(max_fop) * 20L
  total <- prod(as.numeric(psz))
  idx <- if (total <= ITER_CAP) {
    do.call(expand.grid, c(lapply(psz, seq_len), list(KEEP.OUT.ATTRS = FALSE)))
  } else {
    .with_seed(seed, as.data.frame(lapply(psz, function(n) sample.int(n, ITER_CAP, replace = TRUE))))
  }

  # The hot loop below runs up to ITER_CAP (= max_fop * 20) draws per call, and the null
  # calls it once per cycle. The draw index is materialized once as an integer matrix and
  # each Voronoi pool as bare vectors; deduplication goes through a hashed environment and
  # the Dunn test runs on integer indices (overall_dunn_int).
  idx_mat  <- matrix(as.integer(unlist(idx, use.names = FALSE)), ncol = length(psz))
  pool_s1  <- lapply(pools, `[[`, "species1")
  pool_s2  <- lapply(pools, `[[`, "species2")
  pool_pss <- lapply(pools, `[[`, "pss_score")
  pool_dst <- lapply(pools, `[[`, "distance")
  pool_dif <- lapply(pools, `[[`, "abs_diff")
  Kseq <- seq_len(K)

  # Species -> integer row/col index into Dm, once. The hot loop below then runs the
  # Dunn test and the dedup signature on integers: name-indexing of D (as in
  # mod_dunn_lean) is the costly operation there.
  sp_idx  <- setNames(seq_len(nrow(Dm)), rownames(Dm))
  pool_i1 <- lapply(pool_s1, function(s) sp_idx[s])
  pool_i2 <- lapply(pool_s2, function(s) sp_idx[s])

  seen <- new.env(parent = emptyenv())
  h1i <- sp_idx[c(hyps$H1$species1, hyps$H1$species2)]
  assign(paste(sort(h1i), collapse = "|"), TRUE, envir = seen)
  n_it <- nrow(idx_mat)
  # The per-hypothesis data.frame is built only at the final max_fop cut below, because
  # most iterations are duplicates or fail the Dunn gate; until then the plain
  # species/PSS vectors are stored.
  harv_s1 <- vector("list", n_it); harv_s2 <- vector("list", n_it)
  harv_ps <- vector("list", n_it); harv_dst <- vector("list", n_it)
  harv_dif <- vector("list", n_it); harv_sig <- character(n_it)
  minpss <- numeric(n_it); meanpss <- numeric(n_it); hdunn <- numeric(n_it); nh <- 0L
  s1 <- character(K); s2 <- character(K); ps <- numeric(K)
  dst <- numeric(K); dif <- numeric(K)
  i1 <- integer(K);   i2 <- integer(K)
  for (it in seq_len(n_it)) {
    ci <- idx_mat[it, ]
    for (k in Kseq) {
      r <- ci[k]
      s1[k] <- pool_s1[[k]][r]; s2[k] <- pool_s2[[k]][r]; ps[k] <- pool_pss[[k]][r]
      dst[k] <- pool_dst[[k]][r]; dif[k] <- pool_dif[[k]][r]
      i1[k] <- pool_i1[[k]][r]; i2[k] <- pool_i2[[k]][r]
    }
    spvi <- c(i1, i2)
    if (anyDuplicated(spvi) > 0L) next
    sig <- paste(sort.int(spvi), collapse = "|")
    if (!is.null(seen[[sig]])) next
    assign(sig, TRUE, envir = seen)
    dn <- overall_dunn_int(Dm, i1, i2, K)
    if (dn >= 1.0) {
      nh <- nh + 1L
      harv_s1[[nh]] <- s1; harv_s2[[nh]] <- s2; harv_ps[[nh]] <- ps
      harv_dst[[nh]] <- dst; harv_dif[[nh]] <- dif; harv_sig[nh] <- sig
      minpss[nh]  <- suppressWarnings(min(ps, na.rm = TRUE))
      meanpss[nh] <- suppressWarnings(mean(ps, na.rm = TRUE))
      hdunn[nh]   <- dn
    }
  }
  if (nh > 0L) {
    sq <- seq_len(nh)
    .key <- function(x) -round(replace(x[sq], !is.finite(x[sq]), -Inf), 4)
    ord  <- order(.key(minpss), .key(meanpss), .key(hdunn), harv_sig[sq], method = "radix")
    keep <- head(ord, max(0L, as.integer(max_fop) - 1L))
    for (m in seq_along(keep)) {
      j <- keep[m]
      h_id <- paste0("H", m + 1L)
      hyps[[h_id]] <- add_cluster(data.frame(
        species1 = harv_s1[[j]], species2 = harv_s2[[j]],
        distance = harv_dst[[j]], abs_diff = harv_dif[[j]], pss_score = harv_ps[[j]],
        stringsAsFactors = FALSE))
      dunn_out[h_id] <- hdunn[j]
    }
  }
  list(hypotheses = hyps, dunn = dunn_out, species_domain = dom, K = K)
}

#' Lean contrast selection + tiered Dunn validation for one permulated vector.
#'
#' Candidate gate + ranking are the observed selector's (lean_candidate_df +
#' greedy_dunn_select). Without `pss_profile`, the null runs to exactly
#' `target_pairs` and grades independence into tiers, rather than stopping when
#' overall Dunn drops below 1. With `pss_profile`, the canonical pairs follow the
#' PSS of the observed pairs (match_pss_select) and the draw is Tier 1 or rejected.
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
#' @param pss_profile     NULL, or the PSS of the observed canonical pairs in selection order. The canonical pairs
#'                        are then chosen by match_pss_select() from the pool without the PSS top_pct gate
#'                        (the pairs themselves set the PSS range), and the draw is accepted (Tier 1) only if
#'                        every pair is within `pss_tol` of its target. Ignored for count (CI) traits.
#' @param pss_tol         Relative PSS tolerance of the matching.
#' @return list(tier, n_pairs, dunn_min, n_below, fg, bg, canon, mean_pss, mismatch, reason)
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
                                             n_vec = NULL,
                                             pss_profile = NULL,
                                             pss_tol = 0.25) {

  reject <- function(reason, n_pairs = 0L, dunn = 0, n_below = NA_integer_, mismatch = NA_real_) {
    list(tier = 0L, n_pairs = n_pairs, dunn_min = dunn, n_below = n_below,
         fg = NULL, bg = NULL, mismatch = mismatch, reason = reason)
  }

  if (target_pairs <= 0L) return(reject("target_pairs <= 0"))
  D <- as.matrix(D)

  use_match <- !is.null(pss_profile) && is.null(ci_lb)
  if (use_match && length(pss_profile) != target_pairs) return(reject("pss_profile length differs from target_pairs"))
  cc <- lean_candidate_df(trait_vec, D, target_pairs, tree, cov_bm, cov_ou,
                          selected_model, ci_lb, ci_ub, if (use_match) 1 else top_pct, ordinal, n_vec)
  if (is.null(cc$cand_df)) return(reject(cc$reason))
  ranked <- cc$cand_df

  if (use_match) {
    mt <- match_pss_select(ranked, D, pss_profile, pss_tol)
    if (is.null(mt$selected)) return(reject(mt$reason, length(mt$members), mismatch = mt$mismatch))
    sel <- mt$selected; members <- mt$members
    return(list(tier = 1L, n_pairs = length(members), dunn_min = overall_dunn_lean(D, members), n_below = 0L,
                fg = sel$species1, bg = sel$species2,
                canon = sel[, c("species1", "species2", "distance", "abs_diff", "pss_score")],
                fg_values = unname(trait_vec[sel$species1]), bg_values = unname(trait_vec[sel$species2]),
                mean_pd = mean(sel$distance), mean_df = mean(sel$abs_diff), mean_pss = mean(sel$pss_score),
                mismatch = mt$mismatch, mode = cc$mode, reason = "accepted"))
  }

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
       canon = res$selected[, c("species1", "species2", "distance", "abs_diff", "pss_score")],
       fg_values = unname(trait_vec[sel_fg]), bg_values = unname(trait_vec[sel_bg]),
       mean_pd = mean(vapply(members, function(x) D[x[1], x[2]], numeric(1))),
       mean_df = mean(vapply(members, function(x) abs(trait_vec[x[1]] - trait_vec[x[2]]), numeric(1))),
       mean_pss = mean(res$selected$pss_score), mismatch = NA_real_,
       mode = cc$mode,
       reason = "accepted")
}
