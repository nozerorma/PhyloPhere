#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: FOP multi-hypothesis -> domain-pooled position score
# File: subworkflows/SCORING/local/src/fop_pool.R
# =============================================================================
# Collapses the FOP alternative-hypothesis harvest (one disambiguation row per
# (Gene, Position, scheme, hypothesis)) into ONE ASR path score per
# (Gene, Position, scheme), the way scoring_compute.R section 2g needs it.
#
# Why not a plain mean over hypotheses: H1..Hn are NOT independent replicates.
# They are overlapping K-pair designs drawn from the same Voronoi domains, so a
# uniform mean both dilutes a strong canonical signal and lets a position with
# many harvested hypotheses distort every genome-wide rank it enters.
#
# Domain-pooled scoring (see memory project_fop_multihyp_scoring):
#   * pair index i in `mrca_<i>_path_score` == Voronoi domain i, by construction
#     of the harvest (fop_pair_sel.f numbers each hypothesis's pairs in domain
#     order; disambiguation enumerates traitfile_H{m}.tab in that same order).
#   * within domain i, pool s(p, site) over the DISTINCT candidate pairs that
#     domain contributed across the harvest -> c_i(site). Pooling fn:
#     PSS-weighted mean (weights from contrast_hypotheses_pairs.tsv; equal
#     weights when PSS is missing). A pair the harvest kept but that scores a
#     weak isolation is genuine evidence the clade's signal is fragile, so it
#     is not discarded, only down-weighted.
#   * recombine c_1..c_K into the position score with the SAME algebra as
#     path_scores.py::compute_asr_path_score:
#       core        = P(>=2 domains carry a clean change)   [inclusion-exclusion]
#       replication = independence * core
#       strength    = (0.75 + 0.25 * diversity) * derived_agreement
#       asr         = replication * strength * conservation_gate
#     independence / diversity / derived_agreement / conservation_gate are NOT
#     re-derivable from the flat TSV (they need the LAC tree structure), so the
#     per-hypothesis values are pooled by the same PSS weights. Voronoi domains
#     are disjoint, so cross-domain independence is ~1 by construction and the
#     pooled independence only ever mildly discounts.
#   * `core` here treats each c_i as P(genuine change) without the per-side
#     (top/bottom) split path_scores.py uses, because the flat TSV already
#     averaged the two sides into `mrca_<i>_path_score`. Opposite-direction
#     pairs in different domains can therefore be over-credited slightly; the
#     per-hypothesis `core` column (pooled, below) is carried alongside as the
#     conservative cross-check.
#
# Single-hypothesis / non-FOP input is a no-op: the one row passes through with
# its asr_path_score untouched (verified by test_fop_pool.R).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Exact P(>= 2 successes) over independent Bernoullis — verbatim algebra of
# path_scores.py::_p_at_least_2 (inclusion-exclusion on P0 and P1).
.p_at_least_2 <- function(p) {
  p <- p[is.finite(p)]
  if (length(p) < 2) return(0)
  p0 <- prod(1 - p)
  p1 <- sum(vapply(seq_along(p), function(i) p[i] * prod(1 - p[-i]), numeric(1)))
  max(0, min(1, 1 - p0 - p1))
}

# Weighted mean that falls back to the plain mean when no usable weights.
.wmean <- function(x, w) {
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  if (!any(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) return(mean(x))
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) return(mean(x))
  sum(x * w) / sum(w)
}

#' Parse contrast_hypotheses_pairs.tsv into a (hyp_id, domain) -> pss lookup.
#'
#' Returns a data.frame with columns hyp_id, pair (domain id), pss_score.
#' NULL / missing file -> NULL (caller then pools with equal weights).
read_hypothesis_pairs <- function(path) {
  if (is.null(path) || !nzchar(path) || grepl("^NO_", basename(path)) || !file.exists(path)) {
    return(NULL)
  }
  hp <- tryCatch(
    read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(hp) || nrow(hp) == 0) return(NULL)
  needed <- c("hypothesis_id", "pair", "pss_score")
  if (!all(needed %in% names(hp))) return(NULL)
  data.frame(
    hyp_id    = sub(".*?(H[0-9]+).*", "\\1", as.character(hp$hypothesis_id)),
    pair      = suppressWarnings(as.integer(hp$pair)),
    pss_score = suppressWarnings(as.numeric(hp$pss_score)),
    stringsAsFactors = FALSE
  )
}

#' Domain-pool the FOP hypothesis rows of one scheme-level data.frame.
#'
#' @param df data.frame of disambiguation rows for a SINGLE (Gene, Position,
#'   caap_group) group — one row per hypothesis. Must carry: asr_path_score,
#'   independence, mrca_diversity, derived_agreement, conservation_gate, core,
#'   hyp_id, and the mrca_<i>_path_score / mrca_<i>_node columns.
#' @param path_cols character vector of the mrca_<i>_path_score column names.
#' @param node_cols character vector of the mrca_<i>_node column names (same order).
#' @param hyp_pairs (hyp_id, pair, pss_score) lookup from read_hypothesis_pairs(),
#'   or NULL for equal-weight pooling.
#' @param diversity_floor lower bound of the diversity multiplier (path_scores.py
#'   default 0.75).
#' @return one-row data.frame: asr_path_score, independence, mrca_diversity,
#'   derived_agreement, conservation_gate, core, core_perside_pooled,
#'   n_hypotheses, plus c_1..c_K as pooled_domain_<i>.
pool_group <- function(df, path_cols, node_cols, hyp_pairs = NULL,
                       diversity_floor = 0.75) {
  hyps <- unique(df$hyp_id[!is.na(df$hyp_id) & nzchar(df$hyp_id)])

  # No FOP structure (single contrast, or a lone hypothesis) -> pass through.
  if (length(hyps) <= 1) {
    take <- df[which.max(ifelse(is.finite(df$asr_path_score), df$asr_path_score, -Inf)), ]
    if (nrow(take) == 0) take <- df[1, ]
    return(data.frame(
      asr_path_score     = take$asr_path_score[1],
      independence       = take$independence[1],
      mrca_diversity     = take$mrca_diversity[1],
      derived_agreement  = take$derived_agreement[1],
      conservation_gate  = take$conservation_gate[1],
      core               = take$core[1],
      core_perside_pooled = take$core[1],
      n_hypotheses       = length(hyps),
      supporting_hypotheses = if (length(hyps)) paste(sort(hyps), collapse = ",") else "",
      stringsAsFactors = FALSE
    ))
  }

  # Per-hypothesis weight = min pair PSS across its K domains (a hypothesis is
  # only as strong as its weakest contrast — mirrors the FOP ranking key).
  hyp_w <- setNames(rep(NA_real_, length(hyps)), hyps)
  if (!is.null(hyp_pairs)) {
    for (h in hyps) {
      w <- hyp_pairs$pss_score[hyp_pairs$hyp_id == h]
      if (length(w) && any(is.finite(w))) hyp_w[h] <- min(w[is.finite(w)])
    }
  }
  row_w <- unname(hyp_w[df$hyp_id])

  # ── Per-domain pool: c_i over the distinct pairs domain i contributed ──────
  K <- length(path_cols)
  c_vec <- numeric(K)
  for (i in seq_len(K)) {
    s <- suppressWarnings(as.numeric(df[[path_cols[i]]]))
    nd <- if (i <= length(node_cols) && node_cols[i] %in% names(df)) {
      as.character(df[[node_cols[i]]])
    } else {
      as.character(df$hyp_id)  # fall back: treat every hypothesis's pair as distinct
    }
    keep <- is.finite(s)
    if (!any(keep)) { c_vec[i] <- NA_real_; next }
    s <- s[keep]; nd <- nd[keep]; wr <- row_w[keep]
    # Distinct candidate pairs in this domain (dedup by MRCA node); weight each
    # by its hypothesis PSS (min-pair), else equal.
    dd <- data.frame(nd = nd, s = s, w = wr, stringsAsFactors = FALSE) %>%
      group_by(nd) %>%
      summarise(s = dplyr::first(s),
                w = suppressWarnings(max(w, na.rm = TRUE)), .groups = "drop")
    c_vec[i] <- .wmean(dd$s, dd$w)
  }

  # ── Recombine with path_scores.py algebra ────────────────────────────────
  core_pooled <- .p_at_least_2(c_vec)
  indep_pooled <- .wmean(suppressWarnings(as.numeric(df$independence)), row_w)
  div_pooled   <- .wmean(suppressWarnings(as.numeric(df$mrca_diversity)), row_w)
  da_pooled    <- .wmean(suppressWarnings(as.numeric(df$derived_agreement)), row_w)
  cg_pooled    <- .wmean(suppressWarnings(as.numeric(df$conservation_gate)), row_w)
  core_row_pooled <- .wmean(suppressWarnings(as.numeric(df$core)), row_w)

  if (!is.finite(indep_pooled)) indep_pooled <- 1
  if (!is.finite(div_pooled))   div_pooled   <- 0
  if (!is.finite(da_pooled))    da_pooled    <- 1
  if (!is.finite(cg_pooled))    cg_pooled    <- 1

  diversity_mult <- diversity_floor + (1 - diversity_floor) * div_pooled
  replication    <- indep_pooled * core_pooled
  strength       <- diversity_mult * da_pooled
  asr_pooled     <- max(0, min(1, replication * strength * cg_pooled))

  out <- data.frame(
    asr_path_score      = asr_pooled,
    independence        = indep_pooled,
    mrca_diversity      = div_pooled,
    derived_agreement   = da_pooled,
    conservation_gate   = cg_pooled,
    core                = core_pooled,
    core_perside_pooled = core_row_pooled,
    n_hypotheses        = length(hyps),
    supporting_hypotheses = paste(sort(hyps), collapse = ","),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(K)) out[[paste0("pooled_domain_", i)]] <- c_vec[i]
  out
}

#' Apply domain-pooling across a whole disambiguation data.frame.
#'
#' Groups by (Gene, Position, caap_group) and replaces the per-hypothesis
#' asr columns with their pooled value, collapsing to one row per group.
#' Columns not touched by pooling are taken from the highest-asr hypothesis row
#' (deterministic display/gating pick, same spirit as section 2g's
#' scheme_priority first()).
#'
#' @param df full disambiguation data.frame (needs Gene, Position, caap_group,
#'   hyp_id + the asr columns + mrca_<i>_path_score/_node).
#' @param hyp_pairs_path path to contrast_hypotheses_pairs.tsv (or NULL/NO_*).
#' @return df with one row per (Gene, Position, caap_group); adds n_hypotheses,
#'   core_perside_pooled, pooled_domain_<i>.
apply_fop_pooling <- function(df, hyp_pairs_path = NULL) {
  path_cols <- grep("^mrca_\\d+_path_score$", names(df), value = TRUE)
  path_cols <- path_cols[order(as.integer(sub("^mrca_(\\d+)_path_score$", "\\1", path_cols)))]
  node_cols <- sub("_path_score$", "_node", path_cols)

  if (!"hyp_id" %in% names(df)) df$hyp_id <- NA_character_
  for (col in c("asr_path_score", "independence", "mrca_diversity",
                "derived_agreement", "conservation_gate", "core")) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  has_fop <- df %>%
    group_by(Gene, Position, caap_group) %>%
    summarise(n = dplyr::n_distinct(hyp_id[!is.na(hyp_id) & nzchar(hyp_id)]),
              .groups = "drop") %>%
    summarise(any(n > 1)) %>% pull()

  if (!isTRUE(has_fop) || length(path_cols) == 0) {
    # nothing to pool — still add the descriptor columns for a stable schema
    df <- df %>%
      group_by(Gene, Position, caap_group) %>%
      mutate(
        n_hypotheses = dplyr::n_distinct(hyp_id[!is.na(hyp_id) & nzchar(hyp_id)]),
        supporting_hypotheses = paste(sort(unique(hyp_id[!is.na(hyp_id) & nzchar(hyp_id)])), collapse = ",")
      ) %>%
      ungroup()
    if (!"core_perside_pooled" %in% names(df)) df$core_perside_pooled <- df$core
    return(df)
  }

  hyp_pairs <- read_hypothesis_pairs(hyp_pairs_path)

  df <- df %>% arrange(desc(ifelse(is.finite(asr_path_score), asr_path_score, -Inf)))
  df$.grp <- paste(df$Gene, df$Position, df$caap_group, sep = "\r")

  # Only groups with >1 distinct hypothesis need the pooling loop. Everything
  # else (single contrast, or a position only one hypothesis reached) collapses
  # to its highest-asr row with the descriptor columns filled in — vectorised.
  n_hyp_by_grp <- tapply(df$hyp_id, df$.grp,
                         function(h) length(unique(h[!is.na(h) & nzchar(h)])))
  multi_grps <- names(n_hyp_by_grp)[n_hyp_by_grp > 1]

  hyp_by_grp <- tapply(df$hyp_id, df$.grp,
                       function(h) paste(sort(unique(h[!is.na(h) & nzchar(h)])), collapse = ","))
  simple <- df[!df$.grp %in% multi_grps, , drop = FALSE]
  if (nrow(simple) > 0) {
    simple <- simple[!duplicated(simple$.grp), , drop = FALSE]  # keep highest-asr row
    simple$n_hypotheses <- as.integer(n_hyp_by_grp[simple$.grp])
    simple$supporting_hypotheses <- unname(hyp_by_grp[simple$.grp])
    simple$core_perside_pooled <- simple$core
  }

  pooled_list <- list()
  if (length(multi_grps) > 0) {
    split_multi <- split(df[df$.grp %in% multi_grps, , drop = FALSE],
                         df$.grp[df$.grp %in% multi_grps])
    pooled_list <- lapply(split_multi, function(g) {
      base <- g[1, , drop = FALSE]
      pooled <- pool_group(g, path_cols, node_cols, hyp_pairs)
      for (col in names(pooled)) base[[col]] <- pooled[[col]][1]
      base
    })
  }

  out <- dplyr::bind_rows(c(list(simple), pooled_list))
  out$.grp <- NULL
  out
}
