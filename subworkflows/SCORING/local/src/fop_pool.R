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
#
# TWO-JOB PSS WEIGHTING (refined 2026-09-03). One PSS number per hypothesis is
# the wrong instrument for both consumers, because every FOP hypothesis shares
# the harvest's globally weakest domain (min-across-domains is near-constant
# across H1..Hn -> the weighting silently collapses to equal-weight). So the
# PSS split does two distinct jobs:
#   Job A -- the per-domain c_i pool: weight each DISTINCT deduped MRCA node by
#     THAT candidate pair's own PSS in domain i (lookup (hyp_id, domain) -> pss
#     from contrast_hypotheses_pairs.tsv; per node, max over the rows sharing
#     it, since PSS is a property of the species pair). A pair with strong own
#     isolation still drives its domain even when it rides in a hypothesis whose
#     other domain is weak.
#   Job B -- the axis pools (independence / mrca_diversity / derived_agreement /
#     row-level core): per-hypothesis weight = MEAN pair PSS across the
#     hypothesis's K domains (a hypothesis's overall credibility, not its
#     weakest link).
#
# CONSERVATION_GATE (POINT 2, refined 2026-09-03). A conserved pair is a PURE
# conservation_gate input in path_scores.py (scored by conservation-to-root, then
# `continue` — never touches core/independence/diversity/derived_agreement).
# Pooling the per-hypothesis gates counts a conserved pair shared by several
# hypotheses once PER hypothesis (asymmetric with the node-deduped c_i pool) and
# breaks the 0.5+0.5*mean transform when conserved-pair counts differ across
# hypotheses. So when the flat TSV carries a parallel conserved_<j>_node /
# conserved_<j>_cons block, the gate is rebuilt from the DISTINCT conserved pairs
# of the harvest group (dedup by conserved_<j>_node):
#     cg = 0.5 + 0.5 * wmean(cons over distinct conserved pairs, w)
# w = that pair's own PSS in the domain whose mrca_<i>_node it matches on the
# same row (else equal weight); no distinct conserved pairs -> cg = 1.0. Older
# inputs without the conserved_<j>_* columns keep the per-hypothesis-gate wmean.
# Both keep the equal-weight fallback when no PSS file / all PSS non-finite.
#   * recombine c_1..c_K into the position score with the SAME algebra as
#     path_scores.py::compute_asr_path_score:
#       core        = P(>=2 domains carry a clean change)   [inclusion-exclusion]
#       replication = independence * core
#       strength    = (0.75 + 0.25 * diversity) * derived_agreement
#       asr         = replication * strength * conservation_gate
#     independence / diversity / derived_agreement are NOT re-derivable from the
#     flat TSV (they need the LAC tree structure), so the per-hypothesis values
#     are pooled by the same PSS weights. Voronoi domains
#     are disjoint, so cross-domain independence is ~1 by construction and the
#     pooled independence only ever mildly discounts.
#   * DIRECTIONAL core (FIXED 2026-09-04). `mrca_<i>_path_score` is only the
#     per-pair TOP/BOTTOM *average* — pooling that single number per domain
#     into one P(>=2 domains) let a domain whose pair changed on the top
#     phenotype side and an unrelated domain whose pair changed on the bottom
#     side count as "2 domains agreeing", when they are opposite-direction
#     events (exactly what path_scores.py's own core_top/core_bottom split
#     exists to rule out — see disambiguate_single.py's directional core note).
#     The disambiguation writer now also emits `mrca_<i>_top_path_score` /
#     `mrca_<i>_bot_path_score` (the un-averaged per-side scores); when present,
#     `pool_group` pools each side SEPARATELY (dedup by MRCA node, same PSS
#     weighting as above) into c_top_i / c_bot_i and combines them exactly like
#     compute_asr_path_score: `core = 1 - (1-core_top)*(1-core_bottom)`. Older
#     `filtered_discovery.tsv` inputs without those columns fall back to the
#     collapsed, direction-blind pool over `mrca_<i>_path_score` (a known
#     approximation, kept only so a stale file still scores instead of erroring).
#     The per-hypothesis `core` column (pooled, below, as `core_perside_pooled`)
#     is carried alongside as a cross-check either way.
#
# POINT 3 (harvest-wide, scheme-resolved derived-residue reconciliation, added
# 2026-09-03). `derived_agreement` in path_scores.py is measured WITHIN one
# hypothesis: the fraction of that hypothesis's changed pairs landing on its own
# plurality derived residue. Two FOP hypotheses can each be internally unanimous
# (da = 1.0) yet land on DIFFERENT residues from one another; pooling their da
# values (the old `.wmean(df$derived_agreement, row_w)`) never sees that split.
#
# When the flat TSV carries the raw-residue block (mrca_<i>_anc_aa / _top_aa /
# _bot_aa — the un-encoded ancestral and per-side derived residues), da is
# recomputed HARVEST-WIDE: over the pooled, node-deduplicated changed-pair set of
# the group, per side, applying the EXACT path_scores.py plurality logic, but on
# residues encoded under the ROW'S OWN scheme (`encode_aa_r`). A genuinely split
# position (US: V/I/L) then gets a low harvest-wide da under US and 1.0 under a
# scheme that co-encodes those residues — and §2g's mean(caas_row) over schemes
# does the reconciliation with NO new multiplier. Older inputs without the *_aa
# columns keep the per-hypothesis `.wmean` fallback.
#
# apply_fop_pooling also emits four POSITION-LEVEL descriptor columns (raw-AA,
# joined back by (Gene, Position), so they sit OUTSIDE the per-caap_group
# grouping): `derived_residues` (compact "<anc>/<sorted derived>"),
# `top_residue_support` / `bottom_residue_support` (per raw derived residue,
# distinct-pair count, one column per phenotype side), and
# `convergence_schemes` — the SUBSET of US,GS1,GS2,GS3,GS4 whose harvest-wide da
# on the raw pooled pairs is >= tau (default 0.8). `convergence_schemes` is a
# SET/PROFILE, never an ordinal: the 5 schemes are non-nested partitions along
# different physicochemical axes (see aa_grouping.R).
#
# Single-hypothesis / non-FOP input is a no-op: the one row passes through with
# its asr_path_score untouched (verified by test_fop_pool.R).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Grouping schemes for the harvest-wide, per-scheme derived_agreement rebuild.
.aa_grouping_src <- file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
  "aa_grouping.R")
if (!file.exists(.aa_grouping_src)) .aa_grouping_src <- "aa_grouping.R"
if (file.exists(.aa_grouping_src)) {
  source(.aa_grouping_src)
} else if (!exists("encode_aa_r")) {
  # Degrade gracefully: without the scheme maps, harvest-wide da can only be
  # computed under identity (US-equivalent). Fallback path still works.
  AA_SCHEME_NAMES <- c("US", "GS1", "GS2", "GS3", "GS4")
  encode_aa_r <- function(aa, scheme) {
    if (is.null(aa) || length(aa) == 0) return(NA_character_)
    r <- toupper(trimws(as.character(aa)[1]))
    if (is.na(r) || !nzchar(r)) NA_character_ else r
  }
}

#' Harvest-wide, per-scheme derived_agreement over a group's DISTINCT changed pairs.
#'
#' Mirrors path_scores.py::compute_asr_path_score's derived_agreement EXACTLY,
#' but over the pooled, node-deduplicated changed-pair set of the whole harvest
#' group instead of one hypothesis, and on residues encoded under `scheme`.
#'
#' @param changed_df data.frame with columns node, side ("top"/"bot"), raw_aa —
#'   one row per (distinct MRCA node, side) that carried a change. Callers dedup
#'   by (node, side) before passing (a node's derived residue on a side is a
#'   fixed observed fact).
#' @param scheme one of AA_SCHEME_NAMES.
#' @return derived_agreement in (0, 1]. 1.0 when no side has >= 2 distinct
#'   changed pairs (matches path_scores.py).
rebuild_derived_agreement <- function(changed_df, scheme) {
  if (is.null(changed_df) || nrow(changed_df) == 0) return(1.0)
  concentrations <- c()
  for (sd in c("top", "bot")) {
    rs <- changed_df$raw_aa[changed_df$side == sd]
    rs <- rs[!is.na(rs) & nzchar(rs)]
    if (length(rs) < 2) next
    enc <- vapply(rs, encode_aa_r, character(1), scheme = scheme)
    tab <- table(enc)
    concentrations <- c(concentrations, max(tab) / length(enc))
  }
  if (length(concentrations) == 0) return(1.0)
  mean(concentrations)
}

#' Distinct (node, side, raw_aa) changed-pair records from a harvest group's rows.
#'
#' Scans every mrca_<i>_top_aa / mrca_<i>_bot_aa cell across all rows, pairing it
#' with the matching mrca_<i>_node. Deduplicates by (node, side): the same
#' physical pair recurring across hypotheses / schemes is one changed pair.
#' Voronoi domains are disjoint, so a node identifies its pair unambiguously.
.collect_changed_pairs <- function(df, node_cols, top_cols, bot_cols) {
  recs <- list()
  add <- function(node, side, aa) {
    if (is.na(node) || !nzchar(as.character(node))) return(invisible())
    aa <- toupper(trimws(as.character(aa)))
    if (is.na(aa) || !nzchar(aa)) return(invisible())
    recs[[length(recs) + 1L]] <<- data.frame(
      node = as.character(node), side = side, raw_aa = aa, stringsAsFactors = FALSE)
  }
  K <- length(node_cols)
  for (i in seq_len(K)) {
    ncol <- node_cols[i]
    if (!(ncol %in% names(df))) next
    nodes <- as.character(df[[ncol]])
    if (i <= length(top_cols) && top_cols[i] %in% names(df)) {
      tv <- as.character(df[[top_cols[i]]])
      for (r in seq_along(nodes)) add(nodes[r], "top", tv[r])
    }
    if (i <= length(bot_cols) && bot_cols[i] %in% names(df)) {
      bv <- as.character(df[[bot_cols[i]]])
      for (r in seq_along(nodes)) add(nodes[r], "bot", bv[r])
    }
  }
  if (length(recs) == 0) {
    return(data.frame(node = character(0), side = character(0),
                      raw_aa = character(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, recs)
  out[!duplicated(out[, c("node", "side")]), , drop = FALSE]
}

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
#' @param top_path_cols / bot_path_cols character vectors of the
#'   mrca_<i>_top_path_score / mrca_<i>_bot_path_score column names (same order
#'   as path_cols), when present — the un-averaged per-side scores `core` is
#'   actually built from (directional core_top/core_bottom). Absent ->
#'   `core` falls back to the collapsed, direction-blind pool over path_cols.
#' @return one-row data.frame: asr_path_score, independence, mrca_diversity,
#'   derived_agreement, conservation_gate, core, core_perside_pooled,
#'   n_hypotheses, plus c_1..c_K as pooled_domain_<i>.
pool_group <- function(df, path_cols, node_cols, hyp_pairs = NULL,
                       diversity_floor = 0.75,
                       conserved_node_cols = character(0),
                       conserved_cons_cols = character(0),
                       top_aa_cols = character(0),
                       bot_aa_cols = character(0),
                       top_path_cols = character(0),
                       bot_path_cols = character(0),
                       scheme = NA_character_) {
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

  # Job B weight: per-hypothesis MEAN pair PSS across its K domains (overall
  # credibility of the hypothesis, used for the axis pools). Job A weight:
  # each candidate pair's OWN PSS in a given domain, via (hyp_id, domain) -> pss.
  hyp_w <- setNames(rep(NA_real_, length(hyps)), hyps)
  pair_pss <- function(h, dom) NA_real_
  if (!is.null(hyp_pairs)) {
    for (h in hyps) {
      w <- hyp_pairs$pss_score[hyp_pairs$hyp_id == h]
      if (length(w) && any(is.finite(w))) hyp_w[h] <- mean(w[is.finite(w)])
    }
    pair_pss <- function(h, dom) {
      v <- hyp_pairs$pss_score[hyp_pairs$hyp_id == h & hyp_pairs$pair == dom]
      v <- v[is.finite(v)]
      if (length(v)) v[1] else NA_real_
    }
  }
  row_w <- unname(hyp_w[df$hyp_id])

  # ── Per-domain pool: c_i over the distinct pairs domain i contributed ──────
  # Pools ONE score column into a single scalar for domain i, deduped by MRCA
  # node (a physical pair), weighted by that pair's own PSS in domain i (Job
  # A). Shared by the display column below and by the directional core split.
  K <- length(path_cols)
  .pool_domain_col <- function(col_name, node_col, i) {
    if (is.na(col_name) || !(col_name %in% names(df))) return(NA_real_)
    s <- suppressWarnings(as.numeric(df[[col_name]]))
    nd <- if (!is.na(node_col) && node_col %in% names(df)) {
      as.character(df[[node_col]])
    } else {
      as.character(df$hyp_id)  # fall back: treat every hypothesis's pair as distinct
    }
    wr_all <- vapply(as.character(df$hyp_id), function(h) pair_pss(h, i), numeric(1))
    keep <- is.finite(s)
    if (!any(keep)) return(NA_real_)
    s <- s[keep]; nd <- nd[keep]; wr <- wr_all[keep]
    # Distinct candidate pairs in this domain (dedup by MRCA node); weight each
    # by that pair's own PSS in domain i (max over rows sharing the node), else equal.
    dd <- data.frame(nd = nd, s = s, w = wr, stringsAsFactors = FALSE) %>%
      group_by(nd) %>%
      summarise(s = dplyr::first(s),
                w = suppressWarnings(max(w, na.rm = TRUE)), .groups = "drop")
    .wmean(dd$s, dd$w)
  }

  # Display / back-compat column (`pooled_domain_<i>`): pooled from the
  # side-averaged mrca_<i>_path_score, unchanged from before this fix. NOT what
  # `core` is computed from below when the directional columns are available.
  c_vec <- vapply(seq_len(K), function(i) .pool_domain_col(path_cols[i], node_cols[i], i), numeric(1))

  # ── DIRECTIONAL core (the actual path_scores.py algebra) ────────────────────
  # mrca_<i>_path_score is the per-pair TOP/BOTTOM *average* — pooling it
  # directly into one P(>=2 domains) lets a domain whose corroborating pair
  # changed on the top phenotype side and a different domain whose pair
  # changed on the bottom side count as "2 domains agreeing", when they are
  # two unrelated, opposite-direction events. When the flat TSV carries
  # mrca_<i>_top_path_score / _bot_path_score (the un-averaged per-side
  # scores), pool each side SEPARATELY and combine exactly like
  # compute_asr_path_score's core_top/core_bottom/core. Older inputs without
  # those columns fall back to the collapsed, direction-blind pool (documented
  # known approximation) so a stale filtered_discovery.tsv still scores.
  has_side_path <- length(top_path_cols) == K && length(bot_path_cols) == K &&
    (any(top_path_cols %in% names(df)) || any(bot_path_cols %in% names(df)))
  if (has_side_path) {
    c_top <- vapply(seq_len(K), function(i) .pool_domain_col(top_path_cols[i], node_cols[i], i), numeric(1))
    c_bot <- vapply(seq_len(K), function(i) .pool_domain_col(bot_path_cols[i], node_cols[i], i), numeric(1))
    core_top    <- .p_at_least_2(c_top)
    core_bottom <- .p_at_least_2(c_bot)
    core_pooled <- 1 - (1 - core_top) * (1 - core_bottom)
  } else {
    core_pooled <- .p_at_least_2(c_vec)
  }
  indep_pooled <- .wmean(suppressWarnings(as.numeric(df$independence)), row_w)
  div_pooled   <- .wmean(suppressWarnings(as.numeric(df$mrca_diversity)), row_w)
  core_row_pooled <- .wmean(suppressWarnings(as.numeric(df$core)), row_w)

  # ── derived_agreement (POINT 3): harvest-wide, per-scheme ─────────────────────
  # `derived_agreement` in path_scores.py is a within-hypothesis fraction. Two
  # FOP hypotheses can each be unanimous yet land on DIFFERENT residues; pooling
  # their da never sees that. When the raw-residue block is present, recompute da
  # over the group's DISTINCT changed pairs (dedup by node), per side, with the
  # EXACT path_scores.py plurality logic, on residues encoded under THIS row's
  # own scheme. Older inputs (no *_aa columns) keep the per-hypothesis wmean.
  have_aa_cols <- length(top_aa_cols) > 0 &&
    (any(top_aa_cols %in% names(df)) || any(bot_aa_cols %in% names(df)))
  if (have_aa_cols && !is.na(scheme) && nzchar(scheme)) {
    cp_group <- .collect_changed_pairs(df, node_cols, top_aa_cols, bot_aa_cols)
    da_pooled <- rebuild_derived_agreement(cp_group, scheme)
  } else {
    da_pooled <- .wmean(suppressWarnings(as.numeric(df$derived_agreement)), row_w)
  }

  # ── conservation_gate (POINT 2) ─────────────────────────────────────────────
  # A conserved pair hits EXACTLY the conservation_gate branch of
  # path_scores.py::compute_asr_path_score (scored by conservation-to-root, then
  # `continue` — never a `core`/`independence`/`diversity` input). Within one
  # hypothesis: gate = 0.5 + 0.5 * mean(cons over that hypothesis's conserved
  # pairs). Pooling the already-transformed per-hypothesis gates double-counts a
  # conserved pair shared by several hypotheses (asymmetric with the node-deduped
  # c_i pool) and, when hypotheses differ in conserved-pair count,
  # mean(0.5+0.5 cons_h) != 0.5 + 0.5 mean(cons_p). So when the conserved_<j>_*
  # columns are present, rebuild the gate from the DISTINCT conserved pairs of
  # the harvest group, deduped by conserved_<j>_node.
  have_cons_cols <- length(conserved_node_cols) > 0 && length(conserved_cons_cols) > 0 &&
    any(conserved_node_cols %in% names(df)) && any(conserved_cons_cols %in% names(df))
  if (have_cons_cols) {
    cons_recs <- list()
    for (jj in seq_along(conserved_node_cols)) {
      ncol <- conserved_node_cols[jj]; ccol <- conserved_cons_cols[jj]
      if (!(ncol %in% names(df)) || !(ccol %in% names(df))) next
      nd_j <- as.character(df[[ncol]])
      cs_j <- suppressWarnings(as.numeric(df[[ccol]]))
      for (rr in seq_len(nrow(df))) {
        node <- nd_j[rr]
        if (is.na(node) || !nzchar(node) || !is.finite(cs_j[rr])) next
        # Weight: map this conserved node back to a Voronoi domain via the row's
        # own mrca_<i>_node columns; if it matches domain i, use that candidate
        # pair's own PSS in domain i (pair_pss). Otherwise equal weight.
        w_j <- NA_real_
        for (ii in seq_along(node_cols)) {
          if (node_cols[ii] %in% names(df) &&
              identical(as.character(df[[node_cols[ii]]][rr]), node)) {
            w_j <- pair_pss(as.character(df$hyp_id[rr]), ii); break
          }
        }
        cons_recs[[length(cons_recs) + 1L]] <-
          data.frame(nd = node, cons = cs_j[rr], w = w_j, stringsAsFactors = FALSE)
      }
    }
    if (length(cons_recs) == 0) {
      cg_pooled <- 1.0  # FOP structure but no conserved pairs -> novel position
    } else {
      cdd <- do.call(rbind, cons_recs) %>%
        group_by(nd) %>%
        summarise(cons = dplyr::first(cons),
                  w = suppressWarnings(max(w, na.rm = TRUE)), .groups = "drop")
      cg_pooled <- 0.5 + 0.5 * .wmean(cdd$cons, cdd$w)
    }
  } else {
    # Fallback (older inputs, no conserved_<j>_* columns): pool the per-hypothesis
    # transformed gates as before. Residual double-counting noted above.
    cg_pooled <- .wmean(suppressWarnings(as.numeric(df$conservation_gate)), row_w)
  }

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

#' Position-level `convergence_schemes` descriptor for one (Gene, Position).
#'
#' Computed across ALL rows of the position (every scheme, every hypothesis) —
#' scheme convergence is a residue-identity question, so this sits OUTSIDE the
#' per-caap_group pooling and is joined back by (Gene, Position).
#'
#'   convergence_schemes  comma-joined subset of US,GS1,GS2,GS3,GS4 whose
#'                        harvest-wide da on the raw pooled pairs is >= tau.
#'                        A SET/PROFILE (non-nested schemes), never an ordinal.
#'                        "" when none pass / no changed pairs.
#'
#' The raw-AA descriptors `derived_residues` / `top_residue_support` /
#' `bottom_residue_support` are NO LONGER computed here — they are produced once,
#' upstream, in CT_POSTPROC's `residue_descriptors.py` (which follows the
#' `change_side` field rather than scanning both sides), and `apply_fop_pooling`
#' carries them straight through from the input frame.
.position_descriptors <- function(sub_df, node_cols, top_aa_cols, bot_aa_cols, tau) {
  cp <- .collect_changed_pairs(sub_df, node_cols, top_aa_cols, bot_aa_cols)
  if (nrow(cp) == 0) return(data.frame(convergence_schemes = "", stringsAsFactors = FALSE))

  kept <- AA_SCHEME_NAMES[vapply(AA_SCHEME_NAMES, function(sc)
    isTRUE(rebuild_derived_agreement(cp, sc) >= tau), logical(1))]
  data.frame(convergence_schemes = paste(kept, collapse = ","),
             stringsAsFactors = FALSE)
}

# Position-level descriptor columns produced upstream (CT_POSTPROC
# residue_descriptors.py) and carried through here untouched.
FOP_CARRIED_DESCRIPTORS <- c("derived_residues", "top_residue_support",
                             "bottom_residue_support", "n_conserved_pairs")

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
#' @param tau harvest-wide da threshold for `convergence_schemes` (default 0.8).
#' @return df with one row per (Gene, Position, caap_group); adds n_hypotheses,
#'   core_perside_pooled, pooled_domain_<i>, `convergence_schemes` (computed here),
#'   and passes through the upstream `derived_residues` / {top,bottom}_residue_support
#'   columns (empty strings when the input frame lacks them).
apply_fop_pooling <- function(df, hyp_pairs_path = NULL, tau = 0.8) {
  path_cols <- grep("^mrca_\\d+_path_score$", names(df), value = TRUE)
  path_cols <- path_cols[order(as.integer(sub("^mrca_(\\d+)_path_score$", "\\1", path_cols)))]
  node_cols <- sub("_path_score$", "_node", path_cols)
  # Directional per-domain scores (may be absent in older filtered_discovery.tsv
  # files — pool_group falls back to the collapsed pool when so). Same order/
  # index as path_cols by construction (both derived from mrca_<i>_*).
  top_path_cols <- sub("_path_score$", "_top_path_score", path_cols)
  bot_path_cols <- sub("_path_score$", "_bot_path_score", path_cols)

  # Raw derived residue cells — still needed by .collect_changed_pairs (feeds
  # derived_agreement pooling and convergence_schemes).
  idx_of <- function(cols, pat) as.integer(sub(pat, "\\1", cols))
  top_aa_cols <- grep("^mrca_\\d+_top_aa$", names(df), value = TRUE)
  top_aa_cols <- top_aa_cols[order(idx_of(top_aa_cols, "^mrca_(\\d+)_top_aa$"))]
  bot_aa_cols <- grep("^mrca_\\d+_bot_aa$", names(df), value = TRUE)
  bot_aa_cols <- bot_aa_cols[order(idx_of(bot_aa_cols, "^mrca_(\\d+)_bot_aa$"))]

  # Parallel conserved-pair block (j ordered by pair_id ascending upstream).
  conserved_node_cols <- grep("^conserved_\\d+_node$", names(df), value = TRUE)
  conserved_node_cols <- conserved_node_cols[order(as.integer(
    sub("^conserved_(\\d+)_node$", "\\1", conserved_node_cols)))]
  conserved_cons_cols <- sub("_node$", "_cons", conserved_node_cols)

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
    # Stable schema: carried raw-AA descriptors + convergence_schemes, empty when absent.
    for (col in c(FOP_CARRIED_DESCRIPTORS, "convergence_schemes")) {
      if (!col %in% names(df)) df[[col]] <- ""
    }
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
      pooled <- pool_group(g, path_cols, node_cols, hyp_pairs,
                           conserved_node_cols = conserved_node_cols,
                           conserved_cons_cols = conserved_cons_cols,
                           top_aa_cols = top_aa_cols,
                           bot_aa_cols = bot_aa_cols,
                           top_path_cols = top_path_cols,
                           bot_path_cols = bot_path_cols,
                           scheme = as.character(g$caap_group[1]))
      for (col in names(pooled)) base[[col]] <- pooled[[col]][1]
      base
    })
  }

  out <- dplyr::bind_rows(c(list(simple), pooled_list))
  out$.grp <- NULL

  # ── Position-level convergence_schemes pass ─────────────────────────────────
  # Across ALL schemes/hypotheses of a (Gene, Position); joined back by it. The
  # raw-AA descriptors (derived_residues / {top,bottom}_residue_support) are
  # carried through from `out` untouched — they were produced upstream in
  # CT_POSTPROC's residue_descriptors.py.
  out$convergence_schemes <- NULL
  desc <- df %>%
    group_by(Gene, Position) %>%
    group_modify(~ .position_descriptors(.x, node_cols, top_aa_cols, bot_aa_cols, tau)) %>%
    ungroup()
  out <- dplyr::left_join(out, desc, by = c("Gene", "Position"))
  for (col in c(FOP_CARRIED_DESCRIPTORS, "convergence_schemes")) {
    if (!col %in% names(out)) out[[col]] <- ""
    out[[col]][is.na(out[[col]])] <- ""
  }
  out
}
