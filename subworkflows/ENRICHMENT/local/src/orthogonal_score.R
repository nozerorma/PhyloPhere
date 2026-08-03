#!/usr/bin/env Rscript
# =============================================================================
# orthogonal_score.R — Lancaster-weighted "orthogonal support" composite
# =============================================================================
# Ranks a leading-edge gene set not by its own module's p-value, but by how
# much INDEPENDENT evidence backs it: the module's own directional percentile
# concentration (how extreme its own genes rank), plus cross-module/cross-source
# corroboration (FADE, RER, Accumulation, PrimateAI-3D, COSMIC).
#
# Design (validated by hand this session against a real 21-term CAAS `top`
# leading-edge table — "Extracellular matrix organization" -> signal ~= 11.2):
#   * Every component is a one-sided hypergeometric enrichment test: does this
#     term's leading edge carry MORE of a given flag than the flag's own
#     full-universe background rate would predict.
#   * The percentile trio (top10/top5/top1, or the bottom equivalents) is
#     NESTED by construction (top1 subset of top5 subset of top10) — combining
#     all three at full weight would triple-count the same signal. Each gets
#     df=2/3 (Lancaster's weighted Fisher method) so together they contribute
#     exactly one test's worth, while still rewarding deep concentration (a
#     gene in top1 still strengthens all three components).
#   * Every other flag (cross-module or cross-database — genuinely independent
#     data sources) gets full weight, df=2, the standard Fisher unit.
#   * All components are combined via Lancaster's weighted Fisher method:
#     independent partial-chi-squares with df=w_i sum to a chi-square with
#     df=sum(w_i).
#
# gate_sig is deliberately NOT included (gate_fdr no longer exists — see
# SCORING's scoring_compute.R for why it was removed entirely): checked against a real
# leading-edge table this session, they sit at 87-100% for essentially every
# term regardless of significance rank (every leading-edge gene already
# cleared the gate before it could carry a nonzero score at all) — a ceiling
# effect with no discriminating power in this population.
# =============================================================================

# ── Per-gene percentile flags (directional) ──────────────────────────────────
# Verbatim copy of fcs_enrich.R's fcs_percentile_flags() -- duplicated here
# (rather than sourced from fcs_enrich.R) so callers that only need this one
# small function (e.g. 15.Comparison_report.Rmd's Interesting Genes/Positions
# tables) aren't forced to also load fcs_enrich.R's heavier RERconverge
# dependency, which they never otherwise use. Keep both copies in sync if
# this logic ever changes.
# From score_top / score_bottom columns, derive cumulative top-X% membership
# among the NON-ZERO (signal) genes of each directional ranking:
#   pct_top25  = gene is in the top 25% of the top ranking (⊇ top10 ⊇ top5 ⊇ top1)
#   pct_top10  = gene is in the top 10% of the top ranking      (⊇ top5 ⊇ top1)
#   pct_bottom10 = gene is in the top 10% of the bottom ranking (strongly bottom)
# Returns a per-gene data.frame of logical columns (gene + pct_*). Quartet, not
# a trio -- 25% is a real, tested tier elsewhere in this pipeline (posenrich's
# own Fisher test runs at char_fracs 0.25/0.10/0.05/0.01, AMI's DOMINO module
# search runs on top25/bottom25 gene lists), so leaving it out of the
# percentile-membership criteria here (and downstream in the Interesting
# Genes/Positions tables) would be inconsistent with what's actually tested.
#
# gene_lists_dir (optional): SCORING's published gene_lists/slice_<dir><pct>.tsv
# (scoring_compute.R) -- when supplied AND all 4 files for a direction exist,
# membership comes from there instead of re-ranking score_top/score_bottom
# here, so this and every other report/module (AMI already does) agree on
# gene percentile membership by construction. Only valid when score_top/
# score_bottom actually ARE CAAS's gene_caas_score_top_all/bottom_all (true
# for CAAS's own report and for 15.Comparison_report.Rmd's gene_attr_df,
# NEVER for RER's report, which reuses this same function with its own
# ranking under the same column names) -- callers must only pass
# gene_lists_dir when that holds. Falls back to re-ranking (all 4 files for a
# direction, not a mix) when gene_lists_dir is NULL or a file is missing.
fcs_percentile_flags <- function(stats, breaks = c(0.25, 0.10, 0.05, 0.01), gene_lists_dir = NULL) {
  pct_of <- function(v) {
    ok <- !is.na(v) & v > 0
    p <- rep(NA_real_, length(v))
    if (any(ok)) p[ok] <- 1 - (rank(v[ok], ties.method = "max") - 1) / sum(ok)
    p
  }
  .slice_genes <- function(direction, pct) {
    if (is.null(gene_lists_dir)) return(NULL)
    f <- file.path(gene_lists_dir, sprintf("slice_%s%d.tsv", direction, pct))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !"Gene" %in% names(d)) return(NULL)
    unique(as.character(d$Gene))
  }
  .published_quartet <- function(direction) {
    g25 <- .slice_genes(direction, 25); g10 <- .slice_genes(direction, 10)
    g5  <- .slice_genes(direction, 5);  g1  <- .slice_genes(direction, 1)
    if (is.null(g25) || is.null(g10) || is.null(g5) || is.null(g1)) return(NULL)
    list(g25 = g25, g10 = g10, g5 = g5, g1 = g1)
  }

  out <- data.frame(gene = stats$gene, stringsAsFactors = FALSE)
  if ("score_top" %in% names(stats)) {
    pub <- .published_quartet("top")
    if (!is.null(pub)) {
      out$pct_top25 <- out$gene %in% pub$g25
      out$pct_top10 <- out$gene %in% pub$g10
      out$pct_top5  <- out$gene %in% pub$g5
      out$pct_top1  <- out$gene %in% pub$g1
    } else {
      pt <- pct_of(stats$score_top)
      out$pct_top25 <- !is.na(pt) & pt <= breaks[1]
      out$pct_top10 <- !is.na(pt) & pt <= breaks[2]
      out$pct_top5  <- !is.na(pt) & pt <= breaks[3]
      out$pct_top1  <- !is.na(pt) & pt <= breaks[4]
    }
  }
  if ("score_bottom" %in% names(stats)) {
    pub <- .published_quartet("bottom")
    if (!is.null(pub)) {
      out$pct_bottom25 <- out$gene %in% pub$g25
      out$pct_bottom10 <- out$gene %in% pub$g10
      out$pct_bottom5  <- out$gene %in% pub$g5
      out$pct_bottom1  <- out$gene %in% pub$g1
    } else {
      pb <- pct_of(stats$score_bottom)
      out$pct_bottom25 <- !is.na(pb) & pb <= breaks[1]
      out$pct_bottom10 <- !is.na(pb) & pb <= breaks[2]
      out$pct_bottom5  <- !is.na(pb) & pb <= breaks[3]
      out$pct_bottom1  <- !is.na(pb) & pb <= breaks[4]
    }
  }
  out
}

# ── One-sided hypergeometric enrichment, vectorized over rows ────────────────
# k, n_draw: vectors, one entry per (ranking, database, pathway) term.
# rate, n_universe: scalars (background rate and universe size for this flag).
orth_hyperp <- function(k, n_draw, rate, n_universe) {
  K <- round(rate * n_universe)
  phyper(k - 1, K, n_universe - K, n_draw, lower.tail = FALSE)
}

# ── Background (full-universe) TRUE-rate for a set of flag columns ───────────
# df: one row per universe gene (e.g. fcs_stats.tsv, or gene_scores.tsv).
# cols: character vector of column names (logical, or TRUE/"TRUE"/1-coercible).
orth_background_rates <- function(df, cols) {
  istrue <- function(x) x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  rates <- vapply(cols, function(cl) {
    if (!cl %in% names(df)) return(NA_real_)
    mean(istrue(df[[cl]]), na.rm = TRUE)
  }, numeric(1))
  names(rates) <- cols
  rates
}

# ── Composite orthogonal-support score per (ranking, database, pathway) ──────
# le: leading-edge long table, one row per (ranking, database, pathway, gene)
#   membership, with the flag columns present as TRUE/FALSE (or NA/absent).
# bg: named numeric vector of background rates (from orth_background_rates()),
#   must have an entry for every name in percentile_cols/cross_cols that is
#   actually present in `le`.
# percentile_cols: length-4 character vector naming the nested percentile-
#   quartet columns for this ranking's direction, e.g.
#   c("pct_top25","pct_top10","pct_top5","pct_top1") or
#   c("pct_bottom25","pct_bottom10","pct_bottom5","pct_bottom1"). NULL/empty
#   to skip.
# cross_cols: character vector of independent flag columns to combine at full
#   weight, e.g. c("flag_fade_top","flag_rer_acc","flag_rer_decc","flag_accum",
#   "flag_pai3d_pathogenic","flag_cosmic_overlap"). Columns absent from `le`
#   are silently skipped, so per-module reports missing some annotation still
#   work — the composite just combines whatever independent evidence exists.
# n_universe: total universe size used to build `bg`.
# Returns a tibble: ranking, database, pathway, n_le, signal, p_composite.
orthogonal_score <- function(le, bg, percentile_cols = character(0),
                              cross_cols = character(0), n_universe) {
  istrue <- function(x) x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  percentile_cols <- intersect(percentile_cols, names(le))
  cross_cols      <- intersect(cross_cols, names(le))
  all_cols        <- c(percentile_cols, cross_cols)

  base <- le %>%
    dplyr::group_by(ranking, database, pathway) %>%
    dplyr::summarise(n_le = dplyr::n(), .groups = "drop")

  if (length(all_cols) == 0) {
    return(base %>% dplyr::mutate(signal = NA_real_, p_composite = NA_real_))
  }

  counts <- le %>%
    dplyr::group_by(ranking, database, pathway) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(all_cols), ~ sum(istrue(.x)), .names = "k_{.col}"),
      .groups = "drop"
    )

  wide <- dplyr::left_join(base, counts, by = c("ranking", "database", "pathway"))

  Q       <- rep(0, nrow(wide))
  total_w <- 0
  add_component <- function(Q, total_w, cl, w) {
    k <- wide[[paste0("k_", cl)]]
    p <- orth_hyperp(k, wide$n_le, bg[[cl]], n_universe)
    Q <- Q + qchisq(pmax(p, 1e-300), df = w, lower.tail = FALSE)
    list(Q = Q, total_w = total_w + w)
  }
  # Nested-family discount: 4 percentile columns (25/10/5/1%) each get 1/2 the
  # normal df=2 weight, so the family totals 4 x 1/2 = 2 -- the same weight as
  # ONE undiscounted flag, avoiding a single evidence source (CAAS percentile
  # concentration) from quadruple-counting just because it's expressed as 4
  # nested cutoffs.
  for (cl in percentile_cols) {
    r <- add_component(Q, total_w, cl, w = 1/2); Q <- r$Q; total_w <- r$total_w
  }
  for (cl in cross_cols) {
    r <- add_component(Q, total_w, cl, w = 2); Q <- r$Q; total_w <- r$total_w
  }

  wide$p_composite <- pchisq(Q, df = total_w, lower.tail = FALSE)
  wide$signal      <- -log10(pmax(wide$p_composite, 1e-300))
  wide %>% dplyr::select(ranking, database, pathway, n_le, signal, p_composite)
}

# ── Surprisal-weighted composite (single-item analogue of orthogonal_score) ──
# orthogonal_score() runs a proper hypergeometric enrichment test, which
# needs a SET of members drawn from a background to be
# meaningful (a leading edge, a pathway). Ranking individual items -- one
# position, one gene, one network node -- by "how many boolean criteria does
# it meet" has no such set, so a naive unweighted count is the only thing
# that degenerates gracefully to n=1 -- but it also treats every criterion as
# equally informative regardless of how common it is (a lenient significance
# gate that 90% of items pass counts the same as a rare COSMIC overlap that
# 2% do). This weights each TRUE flag by its own surprisal / self-information
# (-log2 of its background TRUE-rate) instead: common flags contribute little,
# rare ones contribute a lot, and it's additive so it stays meaningful for a
# single item. Not a p-value or a significance test -- a ranking heuristic,
# same spirit as the orthogonal composite but usable at n=1.
# df: one row per item. cols: flag column names to weight (columns absent
#   from df or bg_rates are silently skipped). bg_rates: named numeric vector
#   of background TRUE-rates for those columns (e.g. from
#   orth_background_rates()), same names as cols.
# Returns a numeric vector, one score per row of df (0 when no flags apply).
surprisal_weighted_scores <- function(df, cols, bg_rates) {
  istrue <- function(x) x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  cols <- intersect(cols, names(df))
  cols <- intersect(cols, names(bg_rates))
  if (length(cols) == 0 || nrow(df) == 0) return(rep(0, nrow(df)))
  w <- -log2(pmax(bg_rates[cols], 1e-6))
  mat <- vapply(cols, function(cl) istrue(df[[cl]]), logical(nrow(df)))
  if (!is.matrix(mat)) mat <- matrix(mat, nrow = nrow(df), dimnames = list(NULL, cols))
  as.vector(mat %*% w)
}

# ── Standard cross-source flag set (independent of a report's own direction) ─
# Reused verbatim by the Comparison report and both surviving FCS reports so
# the composite means the same thing everywhere it appears.
# flag_accum_top/flag_accum_bottom are the direction-restricted accumulation
# tests (accumulation_top_*/accumulation_bottom_* pooled separately), same
# convention as flag_fade_top/flag_fade_bottom -- NOT the non-directional
# flag_accum (accumulation_all_* pooled), which is a genuinely different test
# and is deliberately excluded here (see orth_cols_for_direction()'s "global"
# fallback below, which is the only place it belongs).
ORTH_CROSS_COLS <- c("flag_fade_top", "flag_fade_bottom", "flag_rer_acc",
                      "flag_rer_decc", "flag_accum_top", "flag_accum_bottom",
                      "flag_pai3d_pathogenic", "flag_cosmic_overlap")

# Direction-matched percentile quartet + cross-source flags for one ranking.
# direction: "top", "bottom", or anything else (treated as non-directional —
#   no percentile quartet, cross-source flags only, and both directions' FADE/
#   Accum flags apply since there's no direction of our own to match against).
orth_cols_for_direction <- function(direction) {
  pct <- switch(direction,
    top    = c("pct_top25", "pct_top10", "pct_top5", "pct_top1"),
    bottom = c("pct_bottom25", "pct_bottom10", "pct_bottom5", "pct_bottom1"),
    character(0)
  )
  cross <- switch(direction,
    top    = setdiff(ORTH_CROSS_COLS, c("flag_fade_bottom", "flag_accum_bottom")),
    bottom = setdiff(ORTH_CROSS_COLS, c("flag_fade_top", "flag_accum_top")),
    ORTH_CROSS_COLS
  )
  list(percentile_cols = pct, cross_cols = cross)
}

# =============================================================================
# Mixed-granularity Lancaster combination (posenrich position-level composite)
# =============================================================================
# orthogonal_score() above assumes ONE shared (n_draw, n_universe) for every
# component -- correct when every flag is measured over
# the same population (e.g. all gene-level, as in the FCS reports). Some
# composites legitimately need to mix components measured at DIFFERENT grains
# within the same term: posenrich's Overall dotplot tests position-native
# evidence (CAAS top-fraction recurrence, FADE-site overlap) against the
# position universe, alongside gene-inherited evidence (RER, Accumulation,
# FCS/AMI membership) tested at gene granularity against the gene universe --
# testing gene-inherited flags at position resolution would be pseudo-
# replication, since every position in one gene shares that gene's flag value
# and is not an independent draw with respect to it (see
# 14.Position_enrichment_report.Rmd's Overall dotplot for the full reasoning).
# Lancaster's method still applies: it combines independent p-values
# regardless of what population/grain each was computed over, as long as each
# component is a properly-formed one-sided hypergeometric test in its own right.
#
# components: a list of lists, each with k (successes drawn), n_draw (draw
#   size), rate (background TRUE-rate for this flag, in the SAME units as
#   n_draw/n_universe), n_universe (population size for this specific
#   component), and weight (chi-square df to contribute -- 2 for a standalone
#   independent flag, a fraction like 2/length(family) for a member of a
#   nested family so the whole family sums to one test's worth, same
#   convention as the percentile trio above).
# Returns list(p_composite, signal, total_w).
lancaster_combine <- function(components) {
  if (length(components) == 0) return(list(p_composite = NA_real_, signal = NA_real_, total_w = 0))
  Q <- 0
  total_w <- 0
  for (comp in components) {
    if (is.na(comp$k) || is.na(comp$n_draw) || is.na(comp$rate) || is.na(comp$n_universe) ||
        comp$n_draw <= 0 || comp$n_universe <= 0) next
    p <- orth_hyperp(comp$k, comp$n_draw, comp$rate, comp$n_universe)
    Q <- Q + qchisq(pmax(p, 1e-300), df = comp$weight, lower.tail = FALSE)
    total_w <- total_w + comp$weight
  }
  if (total_w == 0) return(list(p_composite = NA_real_, signal = NA_real_, total_w = 0))
  p_composite <- pchisq(Q, df = total_w, lower.tail = FALSE)
  list(p_composite = p_composite, signal = -log10(pmax(p_composite, 1e-300)), total_w = total_w)
}

# ── Position-granularity hybrid composite for one posenrich term ─────────────
# le_rk: posenrich_leading_edge.tsv rows for ONE (database, pathway), ONE
#   ranking/direction -- columns top_frac, gene, gene_position (long, one row
#   per position member per top_frac cutoff it was significant at).
# top_fracs: the nested top-fraction cutoffs to treat as the trio-analogue
#   (e.g. c(0.10, 0.05, 0.01) -- 0.25 excluded the same way genes' trio omits
#   the "everything" baseline). Each gets weight 2/length(top_fracs); rate is
#   the cutoff's own literal value (a position's chance of landing in the top
#   X% of the SCORED population is X%, by construction of the cutoff itself --
#   same reasoning as why pct_top10's background rate is ~10%).
# n_universe_scored: the SCORED-position population size for this direction
#   (posenrich_enrich.py's `n_scored` -- positions with CAAS_score > 0, ranked
#   to build top_frac cutoffs). NOT n_universe_positions -- top_frac cutoffs
#   are defined as fractions of the scored subset, which is far smaller than
#   the full honest tested-position background (most tested positions carry no
#   CAAS signal at all), so using the honest background here would test the
#   wrong null (a position's true chance of being top-10%-scored is ~10% of
#   n_scored, not ~10% of the much larger n_universe_positions).
# gene_flag_cols: named character vector, flag_<x> columns in fcs_stats_df to
#   test at GENE granularity (distinct genes among the term's member positions),
#   weight 2 each.
# fcs_stats_df, n_universe_genes: gene-level background + universe for
#   gene_flag_cols.
# fade_df: Gene/Position tibble of this direction's FADE-significant sites (or
#   NULL). Tested at POSITION granularity against the full honest background
#   (FADE sites aren't restricted to CAAS-scored positions), weight 2.
# n_universe_positions, fade_bg_rate: honest tested-position background + FADE's
#   own background rate (n_fade_sites / n_universe_positions) for that test.
# Returns a one-row tibble: n_le (positions), n_genes (distinct genes), signal, p_composite.
posenrich_position_signal <- function(le_dp, top_fracs, n_universe_scored, gene_flag_cols, fcs_stats_df,
                                        n_universe_genes, fade_df, n_universe_positions,
                                        fade_bg_rate) {
  members <- le_dp %>% dplyr::distinct(gene, gene_position)
  n_le    <- nrow(members)
  n_genes <- dplyr::n_distinct(members$gene)
  if (n_le == 0) {
    return(tibble::tibble(n_le = 0L, n_genes = 0L, signal = NA_real_, p_composite = NA_real_))
  }

  components <- list()

  # Nested top-fraction recurrence (percentile-trio analogue), tested against
  # the SCORED population -- see n_universe_scored above.
  for (tf in top_fracs) {
    tf_members <- le_dp %>% dplyr::filter(top_frac == tf) %>% dplyr::distinct(gene_position)
    k <- length(intersect(members$gene_position, tf_members$gene_position))
    components[[paste0("tf_", tf)]] <- list(
      k = k, n_draw = n_le, rate = tf, n_universe = n_universe_scored,
      weight = 2 / length(top_fracs)
    )
  }

  # Gene-inherited cross flags, distinct-gene counted (pseudoreplication guard).
  if (length(gene_flag_cols) > 0 && !is.null(fcs_stats_df)) {
    gene_flags <- fcs_stats_df %>%
      dplyr::filter(gene %in% members$gene) %>%
      dplyr::select(gene, dplyr::any_of(unname(gene_flag_cols)))
    for (nm in names(gene_flag_cols)) {
      cl <- gene_flag_cols[[nm]]
      if (!cl %in% names(gene_flags)) next
      k <- sum(gene_flags[[cl]] %in% TRUE, na.rm = TRUE)
      rate <- mean(fcs_stats_df[[cl]] %in% TRUE, na.rm = TRUE)
      components[[nm]] <- list(k = k, n_draw = n_genes, rate = rate,
                                n_universe = n_universe_genes, weight = 2)
    }
  }

  # Position-native FADE-site overlap, position counted.
  if (!is.null(fade_df) && !is.na(fade_bg_rate)) {
    fade_hits <- members %>%
      dplyr::mutate(Position = suppressWarnings(as.integer(sub(".*:", "", gene_position)))) %>%
      dplyr::inner_join(fade_df %>% dplyr::rename(gene = Gene) %>% dplyr::distinct(gene, Position),
                        by = c("gene", "Position"))
    components[["fade_pos"]] <- list(k = nrow(fade_hits), n_draw = n_le, rate = fade_bg_rate,
                                      n_universe = n_universe_positions, weight = 2)
  }

  lc <- lancaster_combine(components)
  tibble::tibble(n_le = n_le, n_genes = n_genes, signal = lc$signal, p_composite = lc$p_composite)
}
