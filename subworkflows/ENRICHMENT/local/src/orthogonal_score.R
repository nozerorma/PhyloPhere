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
# gate_sig/gate_fdr are deliberately NOT included: checked against a real
# leading-edge table this session, they sit at 87-100% for essentially every
# term regardless of significance rank (every leading-edge gene already
# cleared the gate before it could carry a nonzero score at all) — a ceiling
# effect with no discriminating power in this population.
# =============================================================================

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
# percentile_cols: length-3 character vector naming the nested percentile-trio
#   columns for this ranking's direction, e.g. c("pct_top10","pct_top5","pct_top1")
#   or c("pct_bottom10","pct_bottom5","pct_bottom1"). NULL/empty to skip.
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
  for (cl in percentile_cols) {
    r <- add_component(Q, total_w, cl, w = 2/3); Q <- r$Q; total_w <- r$total_w
  }
  for (cl in cross_cols) {
    r <- add_component(Q, total_w, cl, w = 2); Q <- r$Q; total_w <- r$total_w
  }

  wide$p_composite <- pchisq(Q, df = total_w, lower.tail = FALSE)
  wide$signal      <- -log10(pmax(wide$p_composite, 1e-300))
  wide %>% dplyr::select(ranking, database, pathway, n_le, signal, p_composite)
}

# ── Composite from pre-aggregated percentages (no per-member long table) ─────
# Same statistical machinery as orthogonal_score() (one-sided hypergeometric
# per flag, Lancaster-combined at full weight df=2 each -- no percentile trio
# here, since callers of this entry point have no per-gene percentile-membership
# concept to nest), for reports whose corroboration is already computed as a
# set-level pct_<flag> (0-100) rather than carried as a per-member long table
# -- e.g. posenrich, where posenrich_enrich.py::annotate_overlap() computes
# these percentages once in Python. Recovers integer counts as
# round(pct/100 * n_le) before reusing orth_hyperp(), so this entry point
# agrees exactly with orthogonal_score() when both could apply to the same
# underlying data (same k, n_draw, rate, n_universe -> same p).
# df: one row per term, with an `n_le` column (leading-edge/overlap size) and
#   the columns named in `pct_cols`.
# pct_cols: named character vector, names(df) -> names(bg), e.g.
#   c(pct_fade_top = "flag_fade_top", pct_rer_acc = "flag_rer_acc").
# bg, n_universe: as in orthogonal_score().
# Returns `df` with `signal`/`p_composite` columns appended.
orthogonal_score_pct <- function(df, bg, pct_cols, n_universe) {
  pct_cols <- pct_cols[names(pct_cols) %in% names(df) & pct_cols %in% names(bg)]
  if (length(pct_cols) == 0) {
    return(df %>% dplyr::mutate(signal = NA_real_, p_composite = NA_real_))
  }
  Q <- rep(0, nrow(df))
  for (pc in names(pct_cols)) {
    bgcol <- pct_cols[[pc]]
    k <- round(df[[pc]] / 100 * df$n_le)
    k[is.na(k)] <- 0
    p <- orth_hyperp(k, df$n_le, bg[[bgcol]], n_universe)
    Q <- Q + qchisq(pmax(p, 1e-300), df = 2, lower.tail = FALSE)
  }
  df$p_composite <- pchisq(Q, df = 2 * length(pct_cols), lower.tail = FALSE)
  df$signal      <- -log10(pmax(df$p_composite, 1e-300))
  df
}

# ── Standard cross-source flag set (independent of a report's own direction) ─
# Reused verbatim by the Comparison report and both surviving FCS reports so
# the composite means the same thing everywhere it appears.
ORTH_CROSS_COLS <- c("flag_fade_top", "flag_fade_bottom", "flag_rer_acc",
                      "flag_rer_decc", "flag_accum",
                      "flag_pai3d_pathogenic", "flag_cosmic_overlap")

# Direction-matched percentile trio + cross-source flags for one ranking.
# direction: "top", "bottom", or anything else (treated as non-directional —
#   no percentile trio, cross-source flags only).
orth_cols_for_direction <- function(direction) {
  pct <- switch(direction,
    top    = c("pct_top10", "pct_top5", "pct_top1"),
    bottom = c("pct_bottom10", "pct_bottom5", "pct_bottom1"),
    character(0)
  )
  cross <- switch(direction,
    top    = setdiff(ORTH_CROSS_COLS, "flag_fade_bottom"),
    bottom = setdiff(ORTH_CROSS_COLS, "flag_fade_top"),
    ORTH_CROSS_COLS
  )
  list(percentile_cols = pct, cross_cols = cross)
}
