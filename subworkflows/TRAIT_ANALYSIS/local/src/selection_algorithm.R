# selection_algorithm.R — Selection of independent contrast pairs (H1..Hn) for the observed trait.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: 4.Independent_contrasts.Rmd
#
# Observed-trait entry point to the contrast-selection core shared with the
# permulation null (permulations.R). The core,
# subworkflows/CT/local/scripts/lean_contrast_selector.R, is staged into src/
# by ct_independent-contrasts.nf. Defines fop_pair_sel.f().
# =============================================================================

suppressPackageStartupMessages(library(ape))

# Fallback logger, used only when commons.R has not defined debug_log().
if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# ── Shared selection core ─────────────────────────────────────────────────────

# Locate lean_contrast_selector.R: the staged copy in src/ first, then its home
# in the CT scripts directory (relative to the working directory or to this file).
selector_script_path <- {
  this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  this_dir <- if (!is.null(this_ofile)) dirname(this_ofile) else ""
  cand_paths <- c(
    file.path(getwd(), "src", "lean_contrast_selector.R"),
    file.path(getwd(), "..", "..", "CT", "local", "scripts", "lean_contrast_selector.R"),
    file.path(getwd(), "subworkflows", "CT", "local", "scripts", "lean_contrast_selector.R"),
    file.path(this_dir, "lean_contrast_selector.R"),
    file.path(this_dir, "..", "..", "CT", "local", "scripts", "lean_contrast_selector.R")
  )
  found <- cand_paths[nzchar(cand_paths) & file.exists(cand_paths)]
  if (length(found)) found[1] else ""
}

if (nzchar(selector_script_path) && file.exists(selector_script_path)) {
  source(selector_script_path)
} else {
  stop("selection_algorithm.R: could not locate lean_contrast_selector.R ",
       "(the shared contrast-selection core). Looked in: ",
       paste(cand_paths, collapse = " ; "))
}

# ── Multi-hypothesis selection ────────────────────────────────────────────────

# FOP multi-hypothesis contrast selection for the observed trait.
#
# Canonical contrast (H1): the shared candidate gate and ranking
# (lean_candidate_df) followed by greedy Dunn-gated assembly, which stops when no
# candidate keeps every cluster's modified Dunn index >= 1, or at `max_contrasts`.
# H2..Hn: the shared FOP harvest (lean_fop_harvest) around H1, seeded with the
# pipeline seed.
#
# @param ctx           selection_context() for the observed trait.
# @param ci_lb,ci_ub   per-species Jeffreys bounds (count traits), else NULL.
# @param n_vec         per-species sample sizes (count traits), else NULL.
# @param ordinal       TRUE/FALSE ordinal level gate (NULL = auto).
# @param top_pct       continuous-trait PSS gate (params.pss_top_pct).
# @param max_contrasts cap on canonical pairs (Inf = until Dunn stops it).
# @param max_fop       cap on hypotheses (H1 included).
# @param seed          pipeline seed (params.seed).
# @return list(canon_pairs, hypotheses, summary_df, species_domain); summary_df
#   has one row per hypothesis (pair count, min Dunn, mean distance and PSS,
#   Jaccard overlap of its species with H1, pair composition).
fop_pair_sel.f <- function(ctx, ci_lb = NULL, ci_ub = NULL, n_vec = NULL,
                           ordinal = NULL, top_pct, max_contrasts = Inf,
                           max_fop = 100L, seed) {
  empty <- list(
    canon_pairs = data.frame(species1 = character(), species2 = character(),
                             stringsAsFactors = FALSE),
    hypotheses = list(), summary_df = data.frame(), species_domain = integer(0),
    cand_df = data.frame()
  )

  cc <- lean_candidate_df(ctx$trait_vec, ctx$D, 1L, ctx$tree, ctx$cov_bm, ctx$cov_ou,
                          ctx$selected_model, ci_lb, ci_ub, top_pct, ordinal, n_vec)
  if (is.null(cc$cand_df) || nrow(cc$cand_df) == 0) {
    warning("fop_pair_sel.f: no candidate contrast pairs (", cc$reason, "). Returning 0 pairs.")
    return(empty)
  }
  debug_log("fop_pair_sel.f: %d candidate pairs (gate: %s)", nrow(cc$cand_df), cc$mode)

  canon <- greedy_dunn_select(cc$cand_df, ctx$D, target = max_contrasts, enforce_dunn = TRUE)
  canon_pairs <- canon$selected
  K <- nrow(canon_pairs)
  if (K == 0) return(empty)

  hv <- lean_fop_harvest(ctx$trait_vec, ctx$D, K, ctx$tree, ctx$cov_bm, ctx$cov_ou,
                         ctx$selected_model, ci_lb, ci_ub, top_pct, ordinal, n_vec,
                         max_fop = max_fop, seed = seed, canon_pairs = canon_pairs)
  hypotheses <- hv$hypotheses
  debug_log("fop_pair_sel.f: K=%d canonical pairs, %d hypotheses", K, length(hypotheses))

  h1_sp <- c(canon_pairs$species1, canon_pairs$species2)
  .agg <- function(p, f) if (all(is.na(p))) NA_real_ else round(f(p, na.rm = TRUE), 4)
  summary_df <- do.call(rbind, lapply(names(hypotheses), function(h_id) {
    hdf <- hypotheses[[h_id]]
    all_sp <- c(hdf$species1, hdf$species2)
    data.frame(
      hypothesis_id    = h_id,
      is_canonical     = identical(h_id, "H1"),
      num_pairs        = nrow(hdf),
      min_dunn         = round(hv$dunn[[h_id]], 4),
      mean_distance    = round(mean(hdf$distance), 4),
      mean_abs_diff    = round(mean(hdf$abs_diff), 4),
      mean_pss_score   = .agg(hdf$pss_score, mean),
      min_pss_score    = .agg(hdf$pss_score, min),
      jaccard_to_h1    = round(length(intersect(h1_sp, all_sp)) / length(union(h1_sp, all_sp)), 4),
      pair_composition = paste(paste(hdf$species1, hdf$species2, sep = "~"), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))

  list(
    canon_pairs = canon_pairs,
    hypotheses = hypotheses,
    summary_df = summary_df,
    species_domain = hv$species_domain,
    cand_df = cc$cand_df
  )
}
