# Position-level empirical p and BH over the captured positional universe.
# p.emp = (k_emp + 1)/(N + 1), inclusive upper tail. Null-only positions enter BH at p=1.
# Legacy fit/c_k inputs are accepted but ignored; no pooled conditional-tail estimator is produced.

# obs_max   data frame (Gene, Position, .obs): the observed positive positions and their best-side score
# counters  data frame (Gene, Position, k_emp, nd): one row per observed position
# fit       ignored legacy argument, retained for cached callers
# n_u       size of U among the null-positive positions (distinct positions)
# N         cycles of the null (the roster)
# Returns list(pos = data frame (Gene, Position, p.emp, p.adj_bh, nd, null_unseen),
#              n_matched, n_obs, n_unseen, n_family)
null_position_p <- function(obs_max, counters, fit, n_u, N, TIE_TOL = 1e-12, P_EMP_GUARD_MIN_POS = 10L) {
  stopifnot(all(c("Gene", "Position", ".obs") %in% names(obs_max)))
  at <- match(paste(obs_max$Gene, obs_max$Position), paste(counters$Gene, counters$Position))
  if (anyNA(at)) stop("the null counters lack ", sum(is.na(at)), " observed position(s): counters and observed run disagree")
  ctr <- counters[at, , drop = FALSE]
  # An observed position that no null cycle re-detects has k_emp = 0 (its null statistic is -Inf in every cycle).
  # Observed scores are positive here (position_scores.tsv keeps only CAAS_score > 0), so the guard on .obs <= TIE_TOL
  # only protects against a score read as 0 within the tolerance: such a value is no evidence of convergence, p.emp = 1.
  nd <- as.integer(ctr$nd)
  null_hit <- nd > 0L
  p_emp <- ifelse(obs_max$.obs <= TIE_TOL, 1, (ctr$k_emp + 1) / (N + 1))

  n_obs <- nrow(obs_max)
  n_matched <- sum(null_hit)
  rate <- if (n_obs > 0) n_matched / n_obs else 0
  cat(sprintf("  p.emp: %d/%d observed positions re-detected by the null (%.1f%%), N=%d cycles\n",
              n_matched, n_obs, 100 * rate, N))
  # The overlap rate only diagnoses a coordinate mismatch when there are enough observed positions to estimate it. Below
  # P_EMP_GUARD_MIN_POS an unmatched position is taken at face value: no null cycle re-detects it, so k_emp = 0.
  if (n_obs < P_EMP_GUARD_MIN_POS && n_matched < n_obs) {
    cat(sprintf(paste0("  p.emp: %d of %d observed position(s) never re-detected by the null; scored at k_emp = 0 ",
                       "(too few positions for the coordinate-mismatch check)\n"), n_obs - n_matched, n_obs))
  }
  if (n_obs >= P_EMP_GUARD_MIN_POS && rate < 0.5) {
    # Low overlap means the two tables are on different coordinate systems, not that the unmatched positions are strong:
    # leave them untested (NA).
    p_emp[!null_hit] <- NA_real_
    cat(sprintf(paste0("  WARNING: p.emp join rate %.1f%% < 50%% -- the observed positions (filtered_discovery.tsv) and ",
                       "the positions of the null are likely on different coordinate systems. Unmatched positions left ",
                       "NA; treat p.emp/p.adj_bh as unreliable.\n"), 100 * rate), file = stderr())
  }

  # Diagnostics, informational only: nd is the number of cycles of the null that score the position, null_unseen flags
  # nd == 0. Where the coordinate guard left p.emp NA they are NA too.
  nd_out <- nd
  nd_out[is.na(p_emp)] <- NA_integer_
  null_unseen <- nd_out == 0L

  # BH over the family: U plus the observed positions with a p. The observed positions with a null hit are in U already.
  has_p <- !is.na(p_emp)
  n_only <- n_u - sum(has_p & null_hit)
  n_family <- n_u + sum(has_p & !null_hit)
  fam_e <- c(p_emp[has_p], rep(1, n_only))
  adj_e <- rep(NA_real_, n_obs)
  if (length(fam_e) > 0) {
    adj_e[has_p] <- p.adjust(fam_e, method = "BH")[seq_len(sum(has_p))]
  }
  cat(sprintf("  p.adj_bh: BH over %d positions (%d observed, %d null-only at p = 1)\n",
              length(fam_e), sum(has_p), n_only))

  # columns in the order position_scores.tsv carries them
  list(pos = data.frame(Gene = obs_max$Gene, Position = obs_max$Position, p.emp = p_emp,
                        nd = nd_out, null_unseen = null_unseen, p.adj_bh = adj_e,
                        stringsAsFactors = FALSE),
       n_matched = n_matched, n_obs = sum(!is.na(nd_out)), n_unseen = sum(null_unseen, na.rm = TRUE),
       n_family = length(fam_e))
}

# Rows of null_universe_summary.tsv (metric, cycles, value, note): N, |U|, observed positive positions, the ones the null
# never scored, and the saturation curve of the null. `base` holds what the null itself fixes (null_universe_base.tsv of
# null_position_stats.py); an unseen fraction is a diagnostic of the adequacy of the null (does N reach the positions that
# can be scored?), not a reason to demote a position.
null_universe_summary <- function(base, n_obs, n_unseen) {
  val <- function(m) base$value[base$metric == m & is.na(base$cycles)][1]
  note <- function(m) { v <- base$note[base$metric == m & is.na(base$cycles)][1]; if (is.na(v)) "" else v }
  n_u <- val("null_positive_positions_in_universe")
  metrics <- c("N_cycles", "N_cycles_present_in_null", "null_positive_positions_all", "null_positive_positions_in_universe",
               "observed_positive_positions", "observed_unseen_positions", "fraction_unseen", "universe_size")
  curve <- base[base$metric == "null_positive_positions_cumulative", c("metric", "cycles", "value", "note")]
  rbind(
    data.frame(metric = metrics, cycles = NA_integer_,
               value = c(val("N_cycles"), val("N_cycles_present_in_null"), val("null_positive_positions_all"), n_u, n_obs,
                         n_unseen, if (n_obs > 0) n_unseen / n_obs else NA_real_, n_u + n_unseen),
               note = c(note("N_cycles"), "", "", note("null_positive_positions_in_universe"), "",
                        "positions whose null detections are 0 (p = 1/(N+1), kept)", "",
                        "null-positive positions of U plus the observed positions the null never scored"),
               stringsAsFactors = FALSE),
    curve)
}
