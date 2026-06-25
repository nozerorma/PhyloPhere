#!/usr/bin/env Rscript
# =============================================================================
# CAAS Scoring — Position-Level and Gene-Level Composite Scores
# =============================================================================
#
# Integrates outputs from CT_POSTPROC, FADE, RERConverge, and
# CT_ACCUMULATION into a unified scoring framework.
#
# Usage:
#   Rscript scoring_compute.R \
#     --postproc   <filtered_discovery.tsv> \
#     --fade_top      <fade_summary_top.tsv> \
#     --fade_bottom   <fade_summary_bottom.tsv> \
#     --fade_site_top <fade_site_bf_top.tsv> \
#     --fade_site_bot <fade_site_bf_bottom.tsv> \
#     --rer           <rerconverge_summary.tsv> \
#     --accum_dir     <directory_with_accumulation_CSVs> \
#     --top_pct    0.10 \
#     --top5_pct   0.05 \
#     --top1_pct   0.01 \
#     --gene_top_pct  0.10 \
#     --gene_top5_pct 0.05 \
#     --gene_top1_pct 0.01 \
#     --stress        false \
#     --stress_top_n  25 \
#
# Outputs (in working directory):
#   position_scores.tsv           — per Gene×Position scores
#   gene_scores.tsv               — per Gene scores
#   gene_correlations.tsv         — pairwise correlations between gene scores
#   gene_threshold_enrichment.tsv — gene-level enrichment curve (OR + Fisher) across CAAS thresholds × tools
#   pos_threshold_enrichment.tsv  — position-level FADE enrichment curve across CAAS thresholds
#   gene_lists/slice_*.tsv        — 9 ranked gene lists (Global/Top/Bottom × All/Sig/FDR) for STRING
#   position_score_stress_*.tsv   — leave-one-axis-out / PCA stress diagnostics (only when --stress true)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ─── Argument parsing ────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = "NO_FILE") {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx + 1 > length(args)) return(default)
  args[idx + 1]
}

postproc_file        <- parse_arg("--postproc")
fade_top_file        <- parse_arg("--fade_top")
fade_bottom_file     <- parse_arg("--fade_bottom")
fade_site_top_file   <- parse_arg("--fade_site_top")
fade_site_bot_file   <- parse_arg("--fade_site_bot")
rer_file             <- parse_arg("--rer")
accum_dir            <- parse_arg("--accum_dir")
stress_enabled_raw        <- parse_arg("--stress", "false")
stress_top_n              <- as.integer(parse_arg("--stress_top_n", "25"))
top_pct           <- as.numeric(parse_arg("--top_pct",  "0.10"))
top5_pct          <- as.numeric(parse_arg("--top5_pct", "0.05"))
top1_pct          <- as.numeric(parse_arg("--top1_pct", "0.01"))
gene_top_pct      <- as.numeric(parse_arg("--gene_top_pct",  "0.10"))
gene_top5_pct     <- as.numeric(parse_arg("--gene_top5_pct", "0.05"))
gene_top1_pct     <- as.numeric(parse_arg("--gene_top1_pct", "0.01"))
stress_enabled        <- tolower(as.character(stress_enabled_raw)) %in% c("true", "1", "yes")
if (!is.finite(stress_top_n)      || is.na(stress_top_n)      || stress_top_n < 1)  stress_top_n      <- 25
# direction removed: scoring always runs on the full postproc pool.
# Directional characterisation happens post-scoring via change_side column.

file_exists <- function(f) {
  !is.null(f) && f != "" && !grepl("^NO_", basename(f)) && file.exists(f)
}

# ─── Helper: decile score (0–1 scale) ────────────────────────────────────────
# Lower p-values → higher decile → higher score.
# ntile assigns rank 1 to the smallest values, so we invert: score = 1 - (rank-1)/9
decile_score <- function(x) {
  ranks <- dplyr::ntile(x, 10)
  # For p-values, smaller = better → invert so decile 1 (smallest) gets 1.0
  1 - (ranks - 1) / 9
}

# For Bayes Factors, larger = better → direct mapping
decile_score_higher_better <- function(x) {
  ranks <- dplyr::ntile(x, 10)
  (ranks - 1) / 9
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

safe_top_set <- function(df, score_col, mode = "pct", value = 0.10) {
  vals <- df[[score_col]]
  valid <- !is.na(vals)
  if (!any(valid)) return(character())
  df_valid <- df[valid, , drop = FALSE]
  if (mode == "top_n") {
    n_keep <- min(nrow(df_valid), max(1, as.integer(value)))
    return(df_valid %>% arrange(desc(.data[[score_col]])) %>% slice_head(n = n_keep) %>%
             transmute(id = paste(Gene, Position, sep = "::")) %>% pull(id))
  }
  thr <- quantile(df_valid[[score_col]], 1 - value, na.rm = TRUE)
  df_valid %>%
    filter(.data[[score_col]] >= thr) %>%
    transmute(id = paste(Gene, Position, sep = "::")) %>%
    pull(id) %>%
    unique()
}

pairwise_long <- function(df, cols) {
  expand.grid(score_a = cols, score_b = cols, stringsAsFactors = FALSE) %>%
    filter(score_a < score_b) %>%
    rowwise() %>%
    mutate(
      n_positions = sum(complete.cases(df[[score_a]], df[[score_b]])),
      pearson_r   = safe_cor(df[[score_a]], df[[score_b]], method = "pearson"),
      spearman_r  = safe_cor(df[[score_a]], df[[score_b]], method = "spearman")
    ) %>%
    ungroup()
}

make_rank_matrix <- function(df, cols) {
  out <- expand.grid(score_a = cols, score_b = cols, stringsAsFactors = FALSE) %>%
    filter(score_a < score_b) %>%
    rowwise() %>%
    mutate(
      n_positions = sum(complete.cases(df[[score_a]], df[[score_b]])),
      spearman_rank = safe_cor(rank(df[[score_a]], ties.method = "average", na.last = "keep"),
                               rank(df[[score_b]], ties.method = "average", na.last = "keep"),
                               method = "pearson")
    ) %>%
    ungroup()
  out
}

make_top_overlap <- function(df, cols, top_n = 25) {
  thresholds <- tibble(
    threshold = c("top10pct", "top5pct", "top1pct", "topN"),
    mode = c("pct", "pct", "pct", "top_n"),
    value = c(0.10, 0.05, 0.01, as.numeric(top_n))
  )
  pairs <- expand.grid(score_a = cols, score_b = cols, stringsAsFactors = FALSE) %>%
    filter(score_a < score_b)
  rows <- list()
  k <- 1
  for (i in seq_len(nrow(pairs))) {
    for (j in seq_len(nrow(thresholds))) {
      score_a <- pairs$score_a[i]
      score_b <- pairs$score_b[i]
      threshold <- thresholds$threshold[j]
      mode <- thresholds$mode[j]
      value <- thresholds$value[j]
      set_a <- safe_top_set(df, score_a, mode, value)
      set_b <- safe_top_set(df, score_b, mode, value)
      n_a <- length(set_a)
      n_b <- length(set_b)
      n_intersection <- length(intersect(set_a, set_b))
      n_union <- length(union(set_a, set_b))
      rows[[k]] <- tibble(
        score_a = score_a,
        score_b = score_b,
        threshold = threshold,
        mode = mode,
        value = value,
        n_a = n_a,
        n_b = n_b,
        n_intersection = n_intersection,
        n_union = n_union,
        jaccard = ifelse(n_union > 0, n_intersection / n_union, NA_real_),
        overlap_smaller = ifelse(min(n_a, n_b) > 0, n_intersection / min(n_a, n_b), NA_real_)
      )
      k <- k + 1
    }
  }
  bind_rows(rows)
}

score_to_decile <- function(x) {
  decile_score_higher_better(rank(x, ties.method = "average", na.last = "keep"))
}

variable_numeric_cols <- function(df, cols, min_unique = 2) {
  keep <- vapply(cols, function(col) {
    vals <- suppressWarnings(as.numeric(df[[col]]))
    vals <- vals[is.finite(vals)]
    length(vals) >= 2 && stats::sd(vals) > 0 && length(unique(vals)) >= min_unique
  }, logical(1))
  cols[keep]
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring — Compute\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# 1. LOAD POSTPROC DATA (mandatory)
# =============================================================================
stopifnot(file_exists(postproc_file))
cat("Loading postproc:", postproc_file, "\n")
df <- read_tsv(postproc_file, show_col_types = FALSE)
cat(sprintf("  %d rows, %d unique Gene×Position pairs\n",
            nrow(df), n_distinct(paste(df$Gene, df$Position))))

# =============================================================================
# 2. POSITION-LEVEL SCORING
# =============================================================================

# ── 2a. Pre-scoring exclusions + scheme scope ────────────────────────────────
# Equal weight per scoring scheme (1 each; max achievable sum = 5). Schemes NOT
# in this map (e.g. GS0) are not scoring schemes and are dropped here, BEFORE
# deciles, so the scored pool and every downstream weighted-mean display column
# are consistent. (Previously such rows carried an NA weight, which silently
# dropped them from CAAS_score while NA-poisoning all weighted-mean diagnostics.)
scheme_weights_int <- c(US = 0.2, GS4 = 0.2, GS3 = 0.2, GS2 = 0.2, GS1 = 0.2)

df <- df %>%
  mutate(scheme_weight = scheme_weights_int[CAAP_Group]) %>%
  filter(!is.na(scheme_weight))
cat(sprintf("  %d rows across %d scoring schemes after dropping non-scoring schemes\n",
            nrow(df), n_distinct(df$CAAP_Group)))

# ── 2b. Variability and pattern scores ───────────────────────────────────────
# Deciles computed over the scored pool; lower p-values → higher score.
df <- df %>%
  mutate(
    variability_score = decile_score(Pvalue),
    pattern_score     = decile_score(pvalue_boot)
  )

# Ordinal numeric value for a convergence side label.
# no_change / ambiguous → NA (excluded from weighted mean).
side_value <- function(x) {
  case_when(
    x == "convergent"  ~ 1.00,
    x == "codivergent" ~ 0.50,
    x == "divergent"   ~ 0.25,
    TRUE               ~ NA_real_
  )
}

# Numeric coefficient for a parallel label (categorical → ordinal).
parallel_value <- function(p) {
  case_when(
    is.na(p) | p == "" | tolower(p) == "na" ~ NA_real_,
    tolower(p) == "non-parallel" | p == "False" | p == "FALSE" ~ 1.0,
    tolower(p) == "mixed"        ~ 0.75,
    tolower(p) == "parallel"     | p == "True" | p == "TRUE"  ~ 0.5,
    grepl("^parallel", tolower(p)) ~ 0.5,
    TRUE                          ~ NA_real_
  )
}

# TRUE when a change label indicates an assessable directional event.
assessable_change <- function(x) x %in% c("convergent", "codivergent", "divergent")

# Detect pair-indexed MRCA posterior columns dynamically (mrca_1_posterior, mrca_2_posterior, ...)
mrca_posterior_cols <- grep("^mrca_\\d+_posterior$", names(df), value = TRUE)
n_pairs <- length(mrca_posterior_cols)
cat(sprintf("  Detected %d pairs (%s)\n", n_pairs, paste(mrca_posterior_cols, collapse = ", ")))

# ── 2c. ASR score (per-row) ──────────────────────────────────────────────────
# asr_score sources the unified ASR path score computed upstream in
# CT_DISAMBIGUATION (src/convergence/path_scores.py).
stopifnot("asr_path_score" %in% names(df))
cat("  Using upstream asr_path_score (unified ASR/convergence/parallel signal)\n")
df <- df %>%
  mutate(
    asr_score = suppressWarnings(as.numeric(asr_path_score))
  )

# Ensure the diagnostic path-score columns exist (absent in legacy output).
if (!"mrca_diversity" %in% names(df)) df$mrca_diversity <- NA_real_
df$mrca_diversity <- suppressWarnings(as.numeric(df$mrca_diversity))
if (!"derived_agreement" %in% names(df)) df$derived_agreement <- NA_real_
df$derived_agreement <- suppressWarnings(as.numeric(df$derived_agreement))
if (!"conservation_gate" %in% names(df)) df$conservation_gate <- NA_real_
df$conservation_gate <- suppressWarnings(as.numeric(df$conservation_gate))
if (!"count_factor" %in% names(df)) df$count_factor <- NA_real_
df$count_factor <- suppressWarnings(as.numeric(df$count_factor))
# n_changed_pairs and n_changed_sides removed as they are redundant/misleading diagnostic counts.


# ── 2f. Per-row two-way CAAS score and weighted contribution ─────────────────
# The CAAS score is now the product of TWO orthogonal evidence axes, each already
# on a raw [0,1] scale (no deciles):
#   1. phenotype evolution = permulation evidence (1 - pvalue_boot). Low
#      permutation p-value → strong phenotype-partition signal → high evidence.
#   2. position evolution  = asr_path_score (asr_score). The unified ASR signal,
#      now count-aware (replication re-enters via its count_factor).
# The hypergeometric Pvalue is deliberately absent from the score — it became the
# significance gate (gate_all / gate_sig / gate_fdr, section 2h).
# row_caas:  each scheme contributes at most scheme_weight/5 to the position score
#   (US=GS4=GS3=GS2=GS1=0.2 each; maximum achievable sum = 1.0). The /5 reflects
#   the sum of scheme weights, not the component/axis count, so it is unchanged.
# CAAS_score at position level = sum(row_caas) across all detected schemes.
df <- df %>%
  mutate(
    # ── Display-only diagnostics (deciles): NO LONGER feed CAAS_score ─────────
    # Kept so the report's PCA / stress / correlation sections keep working.
    asr_decile         = decile_score_higher_better(asr_score),
    significance_score  = rowMeans(cbind(variability_score, pattern_score), na.rm = TRUE),
    pre_score           = rowMeans(cbind(significance_score, asr_decile), na.rm = TRUE),

    # ── Two-way CAAS score (per row): phenotype × position evolution ──────────
    # phenotype evolution = permulation evidence (1 - pvalue_boot, the bootstrap
    # permutation p-value); position evolution = asr_path_score (asr_score,
    # already 0..1 and count-aware). Both raw [0,1], higher = better — NO deciles.
    # Their product gives the rank. The hypergeometric Pvalue is NOT here: it is
    # demoted to the significance gate (gate_sig / gate_fdr below).
    phen_score = 1 - pmin(pmax(pvalue_boot, 0), 1),
    caas_row   = phen_score * asr_score,
    row_caas   = caas_row * scheme_weight
  )

# ── 2g. Aggregate to Gene×Position ───────────────────────────────────────────
# Sort descending by scheme_weight so first() picks the most-specific scheme
# for display-only columns (Pvalue, change_top, asr_is_conserved, etc.).
# Score columns use sums or weighted means — no first() for any scored quantity.
df <- df %>% arrange(desc(scheme_weight))

pos_scores <- df %>%
  group_by(Gene, Position) %>%
  summarise(
    CAAS_score         = sum(row_caas,         na.rm = TRUE),
    biochem_weight_sum = sum(scheme_weight,     na.rm = TRUE),
    pre_score          = weighted.mean(pre_score,         w = scheme_weight, na.rm = TRUE),
    significance_score = weighted.mean(significance_score, w = scheme_weight, na.rm = TRUE),
    asr_decile         = weighted.mean(asr_decile,        w = scheme_weight, na.rm = TRUE),
    variability_score  = weighted.mean(variability_score, w = scheme_weight, na.rm = TRUE),
    pattern_score      = weighted.mean(pattern_score,     w = scheme_weight, na.rm = TRUE),
    asr_score          = weighted.mean(asr_score,         w = scheme_weight, na.rm = TRUE),
    mrca_diversity     = weighted.mean(mrca_diversity,   w = scheme_weight, na.rm = TRUE),
    derived_agreement  = weighted.mean(derived_agreement, w = scheme_weight, na.rm = TRUE),
    conservation_gate  = weighted.mean(conservation_gate, w = scheme_weight, na.rm = TRUE),
    count_factor       = weighted.mean(count_factor,     w = scheme_weight, na.rm = TRUE),
    phen_score         = weighted.mean(phen_score,       w = scheme_weight, na.rm = TRUE),
    # Display-only: most-specific scheme via first() after desc(scheme_weight) sort
    Pvalue             = first(Pvalue),
    pvalue_boot        = first(pvalue_boot),
    is_conserved_meta  = first(is_conserved_meta),
    conserved_pair     = first(conserved_pair),
    all_mrca_posterior = first(all_mrca_posterior),
    across(all_of(mrca_posterior_cols), \(x) first(x)),
    change_top         = first(change_top),
    change_bottom      = first(change_bottom),
    has_change_top     = any(assessable_change(change_top)),
    has_change_bottom  = any(assessable_change(change_bottom)),
    .groups = "drop"
  ) %>%
  mutate(
    change_side = case_when(
      has_change_top & has_change_bottom ~ "both",
      has_change_top                     ~ "top",
      has_change_bottom                  ~ "bottom",
      TRUE                               ~ "none"
    )
  ) %>%
  select(-has_change_top, -has_change_bottom)

# ── 2h. Three-way significance gate (hypergeometric Pvalue) ──────────────
# The hypergeometric CAAS p-value no longer feeds CAAS_score (that is now the
# permulation × asr_path product). Instead it gates the result into three tiers,
# so callers can keep all positions for ranking but filter to a defensible set:
#   gate_all : every scored position (the full ranked pool)
#   gate_sig : nominally significant         (Pvalue < 0.05)
#   gate_fdr : significant after BH-FDR       (BH-adjusted Pvalue < 0.05)
# FDR is computed over the whole position pool (one hypergeometric test/position).
pos_scores <- pos_scores %>%
  mutate(
    Pvalue_hyp_fdr = p.adjust(Pvalue, method = "BH"),
    gate_all       = TRUE,
    gate_sig       = !is.na(Pvalue) & Pvalue < 0.05,
    gate_fdr       = !is.na(Pvalue_hyp_fdr) & Pvalue_hyp_fdr < 0.05
  )

cat(sprintf("  %d unique positions after aggregation\n", nrow(pos_scores)))
cat(sprintf("  Significance gate: %d significant (p<0.05), %d FDR-significant (BH<0.05)\n",
            sum(pos_scores$gate_sig, na.rm = TRUE),
            sum(pos_scores$gate_fdr, na.rm = TRUE)))

cat(sprintf("\nPosition-level CAAS_score: min=%.3f, median=%.3f, max=%.3f\n",
            min(pos_scores$CAAS_score, na.rm = TRUE),
            median(pos_scores$CAAS_score, na.rm = TRUE),
            max(pos_scores$CAAS_score, na.rm = TRUE)))

# ── 2i. FADE (gene-level — see section 4d) ──────────────────────────────────
# FADE operates at gene level (max Bayes Factor per gene across sites).
# Significance threshold: BF >= 100 (FADE's standard criterion).
# Top and bottom are loaded separately: fade_top tests acceleration in top
# phenotype direction; fade_bottom tests the bottom direction.
has_fade_top    <- file_exists(fade_top_file)
has_fade_bottom <- file_exists(fade_bottom_file)
has_fade        <- has_fade_top || has_fade_bottom

if (has_fade_top) {
  fade_top_df <- read_tsv(fade_top_file, show_col_types = FALSE) %>%
    select(Gene = gene, fade_max_bf_top = max_bf)
  cat(sprintf("  FADE top: %d genes loaded\n", nrow(fade_top_df)))
} else {
  cat("FADE top: not available, skipping\n")
  fade_top_df <- tibble(Gene = character(), fade_max_bf_top = numeric())
}

if (has_fade_bottom) {
  fade_bottom_df <- read_tsv(fade_bottom_file, show_col_types = FALSE) %>%
    select(Gene = gene, fade_max_bf_bottom = max_bf)
  cat(sprintf("  FADE bottom: %d genes loaded\n", nrow(fade_bottom_df)))
} else {
  cat("FADE bottom: not available, skipping\n")
  fade_bottom_df <- tibble(Gene = character(), fade_max_bf_bottom = numeric())
}

# ── 2j. FADE site-level (position-level BF) ─────────────────────────────────
# Per-site max BF from FADE_report. Coordinate: FADE emits 1-based site index;
# CAAS uses 0-based positions — subtract 1 before any join.
has_fade_site_top <- file_exists(fade_site_top_file)
has_fade_site_bot <- file_exists(fade_site_bot_file)

.load_fade_sites <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    rename(Gene = gene, Position = site) %>%
    mutate(Position = Position - 1L) %>%
    select(Gene, Position, max_site_bf, any_of(c("top_target_aa", "fade_biased")))
}

if (has_fade_site_top) {
  fade_site_top_df <- .load_fade_sites(fade_site_top_file)
  cat(sprintf("  FADE site top: %d positions loaded\n", nrow(fade_site_top_df)))
} else {
  fade_site_top_df <- tibble(Gene = character(), Position = integer(), max_site_bf = numeric())
}

if (has_fade_site_bot) {
  fade_site_bot_df <- .load_fade_sites(fade_site_bot_file)
  cat(sprintf("  FADE site bot: %d positions loaded\n", nrow(fade_site_bot_df)))
} else {
  fade_site_bot_df <- tibble(Gene = character(), Position = integer(), max_site_bf = numeric())
}

# Direction filter removed: all positions are retained for scoring.
# Directional splits (top/bottom) are applied per-analysis after scoring.

# =============================================================================
# 3. POSITION-LEVEL STRESS TESTS (optional)
# =============================================================================
if (stress_enabled) {
  cat("\n─── Position-level scoring stress test ───────────────────────\n")

  stress_df <- pos_scores %>%
    transmute(
      Gene, Position,
      significance_score,            # axis 1: collapsed p-value channels
      asr_decile,                    # axis 2: deciled ASR path score
      variability_score,             # raw inputs kept for diagnostic correlations
      pattern_score,
      asr_score,
      biochem_weight_sum,
      pre_score,
      CAAS_current = CAAS_score
    )

  # Leave-one-out variants: drop one of the TWO scoring axes (significance vs
  # ASR) and rescale. Divisor /5 matches the main row_caas scaling.
  stress_df <- stress_df %>%
    mutate(
      CAAS_unweighted      = pre_score,
      CAAS_no_significance = asr_decile         * biochem_weight_sum / 5,
      CAAS_no_asr          = significance_score * biochem_weight_sum / 5
    )

  latent_loadings <- tibble(model = character(), variable = character(), loading = numeric(),
                            variance_explained = numeric())

  # PCA over the raw inputs (the two p-value deciles + deciled ASR) to show how
  # they load — more informative than a 2-variable PCA on the collapsed axes.
  pca_all_cols <- c("variability_score", "pattern_score", "asr_decile")
  pca_all_df <- stress_df %>% select(all_of(pca_all_cols), CAAS_current) %>% drop_na()
  pca_all_cols_var <- variable_numeric_cols(pca_all_df, pca_all_cols)
  if (length(pca_all_cols_var) < length(pca_all_cols)) {
    cat(sprintf("  PCA_all_components: dropped %d constant/degenerate columns\n",
                length(pca_all_cols) - length(pca_all_cols_var)))
  }
  if (nrow(pca_all_df) >= 3 && length(pca_all_cols_var) >= 2) {
    pca_all <- prcomp(pca_all_df[, pca_all_cols_var, drop = FALSE], center = TRUE, scale. = TRUE)
    scores <- pca_all$x[, 1]
    score_corr <- safe_cor(scores, pca_all_df$CAAS_current, method = "pearson")
    if (!is.na(score_corr) && score_corr < 0) scores <- -scores
    stress_df$PCA_all_components <- NA_real_
    stress_df$PCA_all_components[complete.cases(stress_df[, pca_all_cols_var])] <- score_to_decile(scores)
    latent_loadings <- bind_rows(
      latent_loadings,
      tibble(
        model = "PCA_all_components",
        variable = rownames(pca_all$rotation),
        loading = pca_all$rotation[, 1],
        variance_explained = summary(pca_all)$importance[2, 1]
      )
    )
  } else {
    cat("  PCA_all_components: skipped (fewer than 2 variable columns after filtering)\n")
  }

  factanal_df <- stress_df %>% select(all_of(pca_all_cols), CAAS_current) %>% drop_na()
  factanal_cols_var <- variable_numeric_cols(factanal_df, pca_all_cols)
  if (length(factanal_cols_var) < length(pca_all_cols)) {
    cat(sprintf("  Factor1_all_components: dropped %d constant/degenerate columns\n",
                length(pca_all_cols) - length(factanal_cols_var)))
  }
  if (nrow(factanal_df) > (length(factanal_cols_var) + 2) && length(factanal_cols_var) >= 3) {
    fa_fit <- try(stats::factanal(factanal_df[, factanal_cols_var, drop = FALSE], factors = 1, scores = "regression"), silent = TRUE)
    if (!inherits(fa_fit, "try-error")) {
      scores <- fa_fit$scores[, 1]
      score_corr <- safe_cor(scores, factanal_df$CAAS_current, method = "pearson")
      if (!is.na(score_corr) && score_corr < 0) scores <- -scores
      stress_df$Factor1_all_components <- NA_real_
      stress_df$Factor1_all_components[complete.cases(stress_df[, factanal_cols_var])] <- score_to_decile(scores)
      latent_loadings <- bind_rows(
        latent_loadings,
        tibble(
          model = "Factor1_all_components",
          variable = rownames(fa_fit$loadings[, , drop = FALSE]),
          loading = as.numeric(fa_fit$loadings[, 1]),
          variance_explained = NA_real_
        )
      )
    } else {
      cat("  Factor1_all_components: skipped (factor analysis fit failed)\n")
    }
  } else {
    cat("  Factor1_all_components: skipped (insufficient rows or variable columns)\n")
  }

  analysis_cols <- c(
    "significance_score", "asr_decile",
    "variability_score", "pattern_score", "asr_score",
    "biochem_weight_sum", "pre_score"
  )

  composite_cols <- c(
    "CAAS_current", "CAAS_unweighted",
    "CAAS_no_significance", "CAAS_no_asr"
  )
  latent_cols <- c("PCA_all_components", "Factor1_all_components")
  latent_cols <- latent_cols[latent_cols %in% names(stress_df)]
  latent_cols <- latent_cols[vapply(latent_cols, function(col) any(!is.na(stress_df[[col]])), logical(1))]
  if (length(latent_cols) > 0) composite_cols <- c(composite_cols, latent_cols)

  stress_correlations <- pairwise_long(stress_df, c(analysis_cols, composite_cols))
  stress_rank_agreement <- make_rank_matrix(stress_df, composite_cols)
  stress_top_overlap <- make_top_overlap(stress_df, composite_cols, top_n = stress_top_n)

  # Rank of the reference (current) score is constant across variants — compute
  # once here rather than re-ranking it inside every rowwise cell below.
  rank_current <- rank(-stress_df$CAAS_current, ties.method = "average", na.last = "keep")

  stress_summary <- tibble(variant = composite_cols) %>%
    rowwise() %>%
    mutate(
      pearson_to_current = safe_cor(stress_df[[variant]], stress_df$CAAS_current, method = "pearson"),
      spearman_to_current = safe_cor(stress_df[[variant]], stress_df$CAAS_current, method = "spearman"),
      mean_abs_rank_shift = {
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        mean(abs(rank_current - r1), na.rm = TRUE)
      },
      median_abs_rank_shift = {
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        median(abs(rank_current - r1), na.rm = TRUE)
      },
      max_abs_rank_shift = {
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        max(abs(rank_current - r1), na.rm = TRUE)
      },
      top10_overlap = {
        a <- safe_top_set(stress_df, "CAAS_current", "pct", 0.10)
        b <- safe_top_set(stress_df, variant, "pct", 0.10)
        if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
      },
      top5_overlap = {
        a <- safe_top_set(stress_df, "CAAS_current", "pct", 0.05)
        b <- safe_top_set(stress_df, variant, "pct", 0.05)
        if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
      },
      top1_overlap = {
        a <- safe_top_set(stress_df, "CAAS_current", "pct", 0.01)
        b <- safe_top_set(stress_df, variant, "pct", 0.01)
        if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
      },
      topN_overlap = {
        a <- safe_top_set(stress_df, "CAAS_current", "top_n", stress_top_n)
        b <- safe_top_set(stress_df, variant, "top_n", stress_top_n)
        if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
      }
    ) %>%
    ungroup()

  for (col in composite_cols) {
    stress_df[[paste0("rank_", col)]] <- rank(-stress_df[[col]], ties.method = "average", na.last = "keep")
  }
  rank_cols <- paste0("rank_", composite_cols)
  stress_df$rank_spread <- apply(stress_df[, rank_cols, drop = FALSE], 1, function(x) {
    x <- as.numeric(x[is.finite(x)])
    if (length(x) == 0) NA_real_ else max(x) - min(x)
  })

  write_tsv(stress_summary, "position_score_stress_summary.tsv")
  write_tsv(stress_correlations, "position_score_stress_correlations.tsv")
  write_tsv(stress_rank_agreement, "position_score_stress_rank_agreement.tsv")
  write_tsv(stress_top_overlap, "position_score_stress_top_overlap.tsv")
  write_tsv(stress_df, "position_score_stress_variants.tsv")
  if (nrow(latent_loadings) > 0) {
    write_tsv(latent_loadings, "position_score_stress_latent_loadings.tsv")
  }

  cat(sprintf("  Stress variants: %d columns, %d composite variants\n", ncol(stress_df), length(composite_cols)))
}

# =============================================================================
# 4. GENE-LEVEL SCORING
# =============================================================================

cat("\n─── Gene-level scoring ────────────────────────────────────────\n")

# ── 4a. Gene CAAS Scores: 90th percentile of CAAS_score per gene ────────────
# Three scores computed from different position subsets:
#   gene_caas_score        — all positions (full pool)
#   gene_caas_score_top    — positions with change_side ∈ {top, both}
#   gene_caas_score_bottom — positions with change_side ∈ {bottom, both}
gene_caas <- pos_scores %>%
  group_by(Gene) %>%
  summarise(
    gene_caas_score_global_all = quantile(CAAS_score, 0.90, na.rm = TRUE),
    gene_caas_score_global_sig = {
      vals <- CAAS_score[gate_sig %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_global_fdr = {
      vals <- CAAS_score[gate_fdr %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_top_all    = {
      vals <- CAAS_score[change_side %in% c("top", "both")]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_top_sig    = {
      vals <- CAAS_score[change_side %in% c("top", "both") & gate_sig %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_top_fdr    = {
      vals <- CAAS_score[change_side %in% c("top", "both") & gate_fdr %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_bottom_all = {
      vals <- CAAS_score[change_side %in% c("bottom", "both")]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_bottom_sig = {
      vals <- CAAS_score[change_side %in% c("bottom", "both") & gate_sig %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    gene_caas_score_bottom_fdr = {
      vals <- CAAS_score[change_side %in% c("bottom", "both") & gate_fdr %in% TRUE]
      if (length(vals) > 0) quantile(vals, 0.90, na.rm = TRUE) else NA_real_
    },
    n_positions        = n(),
    n_positions_top    = sum(change_side %in% c("top",    "both"), na.rm = TRUE),
    n_positions_bottom = sum(change_side %in% c("bottom", "both"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    gene_caas_score        = gene_caas_score_global_all,
    gene_caas_score_top    = gene_caas_score_top_all,
    gene_caas_score_bottom = gene_caas_score_bottom_all
  )

cat(sprintf("  gene_caas_score: %d genes (%d with top positions, %d with bottom)\n",
            nrow(gene_caas),
            sum(!is.na(gene_caas$gene_caas_score_top)),
            sum(!is.na(gene_caas$gene_caas_score_bottom))))

# ── 4b. Gene Accumulation Score (optional) ──────────────────────────────────
# Uses accumulation_all_* files: all non-none positions (change_side != 'none').
# Score: accumulation_score = 1 - mean(w_i * p_i) / n_available_schemes
#   where w_i are scheme weights and p_i are PValueEmpirical per scheme.
# Significance: top 5% of accumulation_score (genome-wide empirical threshold).
has_accum <- file_exists(accum_dir) && dir.exists(accum_dir)
if (has_accum) {
  # Check if any accumulation_all_* files exist; skip accumulation if not
  all_accum_files <- list.files(accum_dir, pattern = "^accumulation_all_", full.names = TRUE)
  if (length(all_accum_files) == 0) {
    cat("Accumulation: no accumulation_all_* files found (full pool), skipping\n")
    has_accum <- FALSE
  }
}

if (has_accum) {
  cat("Loading accumulation (all positions) from:", accum_dir, "\n")

  scheme_names <- c("us", "gs4", "gs3", "gs2", "gs1")
  scheme_w     <- c(0.2,  0.2,   0.2,   0.2,   0.2)
  names(scheme_w) <- scheme_names

  accum_pval_df <- NULL  # will hold Gene + one pval col per scheme

  for (scheme in scheme_names) {
    pattern <- paste0("accumulation_all_", scheme, "_aggregated_results.csv")
    f <- list.files(accum_dir, pattern = pattern, full.names = TRUE)
    if (length(f) == 0) {
      cat(sprintf("    %s: file not found, skipping\n", scheme))
      next
    }
    cat(sprintf("    %s: %s\n", scheme, basename(f[1])))
    d <- read_csv(f[1], show_col_types = FALSE)

    pval_col <- grep("PValueEmpirical", names(d), value = TRUE)[1]
    if (is.na(pval_col)) {
      cat(sprintf("    %s: no PValueEmpirical column found, skipping\n", scheme))
      next
    }

    scheme_pvals <- d %>%
      select(Gene, !!paste0("accum_pval_", scheme) := all_of(pval_col))

    if (is.null(accum_pval_df)) {
      accum_pval_df <- scheme_pvals
    } else {
      accum_pval_df <- accum_pval_df %>% full_join(scheme_pvals, by = "Gene")
    }
  }

  if (!is.null(accum_pval_df)) {
    available_schemes <- intersect(scheme_names,
                                   sub("^accum_pval_", "", grep("^accum_pval_", names(accum_pval_df), value = TRUE)))

    # Fisher's combined p-value across available per-group schemes.
    # Each group draws independently, so the per-group p-values are uncorrelated
    # and Fisher's χ²= -2Σln(p) ~ χ²(2k) is valid.
    # Groups absent for a gene contribute p=1 → ln(1)=0, adding no spurious signal.
    gene_rand <- accum_pval_df %>%
      rowwise() %>%
      mutate(
        accum_fisher_p = {
          pvals <- c_across(starts_with("accum_pval_"))
          valid <- !is.na(pvals)
          if (sum(valid) == 0) NA_real_
          else {
            ps   <- pmax(pvals[valid], 1e-300)
            pchisq(-2 * sum(log(ps)), df = 2 * sum(valid), lower.tail = FALSE)
          }
        },
        gene_rand_score = if (is.na(accum_fisher_p)) NA_real_ else -log10(pmax(accum_fisher_p, 1e-300))
      ) %>%
      ungroup() %>%
      select(Gene, gene_rand_score, accum_fisher_p, starts_with("accum_pval_"))

    # BH FDR on genes with at least one observed CAAS (p < 1 strictly).
    # Genes with no CAAS in any group have accum_fisher_p = 1 by construction —
    # they are background members only and must not enter the FDR denominator.
    tested <- !is.na(gene_rand$accum_fisher_p) & gene_rand$accum_fisher_p < 1
    fdr_q  <- rep(NA_real_, nrow(gene_rand))
    if (any(tested))
      fdr_q[tested] <- p.adjust(gene_rand$accum_fisher_p[tested], method = "BH")
    gene_rand$accum_fdr <- fdr_q
    gene_rand <- gene_rand %>%
      mutate(accum_significant = !is.na(accum_fdr) & accum_fdr < 0.05)

    cat(sprintf("  Accumulation: %d genes, %d significant (Fisher p, BH FDR < 0.05)\n",
                nrow(gene_rand),
                sum(gene_rand$accum_significant, na.rm = TRUE)))
  } else {
    gene_rand <- tibble(Gene = character(), gene_rand_score = numeric(),
                        accum_significant = logical())
    has_accum <- FALSE
  }
} else {
  cat("Accumulation: not available, skipping\n")
  gene_rand <- tibble(Gene = character(), gene_rand_score = numeric(),
                      accum_significant = logical())
}

# ── 4c. Gene RERConverge Score (optional) ───────────────────────────────────
has_rer <- file_exists(rer_file)
if (has_rer) {
  cat("Loading RER summary:", rer_file, "\n")
  rer <- read_tsv(rer_file, show_col_types = FALSE)

  # Use p.perm if available, otherwise p.adj
  pval_col <- if ("p.perm" %in% names(rer) && any(!is.na(rer$p.perm) & rer$p.perm > 0)) {
    "p.perm"
  } else {
    "p.adj"
  }
  cat(sprintf("  Using %s for rer_min_pval\n", pval_col))

  gene_rer <- rer %>%
    filter(!is.na(.data[[pval_col]])) %>%
    mutate(rer_min_pval     = .data[[pval_col]],
           rer_significant  = .data[[pval_col]] <= 0.05,
           rer_rho          = Rho,
           rer_acceleration = case_when(
             is.na(Rho) ~ NA_character_,
             Rho > 0    ~ "accelerated",
             Rho < 0    ~ "decelerated",
             TRUE       ~ "neutral"
           )) %>%
    select(Gene = gene, rer_min_pval, rer_significant, rer_rho, rer_acceleration)

  cat(sprintf("  RERConverge: %d genes, %d significant (p <= 0.05)\n",
              nrow(gene_rer), sum(gene_rer$rer_significant, na.rm = TRUE)))
} else {
  cat("RERConverge: not available, skipping\n")
  gene_rer <- tibble(Gene = character(), rer_min_pval = numeric(),
                     rer_significant = logical(), rer_rho = numeric(),
                     rer_acceleration = character())
}

# ── 4d. Gene FADE Significance (optional) ───────────────────────────────────
# BF >= 100 per direction. fade_significant_top / _bottom used separately
# to characterise top-CAAS and bottom-CAAS gene sets respectively.
gene_fade <- tibble(Gene = character())

if (nrow(fade_top_df) > 0) {
  gene_fade <- gene_fade %>%
    full_join(fade_top_df %>% mutate(fade_significant_top = fade_max_bf_top >= 100),
              by = "Gene")
  cat(sprintf("  FADE top: %d significant (BF >= 100)\n",
              sum(gene_fade$fade_significant_top, na.rm = TRUE)))
} else {
  gene_fade$fade_max_bf_top      <- NA_real_
  gene_fade$fade_significant_top <- NA
}

if (nrow(fade_bottom_df) > 0) {
  gene_fade <- gene_fade %>%
    full_join(fade_bottom_df %>% mutate(fade_significant_bottom = fade_max_bf_bottom >= 100),
              by = "Gene")
  cat(sprintf("  FADE bottom: %d significant (BF >= 100)\n",
              sum(gene_fade$fade_significant_bottom, na.rm = TRUE)))
} else {
  gene_fade$fade_max_bf_bottom      <- NA_real_
  gene_fade$fade_significant_bottom <- NA
}

# ── 4e. Assemble gene scores ────────────────────────────────────────────────
gene_scores <- gene_caas

if (nrow(gene_rand) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_rand, by = "Gene")
} else {
  gene_scores$gene_rand_score   <- NA_real_
  gene_scores$accum_significant <- NA
}

if (nrow(gene_rer) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_rer, by = "Gene")
} else {
  gene_scores$rer_min_pval      <- NA_real_
  gene_scores$rer_significant   <- NA
  gene_scores$rer_rho           <- NA_real_
  gene_scores$rer_acceleration  <- NA_character_
}

if (nrow(gene_fade) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_fade, by = "Gene")
} else {
  gene_scores$fade_max_bf_top         <- NA_real_
  gene_scores$fade_significant_top    <- NA
  gene_scores$fade_max_bf_bottom      <- NA_real_
  gene_scores$fade_significant_bottom <- NA
}

# =============================================================================
# 5. CORRELATION ANALYSIS (gene-level)
# =============================================================================

cat("\n─── Correlation analysis ──────────────────────────────────────\n")

# Correlations use gene_caas_score as primary axis; gene_rand_score (accumulation)
# is included as a secondary numeric score. RER and FADE use their
# native significance criteria and are not included as numeric score axes here.
score_cols <- c("gene_caas_score")
if (has_accum && any(!is.na(gene_scores$gene_rand_score)))
  score_cols <- c(score_cols, "gene_rand_score")

if (length(score_cols) >= 2) {
  corr_mat <- gene_scores %>%
    select(all_of(score_cols)) %>%
    drop_na()

  if (nrow(corr_mat) < 3) {
    corr_results <- tibble(
      score_a = character(), score_b = character(),
      pearson_r = numeric(), spearman_r = numeric(), n_genes = integer()
    )
    cat("  Not enough genes with all scores populated, skipping correlations\n")
  } else {

  corr_results <- expand.grid(
    score_a = score_cols, score_b = score_cols,
    stringsAsFactors = FALSE
  ) %>%
    filter(score_a < score_b) %>%
    rowwise() %>%
    mutate(
      pearson_r  = safe_cor(corr_mat[[score_a]], corr_mat[[score_b]], method = "pearson"),
      spearman_r = safe_cor(corr_mat[[score_a]], corr_mat[[score_b]], method = "spearman"),
      n_genes    = sum(complete.cases(corr_mat[[score_a]], corr_mat[[score_b]]))
    ) %>%
    ungroup()
  }

  cat("  Pairwise correlations:\n")
  for (i in seq_len(nrow(corr_results))) {
    cat(sprintf("    %s vs %s: Pearson=%.3f, Spearman=%.3f (n=%d)\n",
                corr_results$score_a[i], corr_results$score_b[i],
                corr_results$pearson_r[i], corr_results$spearman_r[i],
                corr_results$n_genes[i]))
  }
} else {
  corr_results <- tibble(
    score_a = character(), score_b = character(),
    pearson_r = numeric(), spearman_r = numeric(), n_genes = integer()
  )
  cat("  Only one gene score available, skipping correlations\n")
}

# =============================================================================
# 6. WRITE OUTPUTS
# =============================================================================

cat("\n─── Writing outputs ───────────────────────────────────────────\n")

# Position scores — direction-specific top-N% flags (thresholds computed within each direction)
.pos_top_sub <- pos_scores %>% filter(change_side %in% c("top",    "both"))
.pos_bot_sub <- pos_scores %>% filter(change_side %in% c("bottom", "both"))
.thr_top_10 <- quantile(.pos_top_sub$CAAS_score, 0.90, na.rm = TRUE)
.thr_top_5  <- quantile(.pos_top_sub$CAAS_score, 0.95, na.rm = TRUE)
.thr_top_1  <- quantile(.pos_top_sub$CAAS_score, 0.99, na.rm = TRUE)
.thr_bot_10 <- quantile(.pos_bot_sub$CAAS_score, 0.90, na.rm = TRUE)
.thr_bot_5  <- quantile(.pos_bot_sub$CAAS_score, 0.95, na.rm = TRUE)
.thr_bot_1  <- quantile(.pos_bot_sub$CAAS_score, 0.99, na.rm = TRUE)

pos_out <- pos_scores %>%
  select(Gene, Position, any_of(c("Pvalue", "pvalue_boot", "Pvalue_hyp_fdr")),
         asr_score, any_of(c("mrca_diversity", "derived_agreement", "conservation_gate", "count_factor")),
         any_of("phen_score"), biochem_weight_sum, CAAS_score,
         any_of(c("gate_all", "gate_sig", "gate_fdr")), change_side,
         any_of(c("caas", "change_top", "change_bottom"))) %>%
  arrange(desc(CAAS_score))

write_tsv(pos_out, "position_scores.tsv")
cat(sprintf("  position_scores.tsv: %d rows\n", nrow(pos_out)))

# Gene scores
gene_out <- gene_scores %>%
  select(
    Gene,
    n_positions, n_positions_top, n_positions_bottom,
    gene_caas_score_global_all, gene_caas_score_global_sig, gene_caas_score_global_fdr,
    gene_caas_score_top_all, gene_caas_score_top_sig, gene_caas_score_top_fdr,
    gene_caas_score_bottom_all, gene_caas_score_bottom_sig, gene_caas_score_bottom_fdr,
    gene_caas_score, gene_caas_score_top, gene_caas_score_bottom,
    any_of(c("gene_rand_score", "accum_significant",
             "accum_pval_us", "accum_pval_gs4", "accum_pval_gs3",
             "accum_pval_gs2", "accum_pval_gs1")),
    any_of(c("rer_min_pval", "rer_significant", "rer_rho", "rer_acceleration")),
    any_of(c("fade_max_bf_top", "fade_significant_top",
             "fade_max_bf_bottom", "fade_significant_bottom"))
  ) %>%
  arrange(desc(gene_caas_score))

write_tsv(gene_out, "gene_scores.tsv")
cat(sprintf("  gene_scores.tsv: %d rows\n", nrow(gene_out)))

# ── FCS stats table (consumed by FCS_general.Rmd) ─────────────────────────────
# Generic contract: gene + score_<ranking> (zero-floored downstream over the
# cleaned_background universe) + flag_<name> per-gene highlight booleans.
#   score_global/top/bottom = magnitude (the *_all columns); direction via the
#   {top,both}/{bottom,both} membership already baked into the *_top/_bottom cols.
#   flags ride along as leading-edge annotation (NEVER gate the FCS input).
.has  <- function(d, c) c %in% names(d)
.col  <- function(d, c, default = NA) if (.has(d, c)) d[[c]] else rep(default, nrow(d))
.istrue <- function(x) !is.na(x) & x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
.rer_dir <- tolower(as.character(.col(gene_scores, "rer_acceleration")))

fcs_stats <- tibble(
  gene             = gene_scores$Gene,
  score_global     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_global_all"))),
  score_top        = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_top_all"))),
  score_bottom     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_bottom_all"))),
  flag_gate_sig    = !is.na(.col(gene_scores, "gene_caas_score_global_sig")),
  flag_gate_fdr    = !is.na(.col(gene_scores, "gene_caas_score_global_fdr")),
  flag_fade_top    = .istrue(.col(gene_scores, "fade_significant_top")),
  flag_fade_bottom = .istrue(.col(gene_scores, "fade_significant_bottom")),
  flag_rer_acc     = .istrue(.col(gene_scores, "rer_significant")) & grepl("acc", .rer_dir),
  flag_rer_decc    = .istrue(.col(gene_scores, "rer_significant")) & grepl("dec", .rer_dir),
  flag_accum       = .istrue(.col(gene_scores, "accum_significant"))
) %>%
  mutate(flag_fade = flag_fade_top | flag_fade_bottom)
write_tsv(fcs_stats, "fcs_stats.tsv")
cat(sprintf("  fcs_stats.tsv: %d genes (%d top, %d bottom, %d gate_sig)\n",
            nrow(fcs_stats), sum(!is.na(fcs_stats$score_top)),
            sum(!is.na(fcs_stats$score_bottom)), sum(fcs_stats$flag_gate_sig)))

# Correlations
write_tsv(corr_results, "gene_correlations.tsv")

# =============================================================================
# 7. EXPORT RASED/SIGNIFICANCE-GATED GENE LISTS
# =============================================================================

dir.create("gene_lists", showWarnings = FALSE)

# Define the 9 slices
slices_def <- list(
  list(name = "global_all", col = "gene_caas_score_global_all", fade_sig = NA),
  list(name = "global_sig", col = "gene_caas_score_global_sig", fade_sig = NA),
  list(name = "global_fdr", col = "gene_caas_score_global_fdr", fade_sig = NA),

  list(name = "top_all",    col = "gene_caas_score_top_all",    fade_sig = "fade_significant_top"),
  list(name = "top_sig",    col = "gene_caas_score_top_sig",    fade_sig = "fade_significant_top"),
  list(name = "top_fdr",    col = "gene_caas_score_top_fdr",    fade_sig = "fade_significant_top"),

  list(name = "bottom_all", col = "gene_caas_score_bottom_all", fade_sig = "fade_significant_bottom"),
  list(name = "bottom_sig", col = "gene_caas_score_bottom_sig", fade_sig = "fade_significant_bottom"),
  list(name = "bottom_fdr", col = "gene_caas_score_bottom_fdr", fade_sig = "fade_significant_bottom")
)

cat("\n─── Exporting 9 STRING ranked lists ────────────────────────\n")
for (slice in slices_def) {
  col_name <- slice$col
  file_name <- sprintf("gene_lists/slice_%s.tsv", slice$name)

  # Select non-NA genes for this slice
  slice_df <- gene_out %>%
    filter(!is.na(.data[[col_name]])) %>%
    select(Gene, score = all_of(col_name))

  if (nrow(slice_df) == 0) {
    # Write empty file to prevent Nextflow missing output crashes
    write_tsv(tibble(Gene = character(), score = numeric(), is_fade = logical(), is_rer = logical(), is_accum = logical()), file_name)
    cat(sprintf("  %s: empty (0 genes) exported\n", file_name))
    next
  }

  # Add validation columns: is_fade, is_rer, is_accum
  f_col <- slice$fade_sig
  slice_df <- slice_df %>%
    mutate(
      is_fade = if (!is.na(f_col) && f_col %in% names(gene_out)) {
        gene_out[[f_col]][match(Gene, gene_out$Gene)] %in% TRUE
      } else FALSE,
      is_rer = if ("rer_significant" %in% names(gene_out)) {
        is_sig <- gene_out$rer_significant[match(Gene, gene_out$Gene)] %in% TRUE
        rho_val <- gene_out$rer_rho[match(Gene, gene_out$Gene)]
        if (slice$name %in% c("top_all", "top_sig", "top_fdr")) {
          is_sig & !is.na(rho_val) & rho_val > 0
        } else if (slice$name %in% c("bottom_all", "bottom_sig", "bottom_fdr")) {
          is_sig & !is.na(rho_val) & rho_val < 0
        } else {
          is_sig
        }
      } else FALSE,
      is_accum = if ("accum_significant" %in% names(gene_out)) {
        gene_out$accum_significant[match(Gene, gene_out$Gene)] %in% TRUE
      } else FALSE
    ) %>%
    arrange(desc(score))

  write_tsv(slice_df, file_name)
  cat(sprintf("  %s: %d genes exported\n", file_name, nrow(slice_df)))
}

# =============================================================================
# 8. THRESHOLD ENRICHMENT ANALYSIS
# =============================================================================
# For each progressively tighter CAAS threshold, test whether the retained
# set is enriched in independent convergent signals (RER, FADE, accumulation).
# Produces enrichment curves: odds ratio + Fisher p per (direction × threshold × tool).

.fisher_enrich <- function(n_both, n_caas, n_tool, n_total) {
  a  <- n_both
  b  <- n_caas  - n_both
  c_ <- n_tool  - n_both
  d  <- n_total - n_caas - n_tool + n_both
  if (any(c(a, b, c_, d) < 0))
    return(list(or = NA_real_, p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  ft <- tryCatch(
    fisher.test(matrix(c(a, b, c_, d), nrow = 2)),
    error = function(e) list(estimate = NA_real_, p.value = NA_real_,
                              conf.int = c(NA_real_, NA_real_))
  )
  list(or = as.numeric(ft$estimate), p = ft$p.value,
       ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2])
}

enrich_thresholds <- tibble(
  label = c("top100", "top50", "top25", "top10", "top5", "top1"),
  q     = c(0.00,     0.50,   0.75,    0.90,    0.95,   0.99)
)

# ── 8a. Gene-level enrichment ─────────────────────────────────────────────────
cat("\n─── Gene-level threshold enrichment ──────────────────────────\n")
n_total_genes    <- n_distinct(gene_out$Gene)
gene_enrich_rows <- list()
k_ge             <- 1

for (dir in c("top", "bottom", "global")) {
  score_col <- if (dir == "global") "gene_caas_score" else paste0("gene_caas_score_", dir)
  vals  <- gene_out[[score_col]]
  valid <- !is.na(vals)
  if (sum(valid) < 2) next

  tool_sigs <- list()
  if (dir %in% c("top", "global") && has_fade_top &&
      any(gene_out$fade_significant_top %in% TRUE))
    tool_sigs[["fade_top"]] <- gene_out %>%
      filter(fade_significant_top == TRUE) %>% pull(Gene)
  if (dir %in% c("bottom", "global") && has_fade_bottom &&
      any(gene_out$fade_significant_bottom %in% TRUE))
    tool_sigs[["fade_bottom"]] <- gene_out %>%
      filter(fade_significant_bottom == TRUE) %>% pull(Gene)
  if (has_rer && any(gene_out$rer_significant %in% TRUE)) {
    if (dir == "global")
      tool_sigs[["rer_global"]] <- gene_out %>%
        filter(rer_significant == TRUE) %>% pull(Gene)
    if (dir == "top")
      tool_sigs[["rer_accel"]] <- gene_out %>%
        filter(rer_significant == TRUE, !is.na(rer_rho), rer_rho > 0) %>% pull(Gene)
    if (dir == "bottom")
      tool_sigs[["rer_decel"]] <- gene_out %>%
        filter(rer_significant == TRUE, !is.na(rer_rho), rer_rho < 0) %>% pull(Gene)
  }
  if (dir == "global" && has_accum && any(gene_out$accum_significant %in% TRUE))
    tool_sigs[["accum"]] <- gene_out %>%
      filter(accum_significant == TRUE) %>% pull(Gene)

  if (length(tool_sigs) == 0) next

  for (i in seq_len(nrow(enrich_thresholds))) {
    thr_label <- enrich_thresholds$label[i]
    q         <- enrich_thresholds$q[i]
    thr       <- quantile(vals[valid], q, na.rm = TRUE)
    caas_set  <- gene_out %>%
      filter(!is.na(.data[[score_col]]) & .data[[score_col]] >= thr) %>%
      pull(Gene) %>% unique()
    n_caas <- length(caas_set)

    for (tool_name in names(tool_sigs)) {
      tool_set  <- tool_sigs[[tool_name]]
      n_tool    <- length(unique(tool_set))
      n_overlap <- length(intersect(caas_set, tool_set))
      fe        <- .fisher_enrich(n_overlap, n_caas, n_tool, n_total_genes)
      gene_enrich_rows[[k_ge]] <- tibble(
        direction   = dir,
        threshold   = thr_label,
        score_col   = score_col,
        tool        = tool_name,
        n_total     = n_total_genes,
        n_caas_set  = n_caas,
        n_tool_sig  = n_tool,
        n_overlap   = n_overlap,
        enrich_frac = if (n_caas > 0) n_overlap / n_caas else NA_real_,
        odds_ratio  = fe$or,
        or_ci_lo    = fe$ci_lo,
        or_ci_hi    = fe$ci_hi,
        fisher_p    = fe$p
      )
      k_ge <- k_ge + 1
    }
  }
}

if (length(gene_enrich_rows) > 0) {
  gene_threshold_enrichment <- bind_rows(gene_enrich_rows)
  write_tsv(gene_threshold_enrichment, "gene_threshold_enrichment.tsv")
  cat(sprintf("  gene_threshold_enrichment.tsv: %d rows\n", nrow(gene_threshold_enrichment)))
} else {
  cat("  Gene-level enrichment: skipped (no tool data available)\n")
}

# ── 8b. Position-level FADE enrichment ───────────────────────────────────────
has_fade_site_data <- (has_fade_site_top && nrow(fade_site_top_df) > 0) ||
                      (has_fade_site_bot && nrow(fade_site_bot_df) > 0)

if (has_fade_site_data) {
  cat("\n─── Position-level FADE threshold enrichment ─────────────────\n")
  pos_enrich_rows <- list()
  k_pe            <- 1

  for (dir in c("top", "bottom")) {
    pos_subset <- pos_scores %>% filter(change_side %in% c(dir, "both"))
    fade_site  <- if (dir == "top") fade_site_top_df else fade_site_bot_df
    if (nrow(pos_subset) == 0 || nrow(fade_site) == 0) next

    pos_fade <- pos_subset %>%
      left_join(fade_site %>% select(Gene, Position, max_site_bf),
                by = c("Gene", "Position")) %>%
      mutate(fade_sig = !is.na(max_site_bf) & max_site_bf >= 100)

    n_total_pos <- nrow(pos_fade)
    n_fade_sig  <- sum(pos_fade$fade_sig, na.rm = TRUE)
    if (n_fade_sig == 0) {
      cat(sprintf("  FADE site %s: no significant positions (BF >= 100), skipping\n", dir))
      next
    }

    vals  <- pos_fade$CAAS_score
    valid <- !is.na(vals)

    for (i in seq_len(nrow(enrich_thresholds))) {
      thr_label <- enrich_thresholds$label[i]
      q         <- enrich_thresholds$q[i]
      thr       <- quantile(vals[valid], q, na.rm = TRUE)
      caas_pos  <- pos_fade %>% filter(!is.na(CAAS_score) & CAAS_score >= thr)
      n_caas    <- nrow(caas_pos)
      n_overlap <- sum(caas_pos$fade_sig, na.rm = TRUE)
      fe        <- .fisher_enrich(n_overlap, n_caas, n_fade_sig, n_total_pos)
      pos_enrich_rows[[k_pe]] <- tibble(
        direction   = dir,
        threshold   = thr_label,
        n_total     = n_total_pos,
        n_caas_set  = n_caas,
        n_fade_sig  = n_fade_sig,
        n_overlap   = n_overlap,
        enrich_frac = if (n_caas > 0) n_overlap / n_caas else NA_real_,
        odds_ratio  = fe$or,
        or_ci_lo    = fe$ci_lo,
        or_ci_hi    = fe$ci_hi,
        fisher_p    = fe$p
      )
      k_pe <- k_pe + 1
    }
  }

  if (length(pos_enrich_rows) > 0) {
    pos_threshold_enrichment <- bind_rows(pos_enrich_rows)
    write_tsv(pos_threshold_enrichment, "pos_threshold_enrichment.tsv")
    cat(sprintf("  pos_threshold_enrichment.tsv: %d rows\n", nrow(pos_threshold_enrichment)))
  }
} else {
  cat("  Position-level FADE enrichment: skipped (no site-level FADE files)\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring — Complete\n")
cat("═══════════════════════════════════════════════════════════════\n")
