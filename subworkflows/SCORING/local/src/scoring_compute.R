#!/usr/bin/env Rscript
# =============================================================================
# CAAS Scoring - Position-Level and Gene-Level Composite Scores
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
#     --gene_top_pct  0.10 \
#     --stress        false \
#     --stress_top_n  25 \
#
# Outputs (in working directory):
#   position_scores.tsv - per Gene×Position scores
#   gene_scores.tsv - per Gene scores
#   gene_correlations.tsv - pairwise correlations between gene scores
#   gene_threshold_enrichment.tsv - gene-level enrichment curve (OR + Fisher) across CAAS thresholds × tools
#   pos_threshold_enrichment.tsv - position-level FADE enrichment curve across CAAS thresholds
#   gene_lists/slice_*.tsv - 8 ranked gene lists (Top/Bottom × 25/10/5/1%) for STRING
#   position_score_stress_*.tsv - leave-one-axis-out / PCA stress diagnostics (only when --stress true)
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
hyp_pairs_file       <- parse_arg("--hypotheses_pairs")  # contrast_hypotheses_pairs.tsv (FOP); NO_HYP_PAIRS otherwise
stress_enabled_raw        <- parse_arg("--stress", "false")
stress_top_n              <- as.integer(parse_arg("--stress_top_n", "25"))
top_pct           <- as.numeric(parse_arg("--top_pct",  "0.10"))
top25_pct         <- 0.25
top5_pct          <- 0.05
top1_pct          <- 0.01
gene_top_pct      <- as.numeric(parse_arg("--gene_top_pct",  "0.10"))
gene_top25_pct    <- 0.25
gene_top5_pct     <- 0.05
gene_top1_pct     <- 0.01
stress_enabled        <- tolower(as.character(stress_enabled_raw)) %in% c("true", "1", "yes")
if (!is.finite(stress_top_n)      || is.na(stress_top_n)      || stress_top_n < 1)  stress_top_n      <- 25
# direction removed: scoring always runs on the full postproc pool.
# Directional characterisation happens post-scoring via change_side column.

file_exists <- function(f) {
  !is.null(f) && f != "" && !grepl("^NO_", basename(f)) && file.exists(f)
}

# ── Helper: size-adjusted maximum (0–1 scale) ────────────────────────────────
# max(x), like any upper quantile, grows with length(x) as a pure order-
# statistic artifact: more draws from the same distribution push the observed
# extremum further into the tail. Two genes with identical per-position evidence
# therefore score differently if one merely has more detected positions
# (measured on real output: Spearman(n_positions, raw max) = +0.47; see
# section 4a for how gene_caas_score uses this).
#
# size_adj_max compares a gene's observed max against the distribution of the
# max of n draws from the genome-wide position pool, n being that gene's own
# position count. For iid draws that reference distribution is exact and needs
# no resampling, since P(max of n <= m) = F(m)^n with F the pool ECDF:
#
#     size_adj_max(x) = F(max(x)) ^ length(x)
#
# Two properties matter downstream:
#   * Monotone in max(x) at fixed n, so genes with the SAME n keep their exact
#     relative order - the transform only makes different-n genes comparable,
#     it does not reshuffle within a size class.
#   * n-neutral: Spearman(n_positions, .) = -0.09 on real output.
# Validated against a 4000-draw resampled reference (Spearman 0.995, top-1%
# overlap 55/58); the closed form is used because it is exact and seed-free.
#
# pool_sorted MUST be sorted ascending - findInterval then counts pool <= m in
# O(log N) rather than O(N) per gene.
size_adj_max <- function(x, pool_sorted) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || length(pool_sorted) == 0) return(NA_real_)
  (findInterval(max(x), pool_sorted) / length(pool_sorted)) ^ length(x)
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

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring - Compute\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# 1. LOAD POSTPROC DATA (mandatory)
# =============================================================================
stopifnot(file_exists(postproc_file))
cat("Loading postproc:", postproc_file, "\n")
df <- read_tsv(postproc_file, show_col_types = FALSE)
# filtered_discovery.tsv uses disambiguation's canonical lowercase concept names
# (caap_group, pvalue, …) - consumed as-is; position_scores.tsv keeps the same lowercase schema.
cat(sprintf("  %d rows, %d unique Gene×Position pairs\n",
            nrow(df), n_distinct(paste(df$Gene, df$Position))))

# =============================================================================
# 2. POSITION-LEVEL SCORING
# =============================================================================

# ── 2a. Pre-scoring exclusions + scheme scope ────────────────────────────────
# The five scoring schemes.
#
# No per-scheme weight: how many of the five schemes detect a substitution is a
# deterministic property of which amino acids are involved (a discretised
# biochemical distance, see the report's Biochemistry tab), not evidence
# strength. Section 2g aggregates schemes with a mean for this reason.
scoring_schemes <- c("US", "GS4", "GS3", "GS2", "GS1")

# Priority ONLY for picking a representative scheme's display/gating columns
# (pvalue, pvalue_boot, change_side, caap_group, ...) at the Gene×Position
# aggregation below (section 2g). Deliberately separate from the scoring itself,
# which treats all five schemes symmetrically.
scheme_priority_int <- c(US = 5, GS4 = 4, GS3 = 3, GS2 = 2, GS1 = 1)

df <- df %>%
  mutate(
    scheme_priority = scheme_priority_int[caap_group],
    # FOP discovering-hypothesis tag ("H<n>") or NA for a single-contrast run
    # (trait == "post_disambiguation" or similar). Only H-tags drive fop_pool.R.
    hyp_id = if ("trait" %in% names(df)) ifelse(grepl("H[0-9]+", trait), sub(".*(H[0-9]+).*", "\\1", trait), NA_character_) else NA_character_
  ) %>%
  filter(caap_group %in% scoring_schemes)
cat(sprintf("  %d rows across %d scoring schemes after dropping non-scoring schemes\n",
            nrow(df), n_distinct(df$caap_group)))

# ── 2b. FOP domain-pooling (collapse H1..Hn -> one row per Gene×Position×scheme) ─
# A FOP run emits one disambiguation row per (Gene, Position, scheme, hypothesis).
# H1..Hn are overlapping K-pair designs over the same Voronoi domains, NOT
# independent replicates, so §2g's per-scheme mean must not also average over
# them uniformly (it would dilute a strong canonical signal and let a position
# with many harvested hypotheses distort every genome-wide rank). fop_pool.R
# pools s(p,site) within each Voronoi domain (PSS-weighted mean, weights from
# contrast_hypotheses_pairs.tsv) and recombines with the path_scores.py algebra.
# Non-FOP input (single contrast) passes through unchanged.
.fop_pool_src <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fop_pool.R")
if (!file.exists(.fop_pool_src)) .fop_pool_src <- "fop_pool.R"
source(.fop_pool_src)
.n_before <- nrow(df)
df <- apply_fop_pooling(df, hyp_pairs_file)
if (nrow(df) != .n_before) {
  cat(sprintf("  FOP pooling: %d rows -> %d after collapsing hypotheses (max n_hypotheses = %d)\n",
              .n_before, nrow(df),
              if ("n_hypotheses" %in% names(df)) max(df$n_hypotheses, na.rm = TRUE) else 0L))
} else {
  cat("  FOP pooling: no multi-hypothesis positions (single-contrast run) — pass-through\n")
}
df$asr_path_score <- suppressWarnings(as.numeric(df$asr_path_score))

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

# Ensure the diagnostic path-score columns exist even if the input file lacks them.
if (!"mrca_diversity" %in% names(df)) df$mrca_diversity <- NA_real_
df$mrca_diversity <- suppressWarnings(as.numeric(df$mrca_diversity))
if (!"derived_agreement" %in% names(df)) df$derived_agreement <- NA_real_
df$derived_agreement <- suppressWarnings(as.numeric(df$derived_agreement))
if (!"conservation_gate" %in% names(df)) df$conservation_gate <- NA_real_
df$conservation_gate <- suppressWarnings(as.numeric(df$conservation_gate))
if (!"core" %in% names(df)) df$core <- NA_real_
df$core <- suppressWarnings(as.numeric(df$core))


# ── 2f. Per-row two-way CAAS score ────────────────────────────────────────────
# caas_row is the product of two orthogonal [0,1] evidence axes, computed per
# (position, scheme) row:
#   phen_score = 1 - percent_rank(pvalue_boot): permutation confidence.
#   asr_score  = asr_path_score: the unified ASR signal (section above).
# The hypergeometric pvalue is not part of this product; it is the significance
# gate (gate_all / gate_sig, section 2h).
# percent_rank is genome-wide, so it must see each (Gene, Position, scheme)
# exactly once. §2b already collapsed the FOP hypothesis rows, so pvalue_boot
# here is the representative (highest-asr) hypothesis's permutation p and the
# rank is undistorted by how many hypotheses a position happened to harvest.
df <- df %>%
  mutate(
    phen_score = 1 - dplyr::percent_rank(pvalue_boot),
    caas_row   = phen_score * asr_score
  )

# ── 2g. Aggregate to Gene×Position ───────────────────────────────────────────
# Sort descending by scheme_priority (US > GS4 > GS3 > GS2 > GS1) so first()
# deterministically picks the US scheme (falling back to GS4..GS1) for
# display/gating-only columns (pvalue, change_top, asr_is_conserved, etc.).
# Priority is display-only and never enters a scored quantity.
df <- df %>% arrange(desc(scheme_priority))

pos_scores <- df %>%
  group_by(Gene, Position) %>%
  summarise(
    CAAS_score         = mean(caas_row, na.rm = TRUE),
    # FOP descriptors: §2b already pooled H1..Hn, so these are per-position
    # columns now, not a re-count over rows. Recurrence stays descriptor-only —
    # it never multiplies CAAS_score.
    n_hypotheses       = if ("n_hypotheses" %in% names(df)) {
      .nh <- n_hypotheses[is.finite(n_hypotheses)]; if (length(.nh)) max(.nh) else 0L
    } else 0L,
    supporting_hypotheses = if ("supporting_hypotheses" %in% names(df)) {
      .sh <- unique(supporting_hypotheses[!is.na(supporting_hypotheses) & nzchar(supporting_hypotheses)])
      if (length(.sh)) paste(sort(unique(unlist(strsplit(.sh, ",")))), collapse = ",") else ""
    } else "",
    n_schemes          = dplyr::n(),
    scheme_set         = paste(sort(unique(as.character(caap_group))), collapse = "+"),
    asr_score          = mean(asr_score,          na.rm = TRUE),
    mrca_diversity     = mean(mrca_diversity,     na.rm = TRUE),
    derived_agreement  = mean(derived_agreement,  na.rm = TRUE),
    conservation_gate  = mean(conservation_gate,  na.rm = TRUE),
    core               = mean(core,               na.rm = TRUE),
    core_perside_pooled = if ("core_perside_pooled" %in% names(df)) mean(core_perside_pooled, na.rm = TRUE) else NA_real_,
    phen_score         = mean(phen_score,         na.rm = TRUE),
    pvalue_boot        = first(pvalue_boot),
    is_conserved_meta  = first(is_conserved_meta),
    conserved_pair     = first(conserved_pair),
    all_mrca_posterior = first(all_mrca_posterior),
    across(all_of(mrca_posterior_cols), \(x) first(x)),
    change_top         = first(change_top),
    change_bottom      = first(change_bottom),
    caap_group         = first(caap_group),
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

cat(sprintf("  %d unique positions after aggregation\n", nrow(pos_scores)))

cat(sprintf("\nPosition-level CAAS_score: min=%.3f, median=%.3f, max=%.3f\n",
            min(pos_scores$CAAS_score, na.rm = TRUE),
            median(pos_scores$CAAS_score, na.rm = TRUE),
            max(pos_scores$CAAS_score, na.rm = TRUE)))

# ── 2i. FADE (gene-level - see section 4d) ──────────────────────────────────
# FADE operates at gene level (max Bayes Factor per gene across sites).
# Significance threshold: BF >= 100 (FADE's standard criterion).
# Top and bottom are loaded separately: fade_top tests acceleration in top
# phenotype direction; fade_bottom tests the bottom direction.
has_fade_top    <- file_exists(fade_top_file)
has_fade_bottom <- file_exists(fade_bottom_file)
has_fade        <- has_fade_top || has_fade_bottom

.load_fade_summary <- function(path, out_col) {
  df <- read_tsv(path, show_col_types = FALSE)
  # Standardise column names to lowercase for checking
  names(df) <- tolower(names(df))
  if (!"gene" %in% names(df)) {
    names(df)[1] <- "gene"
  }
  bf_col <- intersect(c("max_bf", "max_site_bf", "bayes_factor", "bf"), names(df))[1]
  if (is.na(bf_col)) {
    df[[out_col]] <- 0
  } else {
    df[[out_col]] <- as.numeric(df[[bf_col]])
  }
  df %>%
    select(Gene = gene, !!out_col := !!sym(out_col)) %>%
    mutate(Gene = as.character(Gene))
}

if (has_fade_top) {
  fade_top_df <- .load_fade_summary(fade_top_file, "fade_max_bf_top")
  cat(sprintf("  FADE top: %d genes loaded\n", nrow(fade_top_df)))
} else {
  cat("FADE top: not available, skipping\n")
  fade_top_df <- tibble(Gene = character(), fade_max_bf_top = numeric())
}

if (has_fade_bottom) {
  fade_bottom_df <- .load_fade_summary(fade_bottom_file, "fade_max_bf_bottom")
  cat(sprintf("  FADE bottom: %d genes loaded\n", nrow(fade_bottom_df)))
} else {
  cat("FADE bottom: not available, skipping\n")
  fade_bottom_df <- tibble(Gene = character(), fade_max_bf_bottom = numeric())
}

# ── 2j. FADE site-level (position-level BF) ─────────────────────────────────
# Per-site max BF from FADE_report. Coordinate: FADE emits 1-based site index;
# CAAS uses 0-based positions - subtract 1 before any join.
has_fade_site_top <- file_exists(fade_site_top_file)
has_fade_site_bot <- file_exists(fade_site_bot_file)

.load_fade_sites <- function(path) {
  # FADE_JSON_TO_CSV's real (and only) producer of this file
  # (parse_fade_json_sites.R) writes comma-delimited (write.csv()), not tab --
  # read_tsv() here would silently parse the whole line as one column and
  # crash downstream at df[[2]]. read_csv() matches the actual pipeline output.
  df <- read_csv(path, show_col_types = FALSE)
  names(df) <- tolower(names(df))
  
  if (!"gene" %in% names(df)) names(df)[1] <- "gene"
  if ("site" %in% names(df)) {
    df$position <- as.integer(df$site) - 1L
  } else if ("position" %in% names(df)) {
    df$position <- as.integer(df$position)
  } else {
    df$position <- as.integer(df[[2]]) - 1L
  }
  
  bf_col <- intersect(c("max_site_bf", "max_bf", "bayes_factor", "bf"), names(df))[1]
  if (is.na(bf_col)) {
    df$max_site_bf <- 0
  } else {
    df$max_site_bf <- as.numeric(df[[bf_col]])
  }
  
  aa_col <- intersect(c("top_target_aa", "target_aa", "aa"), names(df))[1]
  if (!is.na(aa_col)) {
    df$top_target_aa <- as.character(df[[aa_col]])
  } else {
    df$top_target_aa <- NA_character_
  }
  
  bias_col <- intersect(c("fade_biased", "biased"), names(df))[1]
  if (!is.na(bias_col)) {
    df$fade_biased <- as.logical(df[[bias_col]])
  } else {
    df$fade_biased <- df$max_site_bf >= 100
  }
  
  df %>%
    select(Gene = gene, Position = position, max_site_bf, top_target_aa, fade_biased)
}

if (has_fade_site_top) {
  fade_site_top_df <- .load_fade_sites(fade_site_top_file)
  cat(sprintf("  FADE site top: %d positions loaded\n", nrow(fade_site_top_df)))
} else {
  fade_site_top_df <- tibble(Gene = character(), Position = integer(), max_site_bf = numeric(), top_target_aa = character(), fade_biased = logical())
}

if (has_fade_site_bot) {
  fade_site_bot_df <- .load_fade_sites(fade_site_bot_file)
  cat(sprintf("  FADE site bot: %d positions loaded\n", nrow(fade_site_bot_df)))
} else {
  fade_site_bot_df <- tibble(Gene = character(), Position = integer(), max_site_bf = numeric(), top_target_aa = character(), fade_biased = logical())
}

# Direction filter removed: all positions are retained for scoring.
# Directional splits (top/bottom) are applied per-analysis after scoring.

# =============================================================================
# 3. POSITION-LEVEL STRESS TESTS (optional)
# =============================================================================
if (stress_enabled) {
  cat("\n─── Position-level scoring stress test ───────────────────────\n")

  # Variants of the ACTUAL production inputs: phen_score, asr_score (the two
  # factors of caas_row, section 2f) and the scheme-aggregation rule (section
  # 2g). phen_score/asr_score here are pos_scores' own mean-over-schemes
  # values, identical to what CAAS_score is built from -- no separate decile
  # transform, no third input.
  scheme_max <- df %>%
    group_by(Gene, Position) %>%
    summarise(CAAS_scheme_max = max(caas_row, na.rm = TRUE), .groups = "drop")

  stress_df <- pos_scores %>%
    select(Gene, Position, phen_score, asr_score, n_schemes, CAAS_current = CAAS_score) %>%
    left_join(scheme_max, by = c("Gene", "Position")) %>%
    mutate(
      CAAS_phen_only = phen_score,  # drop the ASR axis
      CAAS_asr_only  = asr_score    # drop the phenotype axis
    )

  analysis_cols <- c("phen_score", "asr_score", "n_schemes")
  composite_cols <- c("CAAS_current", "CAAS_phen_only", "CAAS_asr_only", "CAAS_scheme_max")

  stress_correlations <- pairwise_long(stress_df, c(analysis_cols, composite_cols))
  stress_rank_agreement <- make_rank_matrix(stress_df, composite_cols)
  stress_top_overlap <- make_top_overlap(stress_df, composite_cols, top_n = stress_top_n)

  # Rank of the reference (current) score is constant across variants - compute
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

  cat(sprintf("  Stress variants: %d columns, %d composite variants\n", ncol(stress_df), length(composite_cols)))
}

# =============================================================================
# 4. GENE-LEVEL SCORING
# =============================================================================

cat("\n─── Gene-level scoring ────────────────────────────────────────\n")

# ── 4a. Gene CAAS Scores: size-adjusted max of CAAS_score per gene ──────────
# Three scores computed from different position subsets:
#   gene_caas_score - all positions (full pool)
#   gene_caas_score_top - positions with change_side in {top, both}
#   gene_caas_score_bottom - positions with change_side in {bottom, both}
#
# Aggregation is size_adj_max (helper at the top of this file): a gene's best
# position, calibrated for how many positions the gene had a chance to draw
# from. Any fixed quantile of a gene's positions (including the 90th) grows
# with the gene's own position count at equal per-position evidence, since a
# higher quantile from more draws pushes further into the tail of the same
# distribution; size_adj_max removes that by comparing the observed max against
# the distribution of the max of n draws from the genome-wide pool, n being the
# gene's own position count.
#
# Reference pools are direction-matched: a gene's top-direction positions are
# ranked against the genome-wide pool of top-direction positions, so each score
# is calibrated against the distribution it is actually drawn from.
.pool_all    <- sort(pos_scores$CAAS_score[!is.na(pos_scores$CAAS_score)])
.pool_top    <- sort(pos_scores$CAAS_score[pos_scores$change_side %in% c("top", "both") &
                                           !is.na(pos_scores$CAAS_score)])
.pool_bottom <- sort(pos_scores$CAAS_score[pos_scores$change_side %in% c("bottom", "both") &
                                           !is.na(pos_scores$CAAS_score)])
cat(sprintf("  size-adjust reference pools: all=%d, top=%d, bottom=%d positions\n",
            length(.pool_all), length(.pool_top), length(.pool_bottom)))

gene_caas <- pos_scores %>%
  group_by(Gene) %>%
  summarise(
    gene_caas_score            = size_adj_max(CAAS_score, .pool_all),
    gene_caas_score_top_all    = {
      vals <- CAAS_score[change_side %in% c("top", "both")]
      if (length(vals) > 0) size_adj_max(vals, .pool_top) else NA_real_
    },
    gene_caas_score_bottom_all = {
      vals <- CAAS_score[change_side %in% c("bottom", "both")]
      if (length(vals) > 0) size_adj_max(vals, .pool_bottom) else NA_real_
    },
    n_positions        = n(),
    n_positions_top    = sum(change_side %in% c("top",    "both"), na.rm = TRUE),
    n_positions_bottom = sum(change_side %in% c("bottom", "both"), na.rm = TRUE),
    max_hypotheses     = if ("n_hypotheses" %in% names(pos_scores)) max(n_hypotheses, na.rm = TRUE) else NA_integer_,
    mean_hypotheses    = if ("n_hypotheses" %in% names(pos_scores)) round(mean(n_hypotheses, na.rm = TRUE), 1) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    gene_caas_score_top    = gene_caas_score_top_all,
    gene_caas_score_bottom = gene_caas_score_bottom_all
  )

cat(sprintf("  gene_caas_score: %d genes (%d with top positions, %d with bottom)\n",
            nrow(gene_caas),
            sum(!is.na(gene_caas$gene_caas_score_top)),
            sum(!is.na(gene_caas$gene_caas_score_bottom))))

# ── 4b. Gene Accumulation Score (optional) ──────────────────────────────────
# Uses accumulation_<direction>_* files: direction in {all, top, bottom}.
# "all" = every non-none position pooled (change_side != 'none'); "top"/"bottom"
# restrict to that direction's positions only (change_side %in% c(dir, 'both')),
# matching how every other cross-module flag (FADE, RER) is direction-aware.
# ct_accumulation.nf runs Channel.of("top","bottom","all"), so all three
# accumulation_{top,bottom,all}_<scheme>_aggregated_results.csv files are
# staged into accum_dir and all three are read here.
# Score: accumulation_score = 1 - mean(w_i * p_i) / n_available_schemes
#   where w_i are scheme weights and p_i are PValueEmpirical per scheme.
# Significance: top 5% of accumulation_score (genome-wide empirical threshold).
#
# Returns a tibble with Gene + accum_cct_p<suffix>/accum_fdr<suffix>/
# accum_significant<suffix>/accum_pval_<scheme><suffix> (suffix = "" for "all",
# "_top"/"_bottom" otherwise) -- suffixed so the three direction-specific
# results can full_join onto gene_scores without colliding.
compute_accum_significance <- function(accum_dir, direction, suffix) {
  scheme_names <- c("us", "gs4", "gs3", "gs2", "gs1")
  empty <- tibble(
    Gene = character(),
    !!paste0("accum_cct_p", suffix) := numeric(),
    !!paste0("accum_fdr", suffix)      := numeric(),
    !!paste0("accum_significant", suffix) := logical()
  )

  files_found <- list.files(accum_dir, pattern = paste0("^accumulation_", direction, "_"), full.names = TRUE)
  if (length(files_found) == 0) {
    cat(sprintf("Accumulation (%s): no accumulation_%s_* files found, skipping\n", direction, direction))
    return(list(df = empty, ok = FALSE))
  }

  cat(sprintf("Loading accumulation (%s) from: %s\n", direction, accum_dir))
  accum_pval_df <- NULL  # will hold Gene + one pval col per scheme

  for (scheme in scheme_names) {
    pattern <- paste0("accumulation_", direction, "_", scheme, "_aggregated_results.csv")
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
      select(Gene, !!paste0("accum_pval_", scheme, suffix) := all_of(pval_col))

    if (is.null(accum_pval_df)) {
      accum_pval_df <- scheme_pvals
    } else {
      accum_pval_df <- accum_pval_df %>% full_join(scheme_pvals, by = "Gene")
    }
  }

  if (is.null(accum_pval_df)) return(list(df = empty, ok = FALSE))

  pval_cols <- grep(paste0("^accum_pval_.*", suffix, "$"), names(accum_pval_df), value = TRUE)

  # Cauchy Combination Test (CCT / ACAT) across the available per-group schemes,
  # collapsing the per-scheme accumulation p-values into one value per gene.
  # Stat: T = sum(w_i * tan((0.5 - p_i) * pi)), p_CCT = pcauchy(T, lower.tail = FALSE).
  #
  # CCT is the right combiner here because the five schemes (US, GS4, GS3, GS2,
  # GS1) are independent partitions of the amino acids but test the same
  # physical positions: a position counted under one scheme is frequently
  # counted under others, making the per-scheme p-values positively correlated.
  # CCT's heavy Cauchy tail keeps its null valid under arbitrary dependence, so
  # no independence assumption is required.
  # The same combiner is applied in accum_gene_lists.nf and 10.Accumulation_report.Rmd.
  out <- accum_pval_df %>%
    rowwise() %>%
    mutate(
      !!paste0("accum_cct_p", suffix) := {
        pvals <- c_across(all_of(pval_cols))
        valid <- !is.na(pvals)
        if (sum(valid) == 0) NA_real_
        else if (all(pvals[valid] >= 1)) 1.0
        else {
          ps <- pmin(pmax(pvals[valid], 1e-15), 1 - 1e-15)
          w <- 1 / length(ps)
          stat <- sum(w * tan((0.5 - ps) * pi))
          pcauchy(stat, lower.tail = FALSE)
        }
      }
    ) %>%
    ungroup() %>%
    select(Gene, !!paste0("accum_cct_p", suffix), all_of(pval_cols))

  fp_col <- paste0("accum_cct_p", suffix)
  # BH FDR on genes with at least one observed CAAS (p < 1 strictly). Genes
  # with no CAAS in any group have accum_cct_p = 1 by construction - they
  # are background members only and must not enter the FDR denominator.
  tested <- !is.na(out[[fp_col]]) & out[[fp_col]] < 1
  fdr_q  <- rep(NA_real_, nrow(out))
  if (any(tested)) fdr_q[tested] <- p.adjust(out[[fp_col]][tested], method = "BH")
  out[[paste0("accum_fdr", suffix)]] <- fdr_q
  out[[paste0("accum_significant", suffix)]] <- !is.na(fdr_q) & fdr_q < 0.05

  cat(sprintf("  Accumulation (%s): %d genes, %d significant (Cauchy CCT p, BH FDR < 0.05)\n",
              direction, nrow(out), sum(out[[paste0("accum_significant", suffix)]], na.rm = TRUE)))
  list(df = out, ok = TRUE)
}

has_accum_dir <- file_exists(accum_dir) && dir.exists(accum_dir)
if (!has_accum_dir) cat("Accumulation: directory not available, skipping\n")

.accum_all    <- if (has_accum_dir) compute_accum_significance(accum_dir, "all",    "")        else list(df = tibble(Gene = character(), accum_cct_p = numeric(), accum_fdr = numeric(), accum_significant = logical()), ok = FALSE)
.accum_top    <- if (has_accum_dir) compute_accum_significance(accum_dir, "top",    "_top")    else list(df = tibble(Gene = character(), accum_cct_p_top = numeric(), accum_fdr_top = numeric(), accum_significant_top = logical()), ok = FALSE)
.accum_bottom <- if (has_accum_dir) compute_accum_significance(accum_dir, "bottom", "_bottom") else list(df = tibble(Gene = character(), accum_cct_p_bottom = numeric(), accum_fdr_bottom = numeric(), accum_significant_bottom = logical()), ok = FALSE)

has_accum <- .accum_all$ok
gene_rand <- .accum_all$df %>%
  full_join(.accum_top$df,    by = "Gene") %>%
  full_join(.accum_bottom$df, by = "Gene")

# ── 4c. Gene RERConverge Score (optional) ───────────────────────────────────
has_rer <- file_exists(rer_file)
if (has_rer) {
  cat("Loading RER summary:", rer_file, "\n")
  rer <- read_tsv(rer_file, show_col_types = FALSE)
  
  # Normalize names to lowercase for case-insensitivity
  names(rer) <- tolower(names(rer))
  if (!"gene" %in% names(rer)) {
    names(rer)[1] <- "gene"
  }

  # Use p.perm if available, otherwise p.adj / p.value / pval
  pval_col <- intersect(c("p.perm", "p.adj", "p.value", "pvalue", "pval", "p"), names(rer))[1]
  if (is.na(pval_col)) {
    pval_col <- names(rer)[grepl("^p", names(rer))][1]
    if (is.na(pval_col)) {
      rer$p.adj <- 1.0
      pval_col <- "p.adj"
    }
  }
  cat(sprintf("  Using %s for rer_min_pval\n", pval_col))

  # Find Rho case-insensitively
  rho_col <- intersect(c("rho", "r", "stat"), names(rer))[1]
  if (is.na(rho_col)) {
    rho_col <- names(rer)[grepl("rho", names(rer))][1]
    if (is.na(rho_col)) {
      rer$rho <- 0.0
      rho_col <- "rho"
    }
  }
  cat(sprintf("  Using %s for rer_rho\n", rho_col))

  gene_rer <- rer %>%
    filter(!is.na(.data[[pval_col]])) %>%
    mutate(
      rer_min_pval     = as.numeric(.data[[pval_col]]),
      rer_significant  = rer_min_pval <= 0.05,
      rer_rho          = as.numeric(.data[[rho_col]]),
      rer_acceleration = case_when(
        is.na(rer_rho) ~ NA_character_,
        rer_rho > 0    ~ "accelerated",
        rer_rho < 0    ~ "decelerated",
        TRUE           ~ "neutral"
      )
    ) %>%
    select(Gene = gene, rer_min_pval, rer_significant, rer_rho, rer_acceleration) %>%
    mutate(Gene = as.character(Gene))

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
# full_join (not left_join): gene_caas only covers genes with a CAAS-detected
# position, but gene_rand/gene_rer/gene_fade each have their own, larger gene
# universe (e.g. RERConverge scores every gene with a sufficiently-populated
# gene tree). A left_join rooted at gene_caas silently drops every gene that
# has RER/FADE/accumulation evidence but no CAAS position - which then never
# reaches gene_scores.tsv or the per-module fcs_stats_*.tsv files, truncating
# their FCS universe to the CAAS gene set instead of each module's own.
gene_scores <- gene_caas

if (nrow(gene_rand) > 0) {
  gene_scores <- gene_scores %>% full_join(gene_rand, by = "Gene")
} else {
  gene_scores$accum_cct_p           <- NA_real_
  gene_scores$accum_fdr                <- NA_real_
  gene_scores$accum_significant        <- NA
  gene_scores$accum_cct_p_top       <- NA_real_
  gene_scores$accum_fdr_top            <- NA_real_
  gene_scores$accum_significant_top    <- NA
  gene_scores$accum_cct_p_bottom    <- NA_real_
  gene_scores$accum_fdr_bottom         <- NA_real_
  gene_scores$accum_significant_bottom <- NA
}

if (nrow(gene_rer) > 0) {
  gene_scores <- gene_scores %>% full_join(gene_rer, by = "Gene")
} else {
  gene_scores$rer_min_pval      <- NA_real_
  gene_scores$rer_significant   <- NA
  gene_scores$rer_rho           <- NA_real_
  gene_scores$rer_acceleration  <- NA_character_
}

if (nrow(gene_fade) > 0) {
  gene_scores <- gene_scores %>% full_join(gene_fade, by = "Gene")
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

# Correlations use gene_caas_score as the primary axis. Accumulation, RER and FADE
# are represented by their native significance (accum_cct_p / RER p / FADE BF),
# not numeric score axes, so they are not included here.
score_cols <- c("gene_caas_score")

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

# Position scores
pos_out <- pos_scores %>%
  select(Gene, Position, any_of("pvalue_boot"),
         asr_score, any_of(c("mrca_diversity", "derived_agreement", "conservation_gate", "core",
                             "core_perside_pooled")),
         any_of("phen_score"), n_schemes, any_of("scheme_set"),
         any_of(c("n_hypotheses", "supporting_hypotheses")), CAAS_score,
         change_side,
         any_of(c("caas", "change_top", "change_bottom"))) %>%
  arrange(desc(CAAS_score))

write_tsv(pos_out, "position_scores.tsv")
cat(sprintf("  position_scores.tsv: %d rows\n", nrow(pos_out)))

# Gene scores
gene_out <- gene_scores %>%
  select(
    Gene,
    n_positions, n_positions_top, n_positions_bottom,
    gene_caas_score_top_all, gene_caas_score_bottom_all,
    gene_caas_score, gene_caas_score_top, gene_caas_score_bottom,
    any_of(c("accum_cct_p", "accum_fdr", "accum_significant",
             "accum_pval_us", "accum_pval_gs4", "accum_pval_gs3",
             "accum_pval_gs2", "accum_pval_gs1",
             "accum_cct_p_top", "accum_fdr_top", "accum_significant_top",
             "accum_pval_us_top", "accum_pval_gs4_top", "accum_pval_gs3_top",
             "accum_pval_gs2_top", "accum_pval_gs1_top",
             "accum_cct_p_bottom", "accum_fdr_bottom", "accum_significant_bottom",
             "accum_pval_us_bottom", "accum_pval_gs4_bottom", "accum_pval_gs3_bottom",
             "accum_pval_gs2_bottom", "accum_pval_gs1_bottom")),
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
  score_global     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score"))),
  score_top        = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_top_all"))),
  score_bottom     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_bottom_all"))),
  flag_fade_top    = .istrue(.col(gene_scores, "fade_significant_top")),
  flag_fade_bottom = .istrue(.col(gene_scores, "fade_significant_bottom")),
  flag_rer_acc     = .istrue(.col(gene_scores, "rer_significant")) & grepl("acc", .rer_dir),
  flag_rer_decc    = .istrue(.col(gene_scores, "rer_significant")) & grepl("dec", .rer_dir),
  flag_accum        = .istrue(.col(gene_scores, "accum_significant")),
  flag_accum_top    = .istrue(.col(gene_scores, "accum_significant_top")),
  flag_accum_bottom = .istrue(.col(gene_scores, "accum_significant_bottom"))
) %>%
  mutate(flag_fade = flag_fade_top | flag_fade_bottom)
write_tsv(fcs_stats, "fcs_stats.tsv")
cat(sprintf("  fcs_stats.tsv: %d genes (%d top, %d bottom)\n",
            nrow(fcs_stats), sum(!is.na(fcs_stats$score_top)),
            sum(!is.na(fcs_stats$score_bottom))))

# ── Per-module FCS RANKING files (rankings only) ─────────────────────────────
# Each carries only its module's score_<ranking> columns. The cross-module
# leading-edge annotation (CAAS gates, top/bottom percentiles, FADE, accum) is
# supplied separately to the FCS report via annot_file = fcs_stats.tsv (above).
# This lets SCORING render an FCS report per module - ranked by that module's own
# statistic, annotated with every module's evidence - downstream of the join.
.nonempty <- function(col) .has(gene_scores, col) && any(!is.na(gene_scores[[col]]))

if (has_rer) {
  .rp <- pmax(suppressWarnings(as.numeric(.col(gene_scores, "rer_min_pval"))), 1e-300)
  .rr <- suppressWarnings(as.numeric(.col(gene_scores, "rer_rho")))
  write_tsv(tibble(
    gene               = gene_scores$Gene,
    score_global       = sign(.rr) * -log10(.rp),
    score_accelerating = ifelse(.rr > 0, -log10(.rp), 0),
    score_decelerating = ifelse(.rr < 0, -log10(.rp), 0)
  ), "fcs_stats_rer.tsv")
  cat("  fcs_stats_rer.tsv written (RER rankings)\n")
}

if (.nonempty("fade_max_bf_top") || .nonempty("fade_max_bf_bottom")) {
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_top    = suppressWarnings(as.numeric(.col(gene_scores, "fade_max_bf_top"))),
    score_bottom = suppressWarnings(as.numeric(.col(gene_scores, "fade_max_bf_bottom")))
  ), "fcs_stats_fade.tsv")
  cat("  fcs_stats_fade.tsv written (FADE BF rankings)\n")
}

if (.nonempty("accum_cct_p")) {
  # FCS ranks descending (higher = more significant) while the CCT p runs the
  # other way (lower = more significant), so the ranking axis is
  # -log10(accum_cct_p), derived inline rather than stored as its own column.
  .afp <- suppressWarnings(as.numeric(.col(gene_scores, "accum_cct_p")))
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_global = ifelse(is.na(.afp), NA_real_, -log10(pmax(.afp, 1e-300)))
  ), "fcs_stats_accum.tsv")
  cat("  fcs_stats_accum.tsv written (accumulation rankings, -log10 CCT p)\n")
}

# Correlations
write_tsv(corr_results, "gene_correlations.tsv")

# =============================================================================
# 7. EXPORT PERCENTILE-RANKED GENE LISTS FOR STRING
# =============================================================================

dir.create("gene_lists", showWarnings = FALSE)

# Define the 12 percentile slices (top/bottom/global x 25/10/5/1%), mirroring
# posenrich_enrich.py's top-fraction-of-scored-positions cutoff at gene
# granularity: rank genes with a defined score desc, keep the top frac% of
# THAT ranked set (not the full background) -- same selection axis as
# 16.Position_enrichment_report.Rmd, just at genes instead of positions.
# Percentile slices, not a significance gate: gate_sig already reports
# significance elsewhere, and a significance-gated foreground close to the
# whole gene universe is a poor input for interaction-density tests.
#
# "global" slices (added alongside top/bottom): ranked on gene_caas_score, the
# undirected 90th-percentile-of-all-positions aggregate -- the gene-level
# analog of CAAS FCS's own undirected "global" ranking (fcs_all_results.tsv's
# ranking=="global" rows). fade_sig is NA for global (FADE has no undirected
# significance flag); is_fade in the export loop below already treats a
# missing/NA fade_sig as FALSE.
slices_def <- list(
  list(name = "top25",    col = "gene_caas_score_top_all",    frac = 0.25, direction = "top",    fade_sig = "fade_significant_top"),
  list(name = "top10",    col = "gene_caas_score_top_all",    frac = 0.10, direction = "top",    fade_sig = "fade_significant_top"),
  list(name = "top5",     col = "gene_caas_score_top_all",    frac = 0.05, direction = "top",    fade_sig = "fade_significant_top"),
  list(name = "top1",     col = "gene_caas_score_top_all",    frac = 0.01, direction = "top",    fade_sig = "fade_significant_top"),

  list(name = "bottom25", col = "gene_caas_score_bottom_all", frac = 0.25, direction = "bottom", fade_sig = "fade_significant_bottom"),
  list(name = "bottom10", col = "gene_caas_score_bottom_all", frac = 0.10, direction = "bottom", fade_sig = "fade_significant_bottom"),
  list(name = "bottom5",  col = "gene_caas_score_bottom_all", frac = 0.05, direction = "bottom", fade_sig = "fade_significant_bottom"),
  list(name = "bottom1",  col = "gene_caas_score_bottom_all", frac = 0.01, direction = "bottom", fade_sig = "fade_significant_bottom"),

  list(name = "global25", col = "gene_caas_score", frac = 0.25, direction = "global", fade_sig = NA_character_),
  list(name = "global10", col = "gene_caas_score", frac = 0.10, direction = "global", fade_sig = NA_character_),
  list(name = "global5",  col = "gene_caas_score", frac = 0.05, direction = "global", fade_sig = NA_character_),
  list(name = "global1",  col = "gene_caas_score", frac = 0.01, direction = "global", fade_sig = NA_character_)
)

cat("\n─── Exporting 12 STRING percentile ranked lists (top/bottom/global x 25/10/5/1%) ──\n")
for (slice in slices_def) {
  col_name <- slice$col
  file_name <- sprintf("gene_lists/slice_%s.tsv", slice$name)

  # Rank all non-NA genes by score (desc, stable Gene tiebreak for determinism).
  ranked <- gene_out %>%
    filter(!is.na(.data[[col_name]])) %>%
    select(Gene, score = all_of(col_name)) %>%
    arrange(desc(score), Gene)

  n_total <- nrow(ranked)
  if (n_total == 0) {
    # Write empty file to prevent Nextflow missing output crashes
    write_tsv(tibble(Gene = character(), score = numeric(), is_fade = logical(), is_rer = logical(), is_accum = logical()), file_name)
    cat(sprintf("  %s: empty (0 genes) exported\n", file_name))
    next
  }

  n_keep <- max(1, round(slice$frac * n_total))
  slice_df <- ranked %>% slice_head(n = n_keep)

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
        if (slice$direction == "top") {
          is_sig & !is.na(rho_val) & rho_val > 0
        } else if (slice$direction == "bottom") {
          is_sig & !is.na(rho_val) & rho_val < 0
        } else {
          is_sig
        }
      } else FALSE,
      # Direction-matched, mirroring is_fade above -- a "top" slice checks
      # accum_significant_top (accumulation_top_* pooled positions), not the
      # non-directional accum_significant (accumulation_all_* pooled).
      is_accum = {
        acc_col <- if (slice$direction == "top") "accum_significant_top"
                   else if (slice$direction == "bottom") "accum_significant_bottom"
                   else "accum_significant"
        if (acc_col %in% names(gene_out)) {
          gene_out[[acc_col]][match(Gene, gene_out$Gene)] %in% TRUE
        } else FALSE
      }
    )
  # slice_df inherits ranked's desc(score) order via slice_head -- already sorted.

  write_tsv(slice_df, file_name)
  cat(sprintf("  %s: %d/%d genes exported (top %.0f%%)\n", file_name, nrow(slice_df), n_total, 100 * slice$frac))
}

# =============================================================================
# 7b. EXPORT PERCENTILE-RANKED POSITION LISTS FOR POSENRICH/REPORTS
# =============================================================================
# Position-level analog of the gene_lists/ slices above, published so
# posenrich_enrich.py and the FCS/comparison reports (12/14/15) all read ONE
# canonical top/bottom/global x 25/10/5/1% position ranking instead of each
# independently re-deriving their own cutoff over position_scores.tsv (which
# invited drift -- see e.g. a single outlier position ranking in a gene's top
# 1% of ALL positions while that gene's own aggregate score doesn't crack the
# gene-level top 5%; both were "correct" under their own re-derivation, just
# never guaranteed to agree). Mirrors posenrich_enrich.py's direction_scores()
# + top-frac cutoff exactly (subworkflows/ENRICHMENT/local/src/posenrich_enrich.py):
# direction filters on change_side (top: {top,both}, bottom: {bottom,both},
# global: unfiltered), ranked over positions with CAAS_score > 0 only, cut at
# round(frac * n_scored), stable sort by (desc score, Gene, Position).
dir.create("position_lists", showWarnings = FALSE)

pos_slices_def <- list(
  list(name = "top25",    frac = 0.25, direction = "top"),
  list(name = "top10",    frac = 0.10, direction = "top"),
  list(name = "top5",     frac = 0.05, direction = "top"),
  list(name = "top1",     frac = 0.01, direction = "top"),
  list(name = "bottom25", frac = 0.25, direction = "bottom"),
  list(name = "bottom10", frac = 0.10, direction = "bottom"),
  list(name = "bottom5",  frac = 0.05, direction = "bottom"),
  list(name = "bottom1",  frac = 0.01, direction = "bottom"),
  list(name = "global25", frac = 0.25, direction = "global"),
  list(name = "global10", frac = 0.10, direction = "global"),
  list(name = "global5",  frac = 0.05, direction = "global"),
  list(name = "global1",  frac = 0.01, direction = "global")
)

cat("\n─── Exporting 12 position-level percentile ranked lists (top/bottom/global x 25/10/5/1%) ──\n")
for (slice in pos_slices_def) {
  file_name <- sprintf("position_lists/slice_%s.tsv", slice$name)

  sub <- if (slice$direction == "top") {
    pos_out %>% dplyr::filter(change_side %in% c("top", "both"))
  } else if (slice$direction == "bottom") {
    pos_out %>% dplyr::filter(change_side %in% c("bottom", "both"))
  } else {
    pos_out
  }
  scored <- sub %>%
    dplyr::filter(!is.na(CAAS_score) & CAAS_score > 0) %>%
    dplyr::arrange(dplyr::desc(CAAS_score), Gene, Position)

  n_scored <- nrow(scored)
  if (n_scored == 0) {
    write_tsv(tibble(Gene = character(), Position = integer(), CAAS_score = numeric()), file_name)
    cat(sprintf("  %s: empty (0 scored positions) exported\n", file_name))
    next
  }

  n_keep <- max(1, round(slice$frac * n_scored))
  slice_df <- scored %>% dplyr::slice_head(n = n_keep) %>% dplyr::select(Gene, Position, CAAS_score)

  write_tsv(slice_df, file_name)
  cat(sprintf("  %s: %d/%d scored positions exported (top %.0f%%)\n", file_name, nrow(slice_df), n_scored, 100 * slice$frac))
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
  # Direction-matched, mirroring fade_top/fade_bottom and rer_accel/rer_decel
  # above -- "top"/"bottom" test against accumulation_top_*/accumulation_bottom_*
  # (positions restricted to that direction), "global" against the
  # accumulation_all_* pool (all non-none positions, non-directional).
  if (dir == "top" && "accum_significant_top" %in% names(gene_out) &&
      any(gene_out$accum_significant_top %in% TRUE))
    tool_sigs[["accum_top"]] <- gene_out %>%
      filter(accum_significant_top == TRUE) %>% pull(Gene)
  if (dir == "bottom" && "accum_significant_bottom" %in% names(gene_out) &&
      any(gene_out$accum_significant_bottom %in% TRUE))
    tool_sigs[["accum_bottom"]] <- gene_out %>%
      filter(accum_significant_bottom == TRUE) %>% pull(Gene)
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
cat("  CAAS Scoring - Complete\n")
cat("═══════════════════════════════════════════════════════════════\n")
