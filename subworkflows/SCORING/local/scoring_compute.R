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
#     --fade_top   <fade_summary_top.tsv> \
#     --fade_bottom <fade_summary_bottom.tsv> \
#     --rer        <rerconverge_summary.tsv> \
#     --accum_dir  <directory_with_accumulation_CSVs> \
#     --top_pct    0.10 \
#     --top5_pct   0.05 \
#     --top1_pct   0.01 \
#     --gene_top_pct  0.10 \
#     --gene_top5_pct 0.05 \
#     --gene_top1_pct 0.01 \
#     --molerate_top    <molerate_summary_top.tsv> \
#     --molerate_bottom <molerate_summary_bottom.tsv>
#
# Outputs (in working directory):
#   position_scores.tsv      — per Gene×Position scores
#   gene_scores.tsv          — per Gene scores
#   gene_correlations.tsv    — pairwise correlations between gene scores
#   pos_gene_lists/*.txt     — gene lists for position-level ORA
#   gene_gene_lists/*.txt    — gene lists for gene-level ORA
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
rer_file             <- parse_arg("--rer")
accum_dir            <- parse_arg("--accum_dir")
molerate_top_file    <- parse_arg("--molerate_top")
molerate_bottom_file <- parse_arg("--molerate_bottom")
stress_enabled_raw   <- parse_arg("--stress", "false")
stress_top_n         <- as.integer(parse_arg("--stress_top_n", "25"))
top_pct           <- as.numeric(parse_arg("--top_pct",  "0.10"))
top5_pct          <- as.numeric(parse_arg("--top5_pct", "0.05"))
top1_pct          <- as.numeric(parse_arg("--top1_pct", "0.01"))
gene_top_pct      <- as.numeric(parse_arg("--gene_top_pct",  "0.10"))
gene_top5_pct     <- as.numeric(parse_arg("--gene_top5_pct", "0.05"))
gene_top1_pct     <- as.numeric(parse_arg("--gene_top1_pct", "0.01"))
stress_enabled    <- tolower(as.character(stress_enabled_raw)) %in% c("true", "1", "yes")
if (!is.finite(stress_top_n) || is.na(stress_top_n) || stress_top_n < 1) stress_top_n <- 25
direction         <- parse_arg("--direction", "both")
stopifnot(direction %in% c("top", "bottom", "both"))

file_exists <- function(f) {
  !is.null(f) && f != "" && !grepl("^NO_", basename(f)) && file.exists(f)
}

dir.create("pos_gene_lists",  showWarnings = FALSE)
dir.create("gene_gene_lists", showWarnings = FALSE)

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

# ── 2a. Variability score (deciles of hypergeometric p-value) ────────────────
# Computed across ALL rows, then collapsed per Gene×Position later.
df <- df %>%
  mutate(
    variability_score = decile_score(Pvalue),
    pattern_score     = decile_score(pvalue_boot)
  )

# ── 2b. Biochem score (cumulative scheme weights × variability × pattern) ────
scheme_weights <- c(US = 0.5, GS4 = 0.2, GS3 = 0.1, GS2 = 0.1, GS1 = 0.1, GS0 = 0)
scheme_weights_flat <- c(US = 0.35, GS4 = 0.2, GS3 = 0.15, GS2 = 0.15, GS1 = 0.15, GS0 = 0)
scheme_weights_us_heavy <- c(US = 0.60, GS4 = 0.2, GS3 = 0.1, GS2 = 0.05, GS1 = 0.05, GS0 = 0)

df <- df %>%
  mutate(
    scheme_weight = scheme_weights[CAAP_Group],
    scheme_weight_flat = scheme_weights_flat[CAAP_Group],
    scheme_weight_us_heavy = scheme_weights_us_heavy[CAAP_Group]
  )

# Compute per-row scheme contribution
df <- df %>%
  mutate(
    scheme_contribution = scheme_weight * variability_score * pattern_score,
    scheme_contribution_mean = scheme_weight * ((variability_score + pattern_score) / 2),
    scheme_contribution_geo = scheme_weight * sqrt(pmax(variability_score * pattern_score, 0)),
    scheme_contribution_flat = scheme_weight_flat * variability_score * pattern_score,
    scheme_contribution_us_heavy = scheme_weight_us_heavy * variability_score * pattern_score
  )

# Detect pair-indexed MRCA posterior columns dynamically (mrca_1_posterior, mrca_2_posterior, ...)
mrca_posterior_cols <- grep("^mrca_\\d+_posterior$", names(df), value = TRUE)
n_pairs <- length(mrca_posterior_cols)
cat(sprintf("  Detected %d pairs (%s)\n", n_pairs, paste(mrca_posterior_cols, collapse = ", ")))

# Aggregate to Gene×Position level: sum scheme contributions (cumulative),
# and take the first value for fields that are constant within Gene×Position.
pos_scores <- df %>%
  group_by(Gene, Position) %>%
  summarise(
    biochem_score         = sum(scheme_contribution, na.rm = TRUE),
    scheme_coverage       = sum(scheme_weight, na.rm = TRUE),
    biochem_mul_current   = sum(scheme_contribution, na.rm = TRUE),
    biochem_mean          = sum(scheme_contribution_mean, na.rm = TRUE),
    biochem_geo           = sum(scheme_contribution_geo, na.rm = TRUE),
    biochem_flat_weights  = sum(scheme_contribution_flat, na.rm = TRUE),
    biochem_us_heavy      = sum(scheme_contribution_us_heavy, na.rm = TRUE),
    biochem_saturated     = {
      contribs <- sort(scheme_contribution[is.finite(scheme_contribution)], decreasing = TRUE)
      if (length(contribs) == 0) 0 else sum(contribs / (2 ^ (seq_along(contribs) - 1)))
    },
    # Take representative values from the row set (constant within Gene×Position)
    Pvalue                = first(Pvalue),
    pvalue_boot           = first(pvalue_boot),
    variability_score     = first(variability_score),
    pattern_score         = first(pattern_score),
    is_conserved_meta     = first(is_conserved_meta),
    asr_is_conserved      = first(asr_is_conserved),
    asr_root_conserved    = first(asr_root_conserved),
    conserved_pair        = first(conserved_pair),
    all_mrca_posterior    = first(all_mrca_posterior),
    across(all_of(mrca_posterior_cols), \(x) first(x)),
    change_top            = first(change_top),
    change_bottom         = first(change_bottom),
    change_side           = first(change_side),
    parallel_top          = first(parallel_top),
    parallel_bottom       = first(parallel_bottom),
    .groups = "drop"
  )

cat(sprintf("  %d unique positions after aggregation\n", nrow(pos_scores)))

# ── 2c. ASR score ────────────────────────────────────────────────────────────
# Pair-weighted formula (w = 1 / n_pairs):
#
#   is_conserved_meta = FALSE                          → 1.0  (novel change)
#   is_conserved_meta = TRUE & asr_is_conserved = FALSE → 0.0  (unconfirmed conservation)
#   is_conserved_meta = TRUE & asr_is_conserved = TRUE  →
#     Σ_{divergent pairs}  w
#   + Σ_{conserved pairs}  min((w × P_MRCA(i) + w × P_root × I(root_cons)) / 2,  0.9 × w)
#
# The 0.9w cap ensures conserved positions never reach the 1.0 ceiling of novel changes.
# conserved_pair may be a single id ("2") or comma-separated ("2,4") for future multi-pair support.
compute_asr_score <- function(is_conserved_meta_vec, asr_is_conserved_vec,
                              asr_root_conserved_vec, conserved_pair_vec,
                              all_mrca_posterior_vec, n_pairs, pos_data) {

  w <- 1.0 / n_pairs

  parse_bool <- function(x) {
    if (is.logical(x)) return(x)
    tolower(as.character(x)) %in% c("true", "1", "yes")
  }

  sapply(seq_along(is_conserved_meta_vec), function(i) {
    is_cons   <- parse_bool(is_conserved_meta_vec[i])
    asr_cons  <- parse_bool(asr_is_conserved_vec[i])
    root_cons <- parse_bool(asr_root_conserved_vec[i])

    if (!is_cons)  return(1.0)
    if (!asr_cons) return(0.0)

    # Parse conserved pair IDs: "2" → [2],  "2,4" → [2, 4]
    conserved_ids <- suppressWarnings(
      as.integer(trimws(strsplit(as.character(conserved_pair_vec[i]), ",")[[1]]))
    )
    conserved_ids <- conserved_ids[!is.na(conserved_ids) & conserved_ids >= 1]
    if (length(conserved_ids) == 0) return(1.0)

    n_divergent <- n_pairs - length(conserved_ids)
    score       <- n_divergent * w

    all_mrca_p <- suppressWarnings(as.numeric(all_mrca_posterior_vec[i]))
    root_p     <- if (root_cons && !is.na(all_mrca_p)) all_mrca_p else 0.0

    for (pid in conserved_ids) {
      mcol    <- paste0("mrca_", pid, "_posterior")
      mrca_p  <- if (mcol %in% names(pos_data)) suppressWarnings(as.numeric(pos_data[[mcol]][i])) else 0.0
      if (is.na(mrca_p)) mrca_p <- 0.0
      score   <- score + min((w * mrca_p + w * root_p) / 2, 0.9 * w)
    }

    return(min(score, 1.0))
  })
}

pos_scores_ref <- pos_scores  # capture before pipe for mrca column lookup inside sapply

pos_scores <- pos_scores %>%
  mutate(
    asr_score = compute_asr_score(
      is_conserved_meta, asr_is_conserved, asr_root_conserved,
      conserved_pair, all_mrca_posterior,
      n_pairs, pos_scores_ref
    )
  )

# ── 2d. Convergence score ───────────────────────────────────────────────────
compute_convergence_score <- function(change_top, change_bottom) {
  side_value <- function(side) {
    case_when(
      side == "convergent"  ~ 1.0,
      side == "codivergent" ~ 0.5,
      side == "divergent"   ~ 0.25,
      TRUE                  ~ NA_real_  # ambiguous / no_change → excluded
    )
  }

  top_val    <- side_value(change_top)
  bottom_val <- side_value(change_bottom)

  mean_change <- ifelse(
    is.na(top_val) & is.na(bottom_val), 0,
    rowMeans(cbind(top_val, bottom_val), na.rm = TRUE)
  )

  completeness_penalty <- ifelse(!is.na(top_val) & !is.na(bottom_val), 1.0, 0.75)

  mean_change * completeness_penalty
}

pos_scores <- pos_scores %>%
  mutate(
    convergence_score = compute_convergence_score(change_top, change_bottom)
  )

# ── 2e. Parallel score ──────────────────────────────────────────────────────
# Parallelism coefficients (categorical, no MRCA-posterior weighting):
#   Non-parallel = 1.0 (highest): independent changes are the most informative
#     signal of convergence — identical substitutions in independent lineages
#     are less surprising and may reflect shared selective constraints.
#   Mixed = 0.75, Parallel = 0.5.
# Completeness penalty: 1.0 when both sides are assessable, 0.75 when only one.
# MRCA reconstruction confidence is already captured by asr_score; applying it
# here again would double-deflate positions with uncertain reconstructions and
# collapse the effective range to ~0–0.375 instead of 0–1.
# FADE and PrimateAI-3D provide complementary position-level characterization.
compute_parallel_score <- function(parallel_top, parallel_bottom,
                                   change_top,   change_bottom) {
  parallel_coeff <- function(p) {
    case_when(
      is.na(p) | p == "" | tolower(p) == "na" ~ NA_real_,
      tolower(p) == "non-parallel" | p == "False" | p == "FALSE" ~ 1.0,
      tolower(p) == "mixed"        ~ 0.75,
      tolower(p) == "parallel"     | p == "True" | p == "TRUE"  ~ 0.5,
      grepl("^parallel", tolower(p)) ~ 0.5,
      TRUE                          ~ NA_real_
    )
  }

  # Side is assessable only if change is convergent/codivergent/divergent
  side_assessable <- function(change) {
    change %in% c("convergent", "codivergent", "divergent")
  }

  top_assessable    <- side_assessable(change_top)
  bottom_assessable <- side_assessable(change_bottom)

  top_coeff    <- parallel_coeff(parallel_top)
  bottom_coeff <- parallel_coeff(parallel_bottom)

  # Use coefficients directly — no MRCA multiplier (see comment above)
  top_side_score    <- ifelse(top_assessable,    top_coeff, NA_real_)
  bottom_side_score <- ifelse(bottom_assessable, bottom_coeff, NA_real_)

  mean_parallel <- ifelse(
    is.na(top_side_score) & is.na(bottom_side_score), 0,
    rowMeans(cbind(top_side_score, bottom_side_score), na.rm = TRUE)
  )

  completeness_penalty <- ifelse(!is.na(top_side_score) & !is.na(bottom_side_score), 1.0, 0.75)

  mean_parallel * completeness_penalty
}

pos_scores <- pos_scores %>%
  mutate(
    parallel_score = compute_parallel_score(
      parallel_top, parallel_bottom, change_top, change_bottom
    )
  )

# ── 2g. FADE (moved to gene-level — see section 3c) ─────────────────────────
# FADE operates at gene level (max Bayes Factor per gene across sites).
# It is loaded here so the data is available, but NOT broadcast to positions.
has_fade <- file_exists(fade_top_file) || file_exists(fade_bottom_file)
if (has_fade) {
  cat(sprintf("Loading FADE summaries (direction='%s')\n", direction))
  # For a directional run use only the relevant file; for 'both' take max across directions.
  fade_file_for_dir <- if (direction == "top")    fade_top_file
                       else if (direction == "bottom") fade_bottom_file
                       else NULL  # 'both' → pool below
  if (!is.null(fade_file_for_dir) && file_exists(fade_file_for_dir)) {
    fade_combined <- read_tsv(fade_file_for_dir, show_col_types = FALSE) %>%
      select(Gene = gene, fade_max_bf = max_bf)
  } else {
    # Pool both directions: keep the gene's highest Bayes Factor across directions
    fade_dfs <- list()
    if (file_exists(fade_top_file))
      fade_dfs[["top"]]    <- read_tsv(fade_top_file,    show_col_types = FALSE) %>%
        select(Gene = gene, max_bf)
    if (file_exists(fade_bottom_file))
      fade_dfs[["bottom"]] <- read_tsv(fade_bottom_file, show_col_types = FALSE) %>%
        select(Gene = gene, max_bf)
    fade_combined <- bind_rows(fade_dfs) %>%
      group_by(Gene) %>%
      summarise(fade_max_bf = max(max_bf, na.rm = TRUE), .groups = "drop")
  }
  cat(sprintf("  FADE: %d genes loaded\n", nrow(fade_combined)))
} else {
  cat("FADE: not available, skipping\n")
  fade_combined <- tibble(Gene = character(), fade_max_bf = numeric())
}

# ── 2h. MOLERATE (moved to gene-level — see section 3e) ─────────────────────
# MoleRate is gene-level. We take min(pval_pp_vs_prop) across top and bottom
# directions per gene (most significant direction, analogous to FADE max_bf).
has_molerate <- file_exists(molerate_top_file) || file_exists(molerate_bottom_file)
if (has_molerate) {
  cat(sprintf("Loading MoleRate summaries (direction='%s')\n", direction))
  mol_file_for_dir <- if (direction == "top")    molerate_top_file
                      else if (direction == "bottom") molerate_bottom_file
                      else NULL  # 'both' → pool below

  if (!is.null(mol_file_for_dir) && file_exists(mol_file_for_dir)) {
    # Single-direction: pval and log2_ratio from the target file only
    mol_dir_col <- paste0("mol_log2ratio_", direction)
    mol_raw <- read_tsv(mol_file_for_dir, show_col_types = FALSE) %>%
      select(Gene = gene, pval_pp_vs_prop, log2_ratio) %>%
      filter(!is.na(pval_pp_vs_prop))
    mol_combined <- mol_raw %>%
      group_by(Gene) %>%
      summarise(mol_min_pval = min(pval_pp_vs_prop, na.rm = TRUE),
                !!mol_dir_col := first(log2_ratio),
                .groups = "drop")
  } else {
    # Pool: min pval across both directions; keep both log2_ratio columns
    mol_dfs <- list()
    if (file_exists(molerate_top_file))
      mol_dfs[["top"]]    <- read_tsv(molerate_top_file,    show_col_types = FALSE) %>%
        select(Gene = gene, pval_pp_vs_prop, mol_log2ratio_top    = log2_ratio) %>%
        filter(!is.na(pval_pp_vs_prop))
    if (file_exists(molerate_bottom_file))
      mol_dfs[["bottom"]] <- read_tsv(molerate_bottom_file, show_col_types = FALSE) %>%
        select(Gene = gene, pval_pp_vs_prop, mol_log2ratio_bottom = log2_ratio) %>%
        filter(!is.na(pval_pp_vs_prop))
    mol_pvals <- bind_rows(
      if (!is.null(mol_dfs[["top"]]))    select(mol_dfs[["top"]],    Gene, pval_pp_vs_prop),
      if (!is.null(mol_dfs[["bottom"]])) select(mol_dfs[["bottom"]], Gene, pval_pp_vs_prop)
    ) %>%
      group_by(Gene) %>%
      summarise(mol_min_pval = min(pval_pp_vs_prop, na.rm = TRUE), .groups = "drop")
    mol_log2ratio <- tibble(Gene = character(), mol_log2ratio_top = numeric(),
                            mol_log2ratio_bottom = numeric())
    if (!is.null(mol_dfs[["top"]]))
      mol_log2ratio <- mol_log2ratio %>%
        full_join(select(mol_dfs[["top"]], Gene, mol_log2ratio_top), by = "Gene")
    if (!is.null(mol_dfs[["bottom"]]))
      mol_log2ratio <- mol_log2ratio %>%
        full_join(select(mol_dfs[["bottom"]], Gene, mol_log2ratio_bottom), by = "Gene")
    mol_combined <- mol_pvals %>% left_join(mol_log2ratio, by = "Gene")
  }
  cat(sprintf("  MoleRate: %d genes loaded\n", nrow(mol_combined)))
} else {
  cat("MoleRate: not available, skipping\n")
  mol_combined <- tibble(Gene = character(), mol_min_pval = numeric())
}

# ── 2i. CAAS Score (composite) ──────────────────────────────────────────────
# Unweighted mean of the four core components: biochem_score, asr_score,
# convergence_score, parallel_score.
# FADE is gene-level and does NOT contribute to position-level CAAS_score.

pos_scores <- pos_scores %>%
  mutate(
    CAAS_score = rowMeans(
      cbind(biochem_score, asr_score, convergence_score, parallel_score),
      na.rm = TRUE
    )
  )

cat(sprintf("\nPosition-level CAAS_score: min=%.3f, median=%.3f, max=%.3f\n",
            min(pos_scores$CAAS_score, na.rm = TRUE),
            median(pos_scores$CAAS_score, na.rm = TRUE),
            max(pos_scores$CAAS_score, na.rm = TRUE)))

# ── Direction filter ──────────────────────────────────────────────────────────
# Applied after all position-level scores are computed so that genome-wide
# decile ranks (variability_score, pattern_score) remain calibrated against the
# full position pool.  Gene-level scoring and ORA then operate on the filtered set.
if (direction %in% c("top", "bottom")) {
  n_before <- nrow(pos_scores)
  pos_scores <- pos_scores %>% filter(change_side %in% c(direction, "both"))
  cat(sprintf("  Direction filter '%s': %d → %d positions retained\n",
              direction, n_before, nrow(pos_scores)))
}

# =============================================================================
# 2j. POSITION-LEVEL STRESS TESTS (optional)
# =============================================================================
if (stress_enabled) {
  cat("\n─── Position-level scoring stress test ───────────────────────\n")

  stress_df <- pos_scores %>%
    transmute(
      Gene, Position,
      variability_score,
      pattern_score,
      scheme_coverage,
      biochem_score,
      biochem_mul_current,
      biochem_mean,
      biochem_geo,
      biochem_saturated,
      biochem_flat_weights,
      biochem_us_heavy,
      asr_score,
      convergence_score,
      parallel_score,
      CAAS_current = CAAS_score
    )

  stress_df <- stress_df %>%
    mutate(
      CAAS_no_biochem = rowMeans(cbind(asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_no_asr = rowMeans(cbind(biochem_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_no_convergence = rowMeans(cbind(biochem_score, asr_score, parallel_score), na.rm = TRUE),
      CAAS_no_parallel = rowMeans(cbind(biochem_score, asr_score, convergence_score), na.rm = TRUE),
      CAAS_biochem_split = rowMeans(cbind(variability_score, pattern_score, scheme_coverage,
                                          asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_biochem_alt_mean = rowMeans(cbind(biochem_mean, asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_biochem_alt_geo = rowMeans(cbind(biochem_geo, asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_biochem_alt_saturated = rowMeans(cbind(biochem_saturated, asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_biochem_alt_flat = rowMeans(cbind(biochem_flat_weights, asr_score, convergence_score, parallel_score), na.rm = TRUE),
      CAAS_biochem_alt_us_heavy = rowMeans(cbind(biochem_us_heavy, asr_score, convergence_score, parallel_score), na.rm = TRUE)
    )

  latent_loadings <- tibble(model = character(), variable = character(), loading = numeric(),
                            variance_explained = numeric())

  pca_all_cols <- c("biochem_score", "asr_score", "convergence_score", "parallel_score")
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

  pca_no_biochem_cols <- c("variability_score", "pattern_score", "scheme_coverage", "asr_score", "convergence_score", "parallel_score")
  pca_no_biochem_df <- stress_df %>% select(all_of(pca_no_biochem_cols), CAAS_current) %>% drop_na()
  pca_no_biochem_cols_var <- variable_numeric_cols(pca_no_biochem_df, pca_no_biochem_cols)
  if (length(pca_no_biochem_cols_var) < length(pca_no_biochem_cols)) {
    cat(sprintf("  PCA_no_biochem: dropped %d constant/degenerate columns\n",
                length(pca_no_biochem_cols) - length(pca_no_biochem_cols_var)))
  }
  if (nrow(pca_no_biochem_df) >= 3 && length(pca_no_biochem_cols_var) >= 2) {
    pca_nb <- prcomp(pca_no_biochem_df[, pca_no_biochem_cols_var, drop = FALSE], center = TRUE, scale. = TRUE)
    scores <- pca_nb$x[, 1]
    score_corr <- safe_cor(scores, pca_no_biochem_df$CAAS_current, method = "pearson")
    if (!is.na(score_corr) && score_corr < 0) scores <- -scores
    stress_df$PCA_no_biochem <- NA_real_
    stress_df$PCA_no_biochem[complete.cases(stress_df[, pca_no_biochem_cols_var])] <- score_to_decile(scores)
    latent_loadings <- bind_rows(
      latent_loadings,
      tibble(
        model = "PCA_no_biochem",
        variable = rownames(pca_nb$rotation),
        loading = pca_nb$rotation[, 1],
        variance_explained = summary(pca_nb)$importance[2, 1]
      )
    )
  } else {
    cat("  PCA_no_biochem: skipped (fewer than 2 variable columns after filtering)\n")
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

  resid_predictors <- c("asr_score", "convergence_score", "parallel_score")
  resid_df <- stress_df %>% select(biochem_score, all_of(resid_predictors)) %>% drop_na()
  if (nrow(resid_df) >= 5) {
    resid_fit <- lm(as.formula(paste("biochem_score ~", paste(resid_predictors, collapse = " + "))), data = resid_df)
    resid_scores <- resid(resid_fit)
    stress_df$residualized_biochem <- NA_real_
    stress_df$residualized_biochem[complete.cases(stress_df[, c("biochem_score", resid_predictors)])] <- score_to_decile(resid_scores)
    stress_df$CAAS_residualized_biochem <- rowMeans(
      cbind(stress_df$residualized_biochem, stress_df$asr_score, stress_df$convergence_score, stress_df$parallel_score),
      na.rm = TRUE
    )
  }

  analysis_cols <- c(
    "variability_score", "pattern_score", "scheme_coverage", "biochem_score",
    "asr_score", "convergence_score", "parallel_score"
  )

  composite_cols <- c(
    "CAAS_current", "CAAS_no_biochem", "CAAS_no_asr", "CAAS_no_convergence",
    "CAAS_no_parallel", "CAAS_biochem_split", "CAAS_biochem_alt_mean",
    "CAAS_biochem_alt_geo", "CAAS_biochem_alt_saturated",
    "CAAS_biochem_alt_flat", "CAAS_biochem_alt_us_heavy"
  )
  latent_cols <- c("PCA_all_components", "PCA_no_biochem", "Factor1_all_components", "CAAS_residualized_biochem")
  latent_cols <- latent_cols[latent_cols %in% names(stress_df)]
  latent_cols <- latent_cols[vapply(latent_cols, function(col) any(!is.na(stress_df[[col]])), logical(1))]
  if (length(latent_cols) > 0) composite_cols <- c(composite_cols, latent_cols)

  stress_correlations <- pairwise_long(stress_df, c(analysis_cols, composite_cols))
  stress_rank_agreement <- make_rank_matrix(stress_df, composite_cols)
  stress_top_overlap <- make_top_overlap(stress_df, composite_cols, top_n = stress_top_n)

  stress_summary <- tibble(variant = composite_cols) %>%
    rowwise() %>%
    mutate(
      pearson_to_current = safe_cor(stress_df[[variant]], stress_df$CAAS_current, method = "pearson"),
      spearman_to_current = safe_cor(stress_df[[variant]], stress_df$CAAS_current, method = "spearman"),
      mean_abs_rank_shift = {
        r0 <- rank(-stress_df$CAAS_current, ties.method = "average", na.last = "keep")
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        mean(abs(r0 - r1), na.rm = TRUE)
      },
      median_abs_rank_shift = {
        r0 <- rank(-stress_df$CAAS_current, ties.method = "average", na.last = "keep")
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        median(abs(r0 - r1), na.rm = TRUE)
      },
      max_abs_rank_shift = {
        r0 <- rank(-stress_df$CAAS_current, ties.method = "average", na.last = "keep")
        r1 <- rank(-stress_df[[variant]], ties.method = "average", na.last = "keep")
        max(abs(r0 - r1), na.rm = TRUE)
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
# 3. GENE-LEVEL SCORING
# =============================================================================

cat("\n─── Gene-level scoring ────────────────────────────────────────\n")

# ── 3a. Gene CAAS Score: 90th percentile of CAAS_score per gene ─────────────
gene_caas <- pos_scores %>%
  group_by(Gene) %>%
  summarise(
    gene_caas_score = quantile(CAAS_score, 0.90, na.rm = TRUE),
    n_positions     = n(),
    .groups = "drop"
  )

cat(sprintf("  gene_caas_score: %d genes\n", nrow(gene_caas)))

# ── 3b. Gene Randomization Score (optional) ─────────────────────────────────
has_accum <- file_exists(accum_dir) && dir.exists(accum_dir)
if (has_accum) {
  cat("Loading accumulation from:", accum_dir, "\n")

  scheme_names <- c("us", "gs4", "gs3", "gs2", "gs1")
  scheme_w     <- c(0.5,  0.2,   0.1,   0.1,   0.1)
  names(scheme_w) <- scheme_names

  accum_scores <- NULL
  accum_meta   <- NULL

  # File prefix includes direction when running a directional report so that
  # accumulation_top_us_*.csv and accumulation_bottom_us_*.csv don't cross-contaminate.
  accum_pfx <- if (direction %in% c("top", "bottom")) paste0("accumulation_", direction) else "accumulation"

  for (scheme in scheme_names) {
    pattern <- paste0(accum_pfx, "_", scheme, "_aggregated_results.csv")
    f <- list.files(accum_dir, pattern = pattern, full.names = TRUE)
    if (length(f) == 0) {
      cat(sprintf("    %s: file not found, skipping\n", scheme))
      next
    }
    cat(sprintf("    %s: %s\n", scheme, basename(f[1])))
    d <- read_csv(f[1], show_col_types = FALSE)

    # Find the PValueEmpirical column (pattern varies by scheme)
    pval_col <- grep("PValueEmpirical", names(d), value = TRUE)[1]
    if (is.na(pval_col)) {
      cat(sprintf("    %s: no PValueEmpirical column found, skipping\n", scheme))
      next
    }

    scheme_df <- d %>%
      select(Gene, pval = all_of(pval_col)) %>%
      filter(!is.na(pval)) %>%
      mutate(
        decile_val = decile_score(pval),
        weighted   = decile_val * scheme_w[scheme]
      ) %>%
      select(Gene, weighted)

    scheme_meta_df <- d %>%
      select(Gene, pval = all_of(pval_col)) %>%
      filter(!is.na(pval)) %>%
      mutate(
        decile_val = decile_score(pval),
        weighted   = decile_val * scheme_w[scheme]
      ) %>%
      select(
        Gene,
        !!paste0("accum_pval_", scheme) := pval,
        !!paste0("accum_score_", scheme) := decile_val,
        !!paste0("accum_weighted_", scheme) := weighted
      )

    if (is.null(accum_scores)) {
      accum_scores <- scheme_df %>% rename(!!scheme := weighted)
    } else {
      accum_scores <- accum_scores %>%
        full_join(scheme_df %>% rename(!!scheme := weighted), by = "Gene")
    }

    if (is.null(accum_meta)) {
      accum_meta <- scheme_meta_df
    } else {
      accum_meta <- accum_meta %>% full_join(scheme_meta_df, by = "Gene")
    }
  }

  if (!is.null(accum_scores)) {
    gene_rand <- accum_scores %>%
      rowwise() %>%
      mutate(gene_rand_score = sum(c_across(-Gene), na.rm = TRUE)) %>%
      ungroup() %>%
      select(Gene, gene_rand_score)
    if (!is.null(accum_meta)) {
      gene_rand <- gene_rand %>% left_join(accum_meta, by = "Gene")
    }
    cat(sprintf("  gene_rand_score: %d genes\n", nrow(gene_rand)))
  } else {
    gene_rand <- tibble(Gene = character(), gene_rand_score = numeric())
    has_accum <- FALSE
  }
} else {
  cat("Accumulation: not available, skipping\n")
  gene_rand <- tibble(Gene = character(), gene_rand_score = numeric())
}

# ── 3c. Gene RERConverge Score (optional) ───────────────────────────────────
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
  cat(sprintf("  Using %s for gene_rer_score\n", pval_col))

  gene_rer <- rer %>%
    filter(!is.na(.data[[pval_col]])) %>%
    mutate(gene_rer_score = decile_score(.data[[pval_col]]),
           rer_min_pval   = .data[[pval_col]]) %>%
    select(Gene = gene, gene_rer_score, rer_min_pval)

  cat(sprintf("  gene_rer_score: %d genes\n", nrow(gene_rer)))
} else {
  cat("RERConverge: not available, skipping\n")
  gene_rer <- tibble(Gene = character(), gene_rer_score = numeric())
}

# ── 3c. Gene FADE Score (optional) ──────────────────────────────────────────
# FADE is inherently gene-level (max BF per gene). Decile rank, higher = better.
if (has_fade && nrow(fade_combined) > 0) {
  gene_fade <- fade_combined %>%
    mutate(gene_fade_score = decile_score_higher_better(fade_max_bf)) %>%
    select(Gene, gene_fade_score, fade_max_bf)
  cat(sprintf("  gene_fade_score: %d genes\n", nrow(gene_fade)))
} else {
  gene_fade <- tibble(Gene = character(), gene_fade_score = numeric())
}

# ── 3e. gene_molerate_score — differential evolutionary rate (MoleRate, min p-value across top/bottom directions)
# MoleRate tests whether foreground branches evolve at a different rate than background.
# Decile rank of min(pval_pp_vs_prop) across directions, lower p-value → higher score.
if (has_molerate && nrow(mol_combined) > 0) {
  gene_molerate <- mol_combined %>%
    mutate(gene_molerate_score = decile_score(mol_min_pval)) %>%
    select(Gene, gene_molerate_score, mol_min_pval,
           any_of(c("mol_log2ratio_top", "mol_log2ratio_bottom")))
  cat(sprintf("  gene_molerate_score: %d genes\n", nrow(gene_molerate)))
} else {
  gene_molerate <- tibble(Gene = character(), gene_molerate_score = numeric(),
                          mol_log2ratio_top = numeric(), mol_log2ratio_bottom = numeric())
}

# ── 3f. Assemble gene scores ────────────────────────────────────────────────
gene_scores <- gene_caas

if (nrow(gene_rand) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_rand, by = "Gene")
} else {
  gene_scores$gene_rand_score <- NA_real_
}

if (nrow(gene_rer) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_rer, by = "Gene")
} else {
  gene_scores$gene_rer_score <- NA_real_
}

if (nrow(gene_fade) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_fade, by = "Gene")
} else {
  gene_scores$gene_fade_score <- NA_real_
}

if (nrow(gene_molerate) > 0) {
  gene_scores <- gene_scores %>% left_join(gene_molerate, by = "Gene")
} else {
  gene_scores$gene_molerate_score <- NA_real_
}

# Gene composite score: unweighted mean of available sub-scores
gene_scores <- gene_scores %>%
  rowwise() %>%
  mutate(
    gene_composite_score = {
      components <- c(gene_caas_score)
      if (!is.na(gene_rand_score))     components <- c(components, gene_rand_score)
      if (!is.na(gene_rer_score))      components <- c(components, gene_rer_score)
      if (!is.na(gene_fade_score))     components <- c(components, gene_fade_score)
      if (!is.na(gene_molerate_score)) components <- c(components, gene_molerate_score)
      mean(components, na.rm = TRUE)
    }
  ) %>%
  ungroup()

cat(sprintf("\nGene composite: min=%.3f, median=%.3f, max=%.3f\n",
            min(gene_scores$gene_composite_score, na.rm = TRUE),
            median(gene_scores$gene_composite_score, na.rm = TRUE),
            max(gene_scores$gene_composite_score, na.rm = TRUE)))

# =============================================================================
# 4. CORRELATION ANALYSIS (gene-level)
# =============================================================================

cat("\n─── Correlation analysis ──────────────────────────────────────\n")

score_cols <- c("gene_caas_score")
if (has_accum    && any(!is.na(gene_scores$gene_rand_score)))
  score_cols <- c(score_cols, "gene_rand_score")
if (has_rer      && any(!is.na(gene_scores$gene_rer_score)))
  score_cols <- c(score_cols, "gene_rer_score")
if (has_fade     && any(!is.na(gene_scores$gene_fade_score)))
  score_cols <- c(score_cols, "gene_fade_score")
if (has_molerate && any(!is.na(gene_scores$gene_molerate_score)))
  score_cols <- c(score_cols, "gene_molerate_score")

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
# 5. WRITE OUTPUTS
# =============================================================================

cat("\n─── Writing outputs ───────────────────────────────────────────\n")

# Position scores
pos_out <- pos_scores %>%
  select(Gene, Position, any_of(c("Pvalue", "pvalue_boot")), variability_score, pattern_score,
         scheme_coverage, biochem_score,
         asr_score, convergence_score, parallel_score,
         CAAS_score, change_side) %>%
  arrange(desc(CAAS_score))

write_tsv(pos_out, "position_scores.tsv")
cat(sprintf("  position_scores.tsv: %d rows\n", nrow(pos_out)))

# Gene scores
gene_out <- gene_scores %>%
  select(Gene, n_positions, gene_caas_score, gene_rand_score,
         gene_rer_score, gene_fade_score,
         gene_molerate_score,
         any_of(c("accum_pval_us", "accum_score_us", "accum_weighted_us",
                  "accum_pval_gs4", "accum_score_gs4", "accum_weighted_gs4",
                  "accum_pval_gs3", "accum_score_gs3", "accum_weighted_gs3",
                  "accum_pval_gs2", "accum_score_gs2", "accum_weighted_gs2",
                  "accum_pval_gs1", "accum_score_gs1", "accum_weighted_gs1")),
         any_of(c("rer_min_pval", "mol_min_pval", "fade_max_bf",
                  "mol_log2ratio_top", "mol_log2ratio_bottom")),
         gene_composite_score) %>%
  arrange(desc(gene_composite_score))

write_tsv(gene_out, "gene_scores.tsv")
cat(sprintf("  gene_scores.tsv: %d rows\n", nrow(gene_out)))

# Correlations
write_tsv(corr_results, "gene_correlations.tsv")

# =============================================================================
# 6. GENE LIST EXTRACTION FOR ORA
# =============================================================================

cat("\n─── Generating ORA gene lists ─────────────────────────────────\n")

# Helper: write gene list file, return count
write_gene_list <- function(genes, dir, filename) {
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) == 0) return(0)
  writeLines(genes, file.path(dir, filename))
  length(genes)
}

# ── Position-level ORA gene lists ────────────────────────────────────────────
# pos_scores is already filtered to the target direction; generate flat lists.
pos_for_ora <- pos_out %>%
  mutate(
    in_top10 = CAAS_score >= quantile(CAAS_score, 1 - top_pct,  na.rm = TRUE),
    in_top5  = CAAS_score >= quantile(CAAS_score, 1 - top5_pct, na.rm = TRUE),
    in_top1  = CAAS_score >= quantile(CAAS_score, 1 - top1_pct, na.rm = TRUE)
  )

for (pct_label in c("top10pct", "top5pct", "top1pct")) {
  col  <- paste0("in_", sub("pct$", "", pct_label))
  genes <- pos_for_ora %>% filter(.data[[col]]) %>% pull(Gene) %>% unique()
  fname <- sprintf("pos_caas_%s.txt", pct_label)
  n <- write_gene_list(genes, "pos_gene_lists", fname)
  cat(sprintf("  %s: %d genes\n", fname, n))
}

# ── Gene-level ORA gene lists ───────────────────────────────────────────────
# gene_out already reflects the direction-filtered position pool.
# Generate lists per score type and threshold — no internal direction loop needed.
gene_score_cols <- c("gene_caas_score", "gene_rand_score",
                     "gene_rer_score", "gene_fade_score", "gene_molerate_score",
                     "gene_composite_score")

for (score_col in gene_score_cols) {
  score_name <- sub("^gene_", "", sub("_score$", "", score_col))
  vals  <- gene_out[[score_col]]
  valid <- !is.na(vals)

  if (sum(valid) < 2) next

  thr_10 <- quantile(vals[valid], 1 - gene_top_pct,  na.rm = TRUE)
  thr_5  <- quantile(vals[valid], 1 - gene_top5_pct, na.rm = TRUE)
  thr_1  <- quantile(vals[valid], 1 - gene_top1_pct, na.rm = TRUE)

  for (pct_label in c("top10", "top5", "top1")) {
    thr   <- switch(pct_label, top10 = thr_10, top5 = thr_5, top1 = thr_1)
    genes <- gene_out %>%
      filter(!is.na(.data[[score_col]]) & .data[[score_col]] >= thr) %>%
      pull(Gene) %>% unique()
    fname <- sprintf("gene_%s_%spct.txt", score_name, pct_label)
    n <- write_gene_list(genes, "gene_gene_lists", fname)
    cat(sprintf("  %s: %d genes\n", fname, n))
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring — Complete\n")
cat("═══════════════════════════════════════════════════════════════\n")
