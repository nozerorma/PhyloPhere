#!/usr/bin/env Rscript
# =============================================================================
# CAAS Scoring — Position-Level and Gene-Level Composite Scores
# =============================================================================
#
# Integrates outputs from CT_POSTPROC, PGLS, FADE, RERConverge, and
# CT_ACCUMULATION into a unified scoring framework.
#
# Usage:
#   Rscript scoring_compute.R \
#     --postproc   <filtered_discovery.tsv> \
#     --pgls       <site_pgls.tsv> \
#     --fade_top   <fade_summary_top.tsv> \
#     --fade_bottom <fade_summary_bottom.tsv> \
#     --rer        <rerconverge_summary.tsv> \
#     --accum_dir  <directory_with_accumulation_CSVs> \
#     --top_pct    0.10 \
#     --top1_pct   0.01 \
#     --gene_top_pct  0.10 \
#     --gene_top1_pct 0.01
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

postproc_file     <- parse_arg("--postproc")
pgls_file         <- parse_arg("--pgls")
fade_top_file     <- parse_arg("--fade_top")
fade_bottom_file  <- parse_arg("--fade_bottom")
rer_file          <- parse_arg("--rer")
accum_dir         <- parse_arg("--accum_dir")
top_pct           <- as.numeric(parse_arg("--top_pct",  "0.10"))
top1_pct          <- as.numeric(parse_arg("--top1_pct", "0.01"))
gene_top_pct      <- as.numeric(parse_arg("--gene_top_pct",  "0.10"))
gene_top1_pct     <- as.numeric(parse_arg("--gene_top1_pct", "0.01"))

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

df <- df %>%
  mutate(scheme_weight = scheme_weights[CAAP_Group])

# Compute per-row scheme contribution
df <- df %>%
  mutate(scheme_contribution = scheme_weight * variability_score * pattern_score)

# Aggregate to Gene×Position level: sum scheme contributions (cumulative),
# and take the first value for fields that are constant within Gene×Position.
pos_scores <- df %>%
  group_by(Gene, Position) %>%
  summarise(
    biochem_score         = sum(scheme_contribution, na.rm = TRUE),
    # Take representative values from the row set (constant within Gene×Position)
    variability_score     = first(variability_score),
    pattern_score         = first(pattern_score),
    is_conserved_meta     = first(is_conserved_meta),
    asr_is_conserved      = first(asr_is_conserved),
    asr_root_conserved    = first(asr_root_conserved),
    conserved_pair        = first(conserved_pair),
    all_mrca_posterior    = first(all_mrca_posterior),
    mrca_1_posterior      = first(mrca_1_posterior),
    mrca_2_posterior      = first(mrca_2_posterior),
    mrca_3_posterior      = first(mrca_3_posterior),
    change_top            = first(change_top),
    change_bottom         = first(change_bottom),
    change_side           = first(change_side),
    parallel_top          = first(parallel_top),
    parallel_bottom       = first(parallel_bottom),
    .groups = "drop"
  )

cat(sprintf("  %d unique positions after aggregation\n", nrow(pos_scores)))

# ── 2c. ASR score ────────────────────────────────────────────────────────────
# Three cases:
#   is_conserved_meta = FALSE                         → 1.0 (novel change, earned)
#   is_conserved_meta = TRUE & asr_is_conserved = TRUE → 0.35 × mrca_pair_posterior
#                                                         + (if root: 0.35 × all_mrca_posterior)
#                                                         max = 0.70
#   is_conserved_meta = TRUE & asr_is_conserved = FALSE → 0.0
compute_asr_score <- function(is_conserved_meta, asr_is_conserved,
                              asr_root_conserved, conserved_pair,
                              mrca_1_posterior, mrca_2_posterior,
                              mrca_3_posterior, all_mrca_posterior) {

  parse_bool <- function(x) {
    if (is.logical(x)) return(x)
    tolower(as.character(x)) %in% c("true", "1", "yes")
  }

  is_cons   <- parse_bool(is_conserved_meta)
  asr_cons  <- parse_bool(asr_is_conserved)
  root_cons <- parse_bool(asr_root_conserved)

  # Look up the MRCA posterior for the conserved pair
  pair_idx  <- suppressWarnings(as.integer(conserved_pair))
  mrca_post <- case_when(
    pair_idx == 1 ~ as.numeric(mrca_1_posterior),
    pair_idx == 2 ~ as.numeric(mrca_2_posterior),
    pair_idx == 3 ~ as.numeric(mrca_3_posterior),
    TRUE          ~ as.numeric(mrca_1_posterior)  # fallback
  )

  case_when(
    !is_cons             ~ 1.0,
    is_cons & asr_cons   ~ pmin(
      0.35 * mrca_post +
        ifelse(root_cons, 0.35 * as.numeric(all_mrca_posterior), 0),
      0.70
    ),
    TRUE                 ~ 0.0
  )
}

pos_scores <- pos_scores %>%
  mutate(
    asr_score = compute_asr_score(
      is_conserved_meta, asr_is_conserved, asr_root_conserved,
      conserved_pair, mrca_1_posterior, mrca_2_posterior,
      mrca_3_posterior, all_mrca_posterior
    )
  )

# ── 2d. Convergence score ───────────────────────────────────────────────────
compute_convergence_score <- function(change_top, change_bottom) {
  side_value <- function(side) {
    case_when(
      side == "convergent"  ~ 1.0,
      side == "codivergent" ~ 0.7,
      side == "divergent"   ~ 0.3,
      TRUE                  ~ NA_real_  # ambiguous / no_change → excluded
    )
  }

  top_val    <- side_value(change_top)
  bottom_val <- side_value(change_bottom)

  # mean_change: mean of assessed sides (non-NA)
  mean_change <- ifelse(
    is.na(top_val) & is.na(bottom_val), 0,
    rowMeans(cbind(top_val, bottom_val), na.rm = TRUE)
  )

  # Bidirectional bonus
  change_both <- !is.na(top_val) & !is.na(bottom_val)
  bidirectional <- ifelse(change_both, 1.0, 0.7)

  (mean_change + bidirectional) / 2
}

pos_scores <- pos_scores %>%
  mutate(
    convergence_score = compute_convergence_score(change_top, change_bottom)
  )

# ── 2e. Parallel score ──────────────────────────────────────────────────────
compute_parallel_score <- function(parallel_top, parallel_bottom, change_top,
                                   change_bottom, mrca_1_posterior,
                                   mrca_2_posterior, mrca_3_posterior) {
  # Parallelism coefficient
  parallel_coeff <- function(p) {
    case_when(
      is.na(p) | p == "" | tolower(p) == "na" ~ NA_real_,
      tolower(p) == "non-parallel" | p == "False" | p == "FALSE" ~ 1.0,
      tolower(p) == "mixed"        ~ 0.7,
      tolower(p) == "parallel"     | p == "True" | p == "TRUE"  ~ 0.5,
      # Handle parallel_{direction} pattern from data
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

  # For each side, compute confidence-weighted score using MRCA posteriors.
  # Average across pairs that changed on that side.
  # Simplified: use mean of available MRCA posteriors for each side.
  mrca_mean <- rowMeans(
    cbind(
      as.numeric(mrca_1_posterior),
      as.numeric(mrca_2_posterior),
      as.numeric(mrca_3_posterior)
    ),
    na.rm = TRUE
  )

  top_side_score    <- ifelse(top_assessable,    top_coeff * mrca_mean, NA_real_)
  bottom_side_score <- ifelse(bottom_assessable, bottom_coeff * mrca_mean, NA_real_)

  mean_parallel <- ifelse(
    is.na(top_side_score) & is.na(bottom_side_score), 0,
    rowMeans(cbind(top_side_score, bottom_side_score), na.rm = TRUE)
  )

  # Bidirectional bonus (same as convergence)
  both_assessed <- !is.na(top_side_score) & !is.na(bottom_side_score)
  bidirectional <- ifelse(both_assessed, 1.0, 0.7)

  (mean_parallel + bidirectional) / 2
}

pos_scores <- pos_scores %>%
  mutate(
    parallel_score = compute_parallel_score(
      parallel_top, parallel_bottom, change_top, change_bottom,
      mrca_1_posterior, mrca_2_posterior, mrca_3_posterior
    )
  )

# ── 2f. PGLS score (optional) ───────────────────────────────────────────────
has_pgls <- file_exists(pgls_file)
if (has_pgls) {
  cat("Loading PGLS:", pgls_file, "\n")
  pgls <- read_tsv(pgls_file, show_col_types = FALSE) %>%
    filter(!is.na(p_pgls)) %>%
    mutate(pgls_score = decile_score(p_pgls)) %>%
    select(Gene, Position, pgls_score)

  pos_scores <- pos_scores %>%
    left_join(pgls, by = c("Gene", "Position"))

  cat(sprintf("  Merged %d PGLS scores\n", sum(!is.na(pos_scores$pgls_score))))
} else {
  cat("PGLS: not available, skipping\n")
  pos_scores$pgls_score <- NA_real_
}

# ── 2g. FADE (moved to gene-level — see section 3c) ─────────────────────────
# FADE operates at gene level (max Bayes Factor per gene across sites).
# It is loaded here so the data is available, but NOT broadcast to positions.
has_fade <- file_exists(fade_top_file) || file_exists(fade_bottom_file)
if (has_fade) {
  cat("Loading FADE summaries (gene-level; will be scored in section 3)\n")
  fade_dfs <- list()
  if (file_exists(fade_top_file)) {
    fade_dfs[["top"]] <- read_tsv(fade_top_file, show_col_types = FALSE) %>%
      select(gene, max_bf) %>%
      rename(Gene = gene)
  }
  if (file_exists(fade_bottom_file)) {
    fade_dfs[["bottom"]] <- read_tsv(fade_bottom_file, show_col_types = FALSE) %>%
      select(gene, max_bf) %>%
      rename(Gene = gene)
  }
  fade_combined <- bind_rows(fade_dfs) %>%
    group_by(Gene) %>%
    summarise(fade_max_bf = max(max_bf, na.rm = TRUE), .groups = "drop")
  cat(sprintf("  FADE: %d genes loaded\n", nrow(fade_combined)))
} else {
  cat("FADE: not available, skipping\n")
  fade_combined <- tibble(Gene = character(), fade_max_bf = numeric())
}

# ── 2h. CAAS Score (composite) ──────────────────────────────────────────────
# Unweighted mean of: biochem_score, asr_score, convergence_score,
# parallel_score, [pgls_score]
# FADE is gene-level and does NOT contribute to position-level CAAS_score.

pos_scores <- pos_scores %>%
  rowwise() %>%
  mutate(
    CAAS_score = {
      components <- c(biochem_score, asr_score, convergence_score, parallel_score)
      if (!is.na(pgls_score)) components <- c(components, pgls_score)
      mean(components, na.rm = TRUE)
    }
  ) %>%
  ungroup()

cat(sprintf("\nPosition-level CAAS_score: min=%.3f, median=%.3f, max=%.3f\n",
            min(pos_scores$CAAS_score, na.rm = TRUE),
            median(pos_scores$CAAS_score, na.rm = TRUE),
            max(pos_scores$CAAS_score, na.rm = TRUE)))

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

  for (scheme in scheme_names) {
    pattern <- paste0("accumulation_", scheme, "_aggregated_results.csv")
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

    if (is.null(accum_scores)) {
      accum_scores <- scheme_df %>% rename(!!scheme := weighted)
    } else {
      accum_scores <- accum_scores %>%
        full_join(scheme_df %>% rename(!!scheme := weighted), by = "Gene")
    }
  }

  if (!is.null(accum_scores)) {
    gene_rand <- accum_scores %>%
      rowwise() %>%
      mutate(gene_rand_score = sum(c_across(-Gene), na.rm = TRUE)) %>%
      ungroup() %>%
      select(Gene, gene_rand_score)
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
    mutate(gene_rer_score = decile_score(.data[[pval_col]])) %>%
    select(Gene = gene, gene_rer_score)

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
    select(Gene, gene_fade_score)
  cat(sprintf("  gene_fade_score: %d genes\n", nrow(gene_fade)))
} else {
  gene_fade <- tibble(Gene = character(), gene_fade_score = numeric())
}

# ── 3d. Assemble gene scores ────────────────────────────────────────────────
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

# Gene composite score: unweighted mean of available sub-scores
gene_scores <- gene_scores %>%
  rowwise() %>%
  mutate(
    gene_composite_score = {
      components <- c(gene_caas_score)
      if (!is.na(gene_rand_score))  components <- c(components, gene_rand_score)
      if (!is.na(gene_rer_score))   components <- c(components, gene_rer_score)
      if (!is.na(gene_fade_score))  components <- c(components, gene_fade_score)
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
if (has_accum && any(!is.na(gene_scores$gene_rand_score)))
  score_cols <- c(score_cols, "gene_rand_score")
if (has_rer && any(!is.na(gene_scores$gene_rer_score)))
  score_cols <- c(score_cols, "gene_rer_score")
if (has_fade && any(!is.na(gene_scores$gene_fade_score)))
  score_cols <- c(score_cols, "gene_fade_score")

if (length(score_cols) >= 2) {
  corr_mat <- gene_scores %>%
    select(all_of(score_cols)) %>%
    drop_na()

  corr_results <- expand.grid(
    score_a = score_cols, score_b = score_cols,
    stringsAsFactors = FALSE
  ) %>%
    filter(score_a < score_b) %>%
    rowwise() %>%
    mutate(
      pearson_r  = cor(corr_mat[[score_a]], corr_mat[[score_b]],
                       method = "pearson", use = "complete.obs"),
      spearman_r = cor(corr_mat[[score_a]], corr_mat[[score_b]],
                       method = "spearman", use = "complete.obs"),
      n_genes    = sum(complete.cases(corr_mat[[score_a]], corr_mat[[score_b]]))
    ) %>%
    ungroup()

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
  select(Gene, Position, variability_score, pattern_score, biochem_score,
         asr_score, convergence_score, parallel_score, pgls_score,
         CAAS_score, change_side) %>%
  arrange(desc(CAAS_score))

write_tsv(pos_out, "position_scores.tsv")
cat(sprintf("  position_scores.tsv: %d rows\n", nrow(pos_out)))

# Gene scores
gene_out <- gene_scores %>%
  select(Gene, n_positions, gene_caas_score, gene_rand_score,
         gene_rer_score, gene_fade_score, gene_composite_score) %>%
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
# Extract genes containing top-scoring positions
pos_for_ora <- pos_out %>%
  mutate(
    in_top10 = CAAS_score >= quantile(CAAS_score, 1 - top_pct,  na.rm = TRUE),
    in_top1  = CAAS_score >= quantile(CAAS_score, 1 - top1_pct, na.rm = TRUE)
  )

for (direction in c("top", "bottom", "both")) {
  subset_df <- if (direction == "both") {
    pos_for_ora
  } else {
    pos_for_ora %>% filter(change_side == direction)
  }

  # Top 10%
  genes10 <- subset_df %>% filter(in_top10) %>% pull(Gene) %>% unique()
  n <- write_gene_list(genes10, "pos_gene_lists",
                       sprintf("pos_caas_top10pct_%s.txt", direction))
  cat(sprintf("  pos_caas_top10pct_%s.txt: %d genes\n", direction, n))

  # Top 1%
  genes1 <- subset_df %>% filter(in_top1) %>% pull(Gene) %>% unique()
  n <- write_gene_list(genes1, "pos_gene_lists",
                       sprintf("pos_caas_top1pct_%s.txt", direction))
  cat(sprintf("  pos_caas_top1pct_%s.txt: %d genes\n", direction, n))
}

# ── Gene-level ORA gene lists ───────────────────────────────────────────────
# For gene-level, direction is determined by which side dominates in that gene.
# We assign each gene a dominant direction from its top-scoring position.
gene_direction <- pos_scores %>%
  group_by(Gene) %>%
  slice_max(CAAS_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Gene, gene_direction = change_side)

gene_for_ora <- gene_out %>%
  left_join(gene_direction, by = "Gene")

# Generate lists for each score type
gene_score_cols <- c("gene_caas_score", "gene_rand_score",
                     "gene_rer_score", "gene_fade_score", "gene_composite_score")

for (score_col in gene_score_cols) {
  if (all(is.na(gene_for_ora[[score_col]]))) next

  score_name <- sub("^gene_", "", sub("_score$", "", score_col))
  vals <- gene_for_ora[[score_col]]
  valid <- !is.na(vals)

  thr_10 <- quantile(vals[valid], 1 - gene_top_pct,  na.rm = TRUE)
  thr_1  <- quantile(vals[valid], 1 - gene_top1_pct, na.rm = TRUE)

  gene_for_ora <- gene_for_ora %>%
    mutate(
      !!paste0(score_col, "_top10") := ifelse(!is.na(.data[[score_col]]) & .data[[score_col]] >= thr_10, TRUE, FALSE),
      !!paste0(score_col, "_top1")  := ifelse(!is.na(.data[[score_col]]) & .data[[score_col]] >= thr_1,  TRUE, FALSE)
    )

  for (direction in c("top", "bottom", "both")) {
    subset_df <- if (direction == "both") {
      gene_for_ora
    } else {
      gene_for_ora %>% filter(gene_direction == direction)
    }

    for (pct_label in c("top10", "top1")) {
      flag_col <- paste0(score_col, "_", pct_label)
      genes <- subset_df %>% filter(.data[[flag_col]]) %>% pull(Gene) %>% unique()
      fname <- sprintf("gene_%s_%spct_%s.txt", score_name, sub("top", "top", pct_label), direction)
      n <- write_gene_list(genes, "gene_gene_lists", fname)
      cat(sprintf("  %s: %d genes\n", fname, n))
    }
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring — Complete\n")
cat("═══════════════════════════════════════════════════════════════\n")
