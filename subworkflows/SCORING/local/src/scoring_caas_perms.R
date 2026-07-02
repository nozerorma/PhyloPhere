#!/usr/bin/env Rscript
# =============================================================================
# scoring_caas_perms.R — CAAS permulation-excess null → genes×N matrices
# =============================================================================
# Turns the long per-(gene, cycle, position, scheme) asr_path_score table emitted
# by disambiguation_perms_main.py (caas_perm_scores.tsv) into three genes×N null
# matrices (global / top / bottom) and writes caas_perms.rds in the shape the FCS
# permulation path consumes:
#
#   saveRDS(list(corStat_byrank = list(global_asr = M, top_asr = M, bottom_asr = M)))
#
# The aggregation MIRRORS scoring_compute.R EXACTLY so the null shares the observed
# statistic (the calibration requirement for a valid empirical p.perm):
#   position asr_score      = weighted.mean(asr_path_score, w = scheme_weight)
#   gene×cycle null score   = quantile(asr_score, 0.90) over the gene's positions
#   directional subsets      = change_side ∈ {top,both} / {bottom,both}
# Ranking names match fcs_stats' score_*_asr columns: global_asr/top_asr/bottom_asr.
#
# NULL SEMANTICS CAVEAT — this is a NAIVE random-relabeling null. The permuted
# labelings (from CAAStools resample) do NOT preserve the observed matched-sister-
# pair design: permuted fg/bg species are paired at RANDOM across the tree, so pair
# MRCAs are deep and the asr_path_score isolation/independence geometry differs from
# the observed (which uses shallow matched pairs). We can't preserve matched pairs
# under permutation — a random labeling almost never yields N independent matched
# contrasts — so random pairing is a necessity, and the up-side is a homogeneous,
# well-defined null. Consequence: p.perm tests "phenotype + matched-pair design"
# jointly (mildly anti-conservative on the MRCA axes), not the phenotype alone. A
# within-pair fg/bg swap would isolate the phenotype but spans only 2^n labelings
# (too coarse unless n is large). POSENRICH therefore drops p.perm entirely and
# relies on the exact XL-mHG p.adj; the gene FCS keeps p.perm with this caveat.
#
# Usage:
#   Rscript scoring_caas_perms.R --perm-scores caas_perm_scores.tsv \
#       [--universe cleaned_background.txt] --output caas_perms.rds
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

# ── minimal flag parser ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
perm_scores_file <- get_arg("--perm-scores")
universe_file    <- get_arg("--universe", "NO_FILE")
out_file         <- get_arg("--output", "caas_perms.rds")
quant            <- as.numeric(get_arg("--quantile", "0.90"))

stopifnot(!is.null(perm_scores_file), file.exists(perm_scores_file))

# Same scheme weights + change classification as scoring_compute.R.
scheme_weights_int <- c(US = 0.2, GS4 = 0.2, GS3 = 0.2, GS2 = 0.2, GS1 = 0.2)
assessable_change <- function(x) x %in% c("convergent", "codivergent", "divergent")

ps <- read_tsv(perm_scores_file, show_col_types = FALSE)
need <- c("Gene", "cycle", "Position", "caap_group", "asr_path_score",
          "change_top", "change_bottom")
miss <- setdiff(need, names(ps))
if (length(miss)) stop("perm-scores missing columns: ", paste(miss, collapse = ", "))

ps <- ps %>%
  mutate(
    asr_path_score = suppressWarnings(as.numeric(asr_path_score)),
    scheme_weight  = scheme_weights_int[caap_group]
  ) %>%
  filter(!is.na(scheme_weight))

if (nrow(ps) == 0) {
  message("[caas_perms] no scored rows; writing empty RDS")
  saveRDS(list(corStat_byrank = list()), out_file)
  quit(status = 0)
}

# ── Position-level aggregation (scheme-weighted), per gene×cycle×position ──────
pos <- ps %>%
  group_by(Gene, cycle, Position) %>%
  summarise(
    asr_score         = weighted.mean(asr_path_score, w = scheme_weight, na.rm = TRUE),
    has_change_top    = any(assessable_change(change_top)),
    has_change_bottom = any(assessable_change(change_bottom)),
    .groups = "drop"
  ) %>%
  mutate(
    change_side = case_when(
      has_change_top & has_change_bottom ~ "both",
      has_change_top                     ~ "top",
      has_change_bottom                  ~ "bottom",
      TRUE                               ~ "none"
    )
  )

# ── Gene×cycle null scores per direction (0.90-quantile of asr_score) ──────────
q90 <- function(v) if (length(v) > 0) as.numeric(quantile(v, quant, na.rm = TRUE)) else NA_real_
gene_cycle <- pos %>%
  group_by(Gene, cycle) %>%
  summarise(
    global_asr = q90(asr_score),
    top_asr    = q90(asr_score[change_side %in% c("top", "both")]),
    bottom_asr = q90(asr_score[change_side %in% c("bottom", "both")]),
    .groups = "drop"
  )

cycle_levels <- sort(unique(gene_cycle$cycle))
n_perms <- length(cycle_levels)

# ── Recovery (detection-frequency) p-value per (gene, position, scheme) ────────
# k = #cycles a (gene,position,scheme) is detected; N = total perm cycles.
# recovery_pval = (N - k + 1)/(N + 1) — the SAME (1+b)/(N+1) form used for every
# p.perm in the pipeline (b = N - k = cycles NOT detected); floor 1/(N+1) when
# detected in every cycle. This is the null analogue of the observed pvalue_boot
# (a bootstrap recovery frequency), so pheno = 1 - recovery_pval mirrors the
# observed pheno = 1 - pvalue_boot used to build CAAS_score.
n_cycles_total <- length(unique(ps$cycle))
recovery <- ps %>%
  distinct(Gene, Position, caap_group, cycle) %>%
  count(Gene, Position, caap_group, name = "n_detected") %>%
  mutate(n_cycles      = n_cycles_total,
         recovery_pval = (n_cycles_total - n_detected + 1) / (n_cycles_total + 1))

# ── Null CAAS score — MIRRORS scoring_compute.R per scheme, then combined ──────
#   row_caas = (1 - recovery_pval) * asr_path_score * scheme_weight   (per scheme)
#   caas_score(position) = sum(row_caas across schemes)               (per gene×cycle×pos)
ps_caas <- ps %>%
  left_join(recovery, by = c("Gene", "Position", "caap_group")) %>%
  mutate(row_caas = (1 - recovery_pval) * asr_path_score * scheme_weight)

pos_caas <- ps_caas %>%
  group_by(Gene, cycle, Position) %>%
  summarise(
    caas_score        = sum(row_caas, na.rm = TRUE),
    has_change_top    = any(assessable_change(change_top)),
    has_change_bottom = any(assessable_change(change_bottom)),
    .groups = "drop"
  ) %>%
  mutate(change_side = case_when(
    has_change_top & has_change_bottom ~ "both",
    has_change_top                     ~ "top",
    has_change_bottom                  ~ "bottom",
    TRUE                               ~ "none"))

gene_cycle_caas <- pos_caas %>%
  group_by(Gene, cycle) %>%
  summarise(
    global_caas = q90(caas_score),
    top_caas    = q90(caas_score[change_side %in% c("top", "both")]),
    bottom_caas = q90(caas_score[change_side %in% c("bottom", "both")]),
    .groups = "drop"
  )

# Gene universe: cleaned_background if given (so the null shares the observed
# zero-floored universe), else the genes present in the null table.
perm_genes <- sort(unique(gene_cycle$Gene))
universe <- perm_genes
if (!is.null(universe_file) && universe_file != "NO_FILE" && file.exists(universe_file)) {
  u <- tryCatch(readLines(universe_file), error = function(e) character(0))
  u <- trimws(u); u <- u[nzchar(u) & u != "Gene"]
  if (length(u)) universe <- sort(unique(c(u, perm_genes)))
}

# ── Build a genes×N matrix for one direction (absent gene/cycle → 0 = no signal)
build_matrix <- function(df, col) {
  m <- df %>%
    select(Gene, cycle, !!col := dplyr::all_of(col)) %>%
    tidyr::pivot_wider(names_from = cycle, values_from = !!col)
  # align rows to the universe, columns to all cycles
  mat <- matrix(0.0, nrow = length(universe), ncol = n_perms,
                dimnames = list(universe, cycle_levels))
  if (nrow(m) > 0) {
    rn <- m$Gene
    sub <- as.matrix(m[, setdiff(names(m), "Gene"), drop = FALSE])
    sub[is.na(sub)] <- 0.0
    common_cols <- intersect(colnames(sub), cycle_levels)
    common_rows <- intersect(rn, universe)
    if (length(common_rows) && length(common_cols)) {
      ri <- match(common_rows, rn)
      mat[common_rows, common_cols] <- sub[ri, common_cols, drop = FALSE]
    }
  }
  mat
}

# corStat_byrank (asr axis) is the FCS p.perm null — UNCHANGED so the FCS path is
# untouched. caas_corStat_byrank (CAAS_score axis) is NEW, for the report's gene-level
# CAAS-score null overlay only.
corStat_byrank <- list(
  global_asr = build_matrix(gene_cycle, "global_asr"),
  top_asr    = build_matrix(gene_cycle, "top_asr"),
  bottom_asr = build_matrix(gene_cycle, "bottom_asr")
)
caas_corStat_byrank <- list(
  global = build_matrix(gene_cycle_caas, "global_caas"),
  top    = build_matrix(gene_cycle_caas, "top_caas"),
  bottom = build_matrix(gene_cycle_caas, "bottom_caas")
)

saveRDS(list(corStat_byrank = corStat_byrank,
             caas_corStat_byrank = caas_corStat_byrank), out_file)
cat(sprintf("[caas_perms] wrote %s — %d genes × %d cycles, asr+caas nulls (%s)\n",
            out_file, length(universe), n_perms,
            paste(names(corStat_byrank), collapse = ", ")))

# ── Lean report inputs (so the report never loads the raw per-cycle table) ─────
out_dir <- dirname(out_file)
# (1) recovery p-value per (gene, position, scheme) — small, full fidelity.
readr::write_tsv(
  recovery %>% select(Gene, Position, caap_group, n_detected, n_cycles, recovery_pval),
  file.path(out_dir, "perm_pos_pval.tsv"))

# (2) capped per-scheme sample of the per-(position,cycle) null for distribution
# plots (asr_path_score + per-scheme CAAS contribution row_caas). Full data scales
# to ~genes×positions×schemes×N; cap keeps the report light at 1K cycles.
set.seed(1)
cap <- 20000L
sample_df <- ps_caas %>%
  select(Gene, Position, caap_group, cycle, asr_path_score, recovery_pval, row_caas) %>%
  group_by(caap_group) %>%
  slice_sample(n = cap) %>%   # dplyr caps to group size when a group has < cap rows
  ungroup()
readr::write_tsv(sample_df, file.path(out_dir, "perm_pos_sample.tsv"))
cat(sprintf("[caas_perms] wrote perm_pos_pval.tsv (%d rows) + perm_pos_sample.tsv (%d rows)\n",
            nrow(recovery), nrow(sample_df)))
