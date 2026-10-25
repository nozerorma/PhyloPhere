#!/usr/bin/env Rscript
# scoring_compute.R — Position-level and gene-level CAAS scores, integrated with FADE, RER and accumulation, with the
# biochemical groupings that the derived side of each position agrees with.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  SCORING_COMPUTE Nextflow process (scoring_compute.nf → Rscript scoring_compute.R ...),
#             after observed_core_scores.py, whose two tables it reads
#
# The position score (CAAS_score) and the gene score (raw maximum) come from core.scores via observed_core_scores.py.
# The biochemical groupings (aa_grouping.R) do not enter either score: they annotate each position with the groupings
# in which its derived residues share one class that differs from the ancestral one. This script adds the empirical permulation p-values of positions (p.emp and
# their BH adjustments), joins the per-gene evidence of
# FADE, RERConverge and accumulation, and writes the tables that the ENRICHMENT
# subworkflow and the reports read (the nominal activity lists are produced by
# scoring_evidence.R). It runs once on the full pool of
# filtered_discovery.tsv; direction is carried by the `side` column (top, bottom, none).
#
# Args (named flags, from task.script). An input whose file name starts with NO_ counts as absent.
#   --postproc            filtered_discovery.tsv (mandatory)
#   --core_positions      observed_core_scores.py positions table (mandatory)
#   --core_genes          observed_core_scores.py genes table (mandatory)
#   --fade_top, --fade_bottom         gene-level FADE summaries (fade_summary_{top,bottom}.tsv)
#   --fade_site_top, --fade_site_bot  per-site FADE tables (comma-delimited)
#   --rer                 rerconverge_summary_<trait>.tsv
#   --accum_dir           directory with accumulation_<direction>_results.csv
#   --null_counters       null_counters.tsv (null_position_stats.py reduce): per observed position the counters of the
#                         position-level null used by p.emp
#   --null_fit            null_fit.tsv: legacy staging metadata (unused)
#   --null_base           null_universe_base.tsv: what the null fixes of the universe summary (sizes, saturation curve)
#   --cycle_roster        cycle_roster.txt of the CAAS_CORE_MERGE run: N of the null, the cycles with no detection included
#   --hypotheses_pairs    contrast_hypotheses_pairs.tsv (FOP); parsed but not used by the computation
#
# Outputs (working directory):
#   position_scores.tsv            one row per Gene, Position and side with a positive score; biochem_groupings and
#                                  is_biochem give the groupings the derived side agrees with; nd (cycles of the null
#                                  that score the position) and null_unseen (nd == 0) next to the position p-values
#   null_universe_summary.tsv      N, size of the position universe, observed positions the null never scored and the
#                                  saturation curve of the null (written when a null is given)
#   gene_scores.tsv                one row per Gene
#   gene_correlations.tsv          pairwise correlations between gene scores (header only while one score exists)
#   fcs_stats.tsv                  gene scores and flags for the FCS reports
#   fcs_stats_{rer,fade,accum}.tsv per-module FCS rankings, written when that evidence exists
#                                  (the nominal activity lists gene_lists/ and position_lists/ are written by
#                                  scoring_evidence.R, not by this script)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})
source("aa_grouping.R")   # derived_groupings(), staged next to this script by scoring_compute.nf
source("null_position_p.R")   # null_position_p(), null_universe_summary(), staged next to this script by scoring_compute.nf

# ── Parse arguments ───────────────────────────────────────────────────────────
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
core_positions_file  <- parse_arg("--core_positions")  # observed_core_scores.py: CAAS_score per (Gene, Position, side)
core_genes_file      <- parse_arg("--core_genes")      # observed_core_scores.py: size-adjusted gene CAAS scores
null_counters_file   <- parse_arg("--null_counters")     # null_counters.tsv (counters of the position-level null); NO_FILE otherwise
null_fit_file        <- parse_arg("--null_fit")          # legacy staging metadata (unused); NO_FILE otherwise
null_base_file       <- parse_arg("--null_base")         # null_universe_base.tsv (sizes and saturation curve of the null); NO_FILE otherwise
cycle_roster_file    <- parse_arg("--cycle_roster")      # cycle_roster.txt: every replayed cycle (N of the null); NO_FILE otherwise
# Rule of the position score: the CAAS score (mean over the domains of the PSS-weighted mean over their species), defined only when positive (core.scores.SCORE_RULE). The observed scores (observed_core_scores.py)
# and the null (counts.json of the null, score_rule; checked by null_position_stats.py) must use the same rule.
score_aggregation <- "caas_pos_species_pss"
# Rows reach this script already pooled over hypotheses by CT_DISAMBIGUATION (one row per Gene,
# Position and side, hypothesis NA), so scoring never pools hypotheses itself.
nominal_threshold <- as.numeric(parse_arg("--nominal_threshold", "0.05"))
fdr_threshold <- as.numeric(parse_arg("--fdr_threshold", "0.05"))
accumulation_fdr_threshold <- as.numeric(parse_arg("--accumulation_fdr_threshold", "0.1"))
fade_bf_threshold <- as.numeric(parse_arg("--fade_bf_threshold", "100"))
if (!is.finite(accumulation_fdr_threshold) || accumulation_fdr_threshold < 0 || accumulation_fdr_threshold > 1) {
  stop("--accumulation_fdr_threshold must be a finite number in [0, 1]")
}
# Direction is not a parameter: scoring runs on the full pool and the `side` column carries it.

# TRUE for a real file; a sentinel (basename starting with NO_) or an empty name counts as absent.
file_exists <- function(f) {
  !is.null(f) && f != "" && !grepl("^NO_", basename(f)) && file.exists(f)
}

# Position and gene CAAS scores are computed once, by core.scores (observed_core_scores.py for
# this run, gene_wrapper.py for the permulation null), and read back here. This script
# integrates them with FADE / RER / accumulation and tests the observed position score
# against the null.
#
# Ties: position scores are means of a few values, so scores that are equal in exact
# arithmetic can differ by rounding noise. Comparisons against the null count values within
# TIE_TOL as ties (same constant as core.scores.TIE_TOL).
TIE_TOL <- 1e-12

# The position-level permulation p (p.emp, its BH family and the universe summary) is computed from the
# counters of the null by null_position_p.R.

# Correlation over the complete pairs; NA when fewer than 3 remain.
safe_cor <- function(x, y, method = "pearson") {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CAAS Scoring - Compute\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── 1. Load postproc data (mandatory) ─────────────────────────────────────────
stopifnot(file_exists(postproc_file))
stopifnot("--core_positions and --core_genes (observed_core_scores.py) are required" =
            file_exists(core_positions_file) && file_exists(core_genes_file))
cat("Loading postproc:", postproc_file, "\n")
df <- read_tsv(postproc_file, show_col_types = FALSE)
# filtered_discovery.tsv carries the canonical lowercase column names of CT_DISAMBIGUATION
# (asr_path_score, side, ...), used as they are; position_scores.tsv keeps the same schema.
cat(sprintf("  %d rows, %d unique Gene×Position pairs\n",
            nrow(df), n_distinct(paste(df$Gene, df$Position))))

# ── 2. Position-level scoring ─────────────────────────────────────────────────

# ── 2b. Hypothesis pooling (already done upstream) ────────────────────────────
# Rows arrive pooled over hypotheses by CT_DISAMBIGUATION (fop_pool.pool_domains). H1..Hn are
# overlapping K-pair designs over the same Voronoi domains, not independent replicates: averaging over
# them in section 2g would dilute a strong canonical signal and let a position with many harvested
# hypotheses distort every genome-wide rank. The pooling averages the per-hypothesis domain scores over
# the K fixed domains, weighted by the PSS of contrast_hypotheses_pairs.tsv; a single contrast reduces
# to the plain PSS-weighted domain mean.
cat("  FOP pooling: done in-tree per (Gene, Position, side) [core v3]\n")
# Backfill the descriptor columns that section 2g and the reports expect, so that a missing column is
# never hit.
if (!"n_hypotheses" %in% names(df))          df$n_hypotheses <- 1L
if (!"participating_hypotheses" %in% names(df)) df$participating_hypotheses <- ""
# Descriptors of the species that carry the change; empty text when the input has none.
for (.c in c("top_species_residues", "bottom_species_residues",
             "n_top_species", "n_bottom_species", "n_conserved_pairs")) {
  if (!.c %in% names(df)) df[[.c]] <- ""
}
df$asr_path_score <- suppressWarnings(as.numeric(df$asr_path_score))


# Detect the posterior columns of the K fixed Voronoi domains (domain_<d>_posterior; the
# mrca_<i>_posterior spelling is accepted too). Only their number is reported.
mrca_posterior_cols <- grep("^(?:mrca|domain)_\\d+_posterior$", names(df),
                            value = TRUE, perl = TRUE)
n_pairs <- length(mrca_posterior_cols)
cat(sprintf("  Detected %d domains (%s)\n", n_pairs, paste(mrca_posterior_cols, collapse = ", ")))

# ── 2c. ASR score (per-row) ───────────────────────────────────────────────────
# asr_score is the unified ASR path score computed upstream in CT_DISAMBIGUATION
# (src/convergence/path_scores.py).
stopifnot("asr_path_score" %in% names(df))
cat("  Using upstream asr_path_score (unified ASR/convergence/parallel signal)\n")
df <- df %>%
  mutate(
    asr_score = suppressWarnings(as.numeric(asr_path_score))
  )

# Keep the diagnostic column derived_agreement present even when the input lacks it.
if (!"derived_agreement" %in% names(df)) df$derived_agreement <- NA_real_
df$derived_agreement <- suppressWarnings(as.numeric(df$derived_agreement))


# ── 2f. Per-row CAAS score ────────────────────────────────────────────────────
# A row's CAAS score is its asr_path_score (the unified ASR path score computed upstream);
# the position score aggregating rows is computed by core.scores (section 2g).

# ── 2g. Aggregate to Gene×Position ────────────────────────────────────────────
# Position as integer: §2f-ter joins pos_scores to the per-cycle null on (Gene, Position), and the
# null reads Position as integer, so the key types must agree.
df <- df %>% mutate(Position = suppressWarnings(as.integer(Position)))
# `side` belongs to the aggregation key: a position detected on both sides has two rows, and
# CAAS_score is taken per side.
.pos_grp_keys <- c("Gene", "Position", "side")

pos_scores <- df %>%
  group_by(across(all_of(.pos_grp_keys))) %>%
  summarise(
    CAAS_score         = NA_real_,  # filled from the core table below
    # FOP descriptors: the hypotheses were pooled upstream, so these are per-position columns, not a
    # count over rows. They never enter CAAS_score.
    n_hypotheses       = if ("n_hypotheses" %in% names(df)) {
      .nh <- n_hypotheses[is.finite(n_hypotheses)]; if (length(.nh)) max(.nh) else 0L
    } else 0L,
    participating_hypotheses = if ("participating_hypotheses" %in% names(df)) {
      .ph <- unique(participating_hypotheses[!is.na(participating_hypotheses) & nzchar(participating_hypotheses)])
      if (length(.ph)) paste(sort(unique(unlist(strsplit(.ph, ",")))), collapse = ",") else ""
    } else "",
    # Position-level descriptors: first() carries them, falling back to caas pattern when empty.
    caas                    = if ("caas" %in% names(df)) dplyr::first(caas) else "",
    top_species_residues    = {
      .tsr <- if ("top_species_residues" %in% names(df)) dplyr::first(top_species_residues) else ""
      .c   <- if ("caas" %in% names(df)) dplyr::first(caas) else ""
      if (is.na(.tsr) || !nzchar(.tsr) || identical(.tsr, "NA")) ifelse(grepl("/", .c), sub("/.*", "", .c), .c) else .tsr
    },
    bottom_species_residues = {
      .bsr <- if ("bottom_species_residues" %in% names(df)) dplyr::first(bottom_species_residues) else ""
      .c   <- if ("caas" %in% names(df)) dplyr::first(caas) else ""
      if (is.na(.bsr) || !nzchar(.bsr) || identical(.bsr, "NA")) ifelse(grepl("/", .c), sub(".*/", "", .c), "") else .bsr
    },
    n_top_species           = if ("n_top_species" %in% names(df)) dplyr::first(n_top_species) else "",
    n_bottom_species        = if ("n_bottom_species" %in% names(df)) dplyr::first(n_bottom_species) else "",
    n_conserved_pairs       = if ("n_conserved_pairs" %in% names(df)) dplyr::first(n_conserved_pairs) else "",
    .groups = "drop"
  )

# CAAS_score = asr_path_score of the (Gene, Position, side) row, computed by core.scores.
core_pos <- read_tsv(core_positions_file, show_col_types = FALSE,
                     col_types = cols(Gene = col_character(), Position = col_integer(),
                                      side = col_character(), CAAS_score = col_double()))
.pk  <- function(g, p, sd) paste(g, p, sd, sep = "\r")
.hit <- match(.pk(pos_scores$Gene, pos_scores$Position, pos_scores$side),
              .pk(core_pos$Gene, core_pos$Position, core_pos$side))
if (anyNA(.hit) || nrow(core_pos) != nrow(pos_scores)) {
  stop(sprintf("core positions (%d) and detected positions (%d) disagree: observed_core_scores.py and this script read different rows",
               nrow(core_pos), nrow(pos_scores)))
}
pos_scores$CAAS_score <- core_pos$CAAS_score[.hit]
rm(core_pos, .hit)
cat(sprintf("  %d detected position sides\n", nrow(pos_scores)))

# Keep only positions with true convergence (CAAS_score > 0 and side top/bottom).
# Rows with CAAS_score <= 0 or NA (shadow single-domain rows, divergent non-identical
# residues, or side "none") are filtered from position_scores.tsv; they stay in filtered_discovery.tsv.
pos_scores <- pos_scores %>%
  filter(!is.na(CAAS_score) & CAAS_score > 0 & side %in% c("top", "bottom"))

# CAAS_score is the score of that side; asr_path_score repeats it.
pos_scores <- pos_scores %>% mutate(asr_path_score = CAAS_score)

# Ancestral and derived residues of each (Gene, Position, side), as the ASR inferred them. They are read from
# the K fixed domains of the row: domain_<d>_anc_aa is the ancestral residue of domain d and
# domain_<d>_top_aa / domain_<d>_bot_aa the derived one on the side that carries the change. Residues are pooled
# over the domains, most frequent first, comma-separated. Empty text when the position has no resolved domain,
# or side "none" for the derived set. The `caas` pattern (top/bottom residues) has no direction.
.anc_cols <- sort(grep("^domain_\\d+_anc_aa$", names(df), value = TRUE))
.pool_aa <- function(m) {
  vapply(seq_len(nrow(m)), function(i) {
    v <- m[i, ]
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) return("")
    n <- table(v)
    paste(names(n)[order(-as.integer(n), names(n))], collapse = ",")
  }, character(1))
}
.rep <- df %>%
  filter(!duplicated(paste(Gene, Position, side, sep = "\r"))) %>%
  select(Gene, Position, side, any_of(c(.anc_cols, sub("_anc_aa$", "_top_aa", .anc_cols), sub("_anc_aa$", "_bot_aa", .anc_cols))))
if (length(.anc_cols) && nrow(.rep)) {
  .mat <- function(cols) as.matrix(.rep[, cols, drop = FALSE])
  .der <- .mat(sub("_anc_aa$", "_top_aa", .anc_cols))
  .bot <- .mat(sub("_anc_aa$", "_bot_aa", .anc_cols))
  .der[.rep$side == "bottom", ] <- .bot[.rep$side == "bottom", ]
  .der[!.rep$side %in% c("top", "bottom"), ] <- NA_character_
  .rep <- .rep %>% transmute(Gene, Position, side, ancestral_aa = .pool_aa(.mat(.anc_cols)), derived_aa = .pool_aa(.der))
  pos_scores <- pos_scores %>% left_join(.rep, by = c("Gene", "Position", "side"))
  pos_scores$ancestral_aa[is.na(pos_scores$ancestral_aa)] <- ""
  pos_scores$derived_aa[is.na(pos_scores$derived_aa)] <- ""
} else {
  pos_scores$ancestral_aa <- ""
  pos_scores$derived_aa <- ""
}
rm(.rep, .anc_cols)

# Groupings (aa_grouping.R) in which the derived residues of the position share one class that is not an ancestral class.
pos_scores$biochem_groupings <- derived_groupings(pos_scores$ancestral_aa, pos_scores$derived_aa, pos_scores$caas, pos_scores$side)
pos_scores$is_biochem <- nzchar(pos_scores$biochem_groupings)

cat(sprintf("  %d unique positions after aggregation\n", nrow(pos_scores)))

cat(sprintf("\nPosition-level CAAS_score: min=%.3f, median=%.3f, max=%.3f\n",
            min(pos_scores$CAAS_score, na.rm = TRUE),
            median(pos_scores$CAAS_score, na.rm = TRUE),
            max(pos_scores$CAAS_score, na.rm = TRUE)))

# ── 2f-ter. p.emp: position-level "detects AND exceeds" permulation p ─────────
# Design: docs/scoring_v2_p_emp.md §1-§2. Per (Gene, Position):
#   k_emp = #{null cycle : re-detects the position on ANY side
#                          AND  max_side(caas_score) >= max_side(CAAS_obs)}
#   p.emp = (k_emp + 1) / (N + 1)                          add-one, right-tailed
# The pooled statistic is the max-over-sides "all" axis -- identical to .pos_undirected (§4a) on the observed side and
# to the position score of the null. Values within TIE_TOL of the observed score count as ties (>=). The null is not read
# as a table of every (position, cycle): null_position_stats.py sums, over the blocks of cycles of the null, the integers
# used by p.emp (nd and k_emp); null_position_p.R turns them into the p and
# the BH families (see its header).
pos_scores$p.emp <- NA_real_
pos_scores$nd <- NA_integer_          # cycles of the null that score the position (informational)
pos_scores$null_unseen <- NA          # nd == 0: the null never produced the position (informational, no demotion)
pos_scores$p.adj_bh  <- NA_real_
has_null_counters <- file_exists(null_counters_file) && file_exists(null_base_file) &&
  file_exists(cycle_roster_file)
# A null in which no cycle has a row (N = 0, or no permuted cycle re-detected any position) has no cycle to count:
# (k + 1) / (N + 1) would read 1 for every position, a value rather than the absence of one. p.emp and p.adj_bh stay NA,
# as when no null is given.
if (has_null_counters) {
  .null_base <- read_tsv(null_base_file, show_col_types = FALSE, col_types = "cddc") %>% mutate(cycles = suppressWarnings(as.integer(cycles)))
  if (!isTRUE(.null_base$value[.null_base$metric == "N_cycles_present_in_null" & is.na(.null_base$cycles)][1] > 0)) {
    cat("  the null has no row: no null cycle, so p.emp and p.adj_bh stay NA\n")
    has_null_counters <- FALSE
  }
}
if (has_null_counters) {
  cat("Loading the counters of the position-level null (p.emp):", null_counters_file, "\n")
  # N = cycles replayed, not cycles with a detection: a cycle that re-detects no position counts in the denominator.
  .roster <- trimws(readLines(cycle_roster_file, warn = FALSE)); .roster <- .roster[nzchar(.roster)]
  if (anyDuplicated(.roster)) stop("Duplicate null cycle IDs in cycle_roster.txt")
  N_emp <- length(.roster)
  cat(sprintf("  N = %d cycles (source: cycle_roster.txt)\n", N_emp))
  .counters <- read_tsv(null_counters_file, show_col_types = FALSE, col_types = cols(Gene = col_character(), Position = col_integer(), .default = col_integer()))
  .fit <- NULL  # legacy staging input is no longer used

  # observed pooled statistic = best side per (Gene, Position) (== .pos_undirected)
  obs_max <- pos_scores %>%
    filter(!is.na(CAAS_score)) %>%
    group_by(Gene, Position) %>%
    summarise(.obs = max(CAAS_score), .groups = "drop")

  .n_u <- .null_base$value[.null_base$metric == "null_positive_positions_in_universe" & is.na(.null_base$cycles)][1]
  .universe_note <- .null_base$note[.null_base$metric == "null_positive_positions_in_universe" & is.na(.null_base$cycles)][1]
  cat(sprintf("  universe: %s; %d null-positive positions in U\n", .universe_note, as.integer(.n_u)))
  .res <- null_position_p(obs_max, .counters, .fit, .n_u, N_emp)

  pos_scores <- pos_scores %>%
    select(-p.emp, -nd, -null_unseen, -p.adj_bh) %>%
    left_join(.res$pos, by = c("Gene", "Position"))
  cat(sprintf("  null adequacy: %d of %d observed positive positions never scored by the null (%.1f%%); kept at p = 1/(N+1)\n",
              .res$n_unseen, .res$n_obs, if (.res$n_obs > 0) 100 * .res$n_unseen / .res$n_obs else 0))
  write_tsv(null_universe_summary(.null_base, .res$n_obs, .res$n_unseen), "null_universe_summary.tsv")
  cat("  null_universe_summary.tsv written\n")
} else if (!file_exists(null_counters_file)) {
  cat("  no --null_counters provided, skipping p.emp\n")
}

# ── 2h. Position-level multiple testing: p.adj_bh ─────────────────────────────
# p.adj_bh: BH over the permutation family, one test per (Gene, Position) (the side rows of a position share one pooled
# p.emp and enter BH once). The family is the positional universe U of section 2f-ter: the positions the null scores in
# >= 1 cycle inside the tested set, plus every observed position with a p.emp (the ones the null never scored included,
# at p = 1/(N+1)). p.emp's statistic is "max-side CAAS_score if detected, -Inf otherwise", so a null-detectable position
# the observed data did not detect is a tested position with p = 1; restricting BH to observed-detected positions would
# select on the statistic itself. Positions detected neither by the null nor by the observed data are left out, so m
# counts only columns the finite null sample happened to reach. p.adj_bh is computed with p.emp in null_position_p.R.

# ── 2i. FADE (gene-level - see section 4d) ────────────────────────────────────
# FADE is read at gene level here (max Bayes Factor per gene across sites); significance is BF >= fade_bf_threshold.
# fade_top and fade_bottom hold the FADE results of the top and the bottom direction.
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

# ── 2j. FADE site-level (position-level BF) ───────────────────────────────────
# Per-site max BF from FADE_JSON_TO_CSV (parse_fade_json_sites.R): comma-delimited, with columns gene,
# position (0-based, as in position_scores.tsv), max_bf and target_aa. A 1-based `site` column, when
# the file has one instead of `position`, is shifted by one.
has_fade_site_top <- file_exists(fade_site_top_file)
has_fade_site_bot <- file_exists(fade_site_bot_file)

.load_fade_sites <- function(path) {
  # The file is comma-delimited (write.csv in parse_fade_json_sites.R); read_tsv() would parse each
  # line as a single column.
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
    df$fade_biased <- df$max_site_bf >= fade_bf_threshold
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

# All positions are retained; the directional splits (top, bottom) are applied per analysis after scoring.

# ── 4. Gene-level scoring ─────────────────────────────────────────────────────

cat("\n─── Gene-level scoring ────────────────────────────────────────\n")

# ── 4a. Gene CAAS Scores: size-adjusted max of CAAS_score per gene ────────────
# Three scores computed from different position subsets:
#   gene_caas_score - all positions (full pool)
#   gene_caas_score_top - positions with side == "top"
#   gene_caas_score_bottom - positions with side == "bottom"
#
# The scores (raw maximum, a "both" position counted
# once in the undirected score) come from core.scores via observed_core_scores.py; NA when
# the gene has no scored position in the direction. n_positions* are descriptors.
core_gene <- read_tsv(core_genes_file, show_col_types = FALSE,
                      col_types = cols(Gene = col_character(), .default = col_double()))
.core_gene <- function(col, g) core_gene[[col]][match(g, core_gene$Gene)]
.pos_undirected <- pos_scores %>%
  filter(!is.na(CAAS_score)) %>%
  group_by(Gene, Position) %>%
  summarise(CAAS_score = max(CAAS_score), .groups = "drop")
cat(sprintf("  scored position counts: all=%d, top=%d, bottom=%d positions\n",
            nrow(.pos_undirected),
            sum(pos_scores$side == "top"    & !is.na(pos_scores$CAAS_score)),
            sum(pos_scores$side == "bottom" & !is.na(pos_scores$CAAS_score))))

.gene_undirected <- .pos_undirected %>%
  group_by(Gene) %>%
  summarise(
    gene_caas_score = .core_gene("gene_caas_score", dplyr::first(Gene)),
    n_positions     = dplyr::n_distinct(Position),
    .groups = "drop"
  )
gene_caas <- pos_scores %>%
  group_by(Gene) %>%
  summarise(
    gene_caas_score_top    = .core_gene("gene_caas_score_top",    dplyr::first(Gene)),
    gene_caas_score_bottom = .core_gene("gene_caas_score_bottom", dplyr::first(Gene)),
    n_positions_top    = sum(side == "top",    na.rm = TRUE),
    n_positions_bottom = sum(side == "bottom", na.rm = TRUE),
    max_hypotheses     = if ("n_hypotheses" %in% names(pos_scores)) max(n_hypotheses, na.rm = TRUE) else NA_integer_,
    mean_hypotheses    = if ("n_hypotheses" %in% names(pos_scores)) round(mean(n_hypotheses, na.rm = TRUE), 1) else NA_real_,
    .groups = "drop"
  ) %>%
  left_join(.gene_undirected, by = "Gene")

cat(sprintf("  gene_caas_score: %d genes (%d with top positions, %d with bottom)\n",
            nrow(gene_caas),
            sum(!is.na(gene_caas$gene_caas_score_top)),
            sum(!is.na(gene_caas$gene_caas_score_bottom))))

# ── 4b. Gene Accumulation Score (optional) ────────────────────────────────────
# Reads accumulation_<direction>_results.csv for direction in {all, top, bottom};
# ct_accumulation.nf runs the three directions and stages all of them in accum_dir. "all" pools every
# position with side != "none"; "top" / "bottom" restrict to that side, so the flag is
# direction-aware like FADE and RER.
#
# Returns a tibble with Gene + accum_p<suffix> (the empirical p),
# accum_fdr<suffix> / accum_significant<suffix> (suffix = "" for "all",
# "_top" / "_bottom" otherwise); the suffixes let the three directions full_join onto gene_scores without colliding.
compute_accum_significance <- function(accum_dir, direction, suffix) {
  empty <- tibble(
    Gene = character(),
    !!paste0("accum_p", suffix)       := numeric(),
    !!paste0("accum_fdr", suffix)     := numeric(),
    !!paste0("accum_significant", suffix) := logical()
  )

  pattern <- paste0("accumulation_", direction, "_results.csv")
  f <- list.files(accum_dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) {
    cat(sprintf("Accumulation (%s): file not found (%s), skipping\n", direction, pattern))
    return(list(df = empty, ok = FALSE))
  }

  cat(sprintf("Loading accumulation (%s) from: %s\n", direction, basename(f[1])))
  d <- read_csv(f[1], show_col_types = FALSE)

  gcol <- intersect(c("Gene", "gene"), names(d))[1]
  pval_col <- grep("^PValue$", names(d), value = TRUE)[1]
  if (is.na(gcol) || is.na(pval_col)) {
    cat(sprintf("  Accumulation (%s): missing Gene or PValue column, skipping\n", direction))
    return(list(df = empty, ok = FALSE))
  }

  pvals <- as.numeric(d[[pval_col]])
  out <- tibble(
    Gene = d[[gcol]],
    !!paste0("accum_p", suffix)       := pvals
  )

  fp_col <- paste0("accum_p", suffix)
  # BH FDR on genes with at least one observed CAAS (ActualCount > 0). Genes
  # with no CAAS have p = 1 by construction - they are background members
  # only and must not enter the FDR denominator.
  act_col <- grep("ActualCount", names(d), value = TRUE)[1]
  if (is.na(act_col)) {
    cat(sprintf("  Accumulation (%s): missing ActualCount column, skipping\n", direction))
    return(list(df = empty, ok = FALSE))
  }
  tested <- !is.na(out[[fp_col]]) & !is.na(d[[act_col]]) & as.numeric(d[[act_col]]) > 0
  fdr_q  <- rep(NA_real_, nrow(out))
  if (any(tested)) fdr_q[tested] <- p.adjust(out[[fp_col]][tested], method = "BH")
  out[[paste0("accum_fdr", suffix)]] <- fdr_q
  out[[paste0("accum_significant", suffix)]] <- !is.na(fdr_q) & fdr_q < accumulation_fdr_threshold

  cat(sprintf("  Accumulation (%s): %d genes, %d significant (exact hypergeometric p, BH FDR < %g)\n",
              direction, nrow(out), sum(out[[paste0("accum_significant", suffix)]], na.rm = TRUE), accumulation_fdr_threshold))
  list(df = out, ok = TRUE)
}

has_accum_dir <- file_exists(accum_dir) && dir.exists(accum_dir)
if (!has_accum_dir) cat("Accumulation: directory not available, skipping\n")

.accum_all    <- if (has_accum_dir) compute_accum_significance(accum_dir, "all",    "")        else list(df = tibble(Gene = character(), accum_p = numeric(), accum_fdr = numeric(), accum_significant = logical()), ok = FALSE)
.accum_top    <- if (has_accum_dir) compute_accum_significance(accum_dir, "top",    "_top")    else list(df = tibble(Gene = character(), accum_p_top = numeric(), accum_fdr_top = numeric(), accum_significant_top = logical()), ok = FALSE)
.accum_bottom <- if (has_accum_dir) compute_accum_significance(accum_dir, "bottom", "_bottom") else list(df = tibble(Gene = character(), accum_p_bottom = numeric(), accum_fdr_bottom = numeric(), accum_significant_bottom = logical()), ok = FALSE)

has_accum <- .accum_all$ok
gene_rand <- .accum_all$df %>%
  full_join(.accum_top$df,    by = "Gene") %>%
  full_join(.accum_bottom$df, by = "Gene")

# ── 4c. Gene RERConverge Score (optional) ─────────────────────────────────────
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
  if (!"p.perm" %in% names(rer)) rer$p.perm <- NA_real_
  pval_col <- "p.perm"
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

  # rer_significant is the nominal call (p.perm <= 0.05, not corrected for multiple testing), which RER gene lists
  # and the AMI flags use. rer_perm_padj is its BH adjustment over all genes (Saputra et al. 2021 correct the
  # permulation p-values before calling genes) and rer_significant_fdr the call on it. With few permulations the
  # floor of p.perm keeps rer_perm_padj high: it is about 1 / (N/2 + 1) per gene, so BH cannot go below
  # floor * genes / (genes at the floor).
  if (pval_col == "p.perm") {
    rer$rer_perm_padj <- if ("p.perm.adj" %in% names(rer)) as.numeric(rer[["p.perm.adj"]]) else
      p.adjust(as.numeric(rer[["p.perm"]]), method = "BH")
  } else {
    rer$rer_perm_padj <- NA_real_
  }

  # rer_P is the parametric P of the correlation (RERconverge's getStat takes sign(Rho) * -log10(P) as the gene statistic).
  rer$rer_P <- if ("p" %in% names(rer)) as.numeric(rer[["p"]]) else NA_real_

  gene_rer <- rer %>%
    mutate(
      rer_min_pval     = as.numeric(.data[[pval_col]]),
      rer_significant  = rer_min_pval < nominal_threshold,
      rer_significant_fdr = !is.na(rer_perm_padj) & rer_perm_padj < fdr_threshold,
      rer_rho          = as.numeric(.data[[rho_col]]),
      rer_acceleration = case_when(
        is.na(rer_rho) ~ NA_character_,
        rer_rho > 0    ~ "accelerated",
        rer_rho < 0    ~ "decelerated",
        TRUE           ~ "neutral"
      )
    ) %>%
    select(Gene = gene, rer_min_pval, rer_P, rer_significant, rer_perm_padj, rer_significant_fdr, rer_rho, rer_acceleration) %>%
    mutate(Gene = as.character(Gene))

  cat(sprintf("  RERConverge: %d genes, %d significant (p <= 0.05, uncorrected), %d with BH-adjusted p <= 0.05\n",
              nrow(gene_rer), sum(gene_rer$rer_significant, na.rm = TRUE),
              sum(gene_rer$rer_significant_fdr, na.rm = TRUE)))
} else {
  cat("RERConverge: not available, skipping\n")
  gene_rer <- tibble(Gene = character(), rer_min_pval = numeric(), rer_P = numeric(),
                     rer_significant = logical(), rer_perm_padj = numeric(),
                     rer_significant_fdr = logical(), rer_rho = numeric(),
                     rer_acceleration = character())
}

# ── 4d. Gene FADE Significance (optional) ─────────────────────────────────────
# BF >= fade_bf_threshold per direction. fade_significant_top / _bottom used separately
# to characterise top-CAAS and bottom-CAAS gene sets respectively.
gene_fade <- tibble(Gene = character())

if (nrow(fade_top_df) > 0) {
  gene_fade <- gene_fade %>%
    full_join(fade_top_df %>% mutate(fade_significant_top = fade_max_bf_top >= fade_bf_threshold),
              by = "Gene")
  cat(sprintf("  FADE top: %d significant (BF >= fade_bf_threshold)\n",
              sum(gene_fade$fade_significant_top, na.rm = TRUE)))
} else {
  gene_fade$fade_max_bf_top      <- NA_real_
  gene_fade$fade_significant_top <- NA
}

if (nrow(fade_bottom_df) > 0) {
  gene_fade <- gene_fade %>%
    full_join(fade_bottom_df %>% mutate(fade_significant_bottom = fade_max_bf_bottom >= fade_bf_threshold),
              by = "Gene")
  cat(sprintf("  FADE bottom: %d significant (BF >= fade_bf_threshold)\n",
              sum(gene_fade$fade_significant_bottom, na.rm = TRUE)))
} else {
  gene_fade$fade_max_bf_bottom      <- NA_real_
  gene_fade$fade_significant_bottom <- NA
}

# ── 4e. Assemble gene scores ──────────────────────────────────────────────────
# full_join, not left_join: gene_caas only covers genes with a CAAS-detected position, whereas
# gene_rand, gene_rer and gene_fade each have their own, larger gene universe (RERConverge scores every
# gene with a sufficiently populated gene tree). A left_join rooted at gene_caas would drop every gene
# with RER, FADE or accumulation evidence but no CAAS position, and the per-module fcs_stats_*.tsv files
# would then be truncated to the CAAS gene set instead of the universe of each module.
gene_scores <- gene_caas

if (nrow(gene_rand) > 0) {
  gene_scores <- gene_scores %>% full_join(gene_rand, by = "Gene")
} else {
  gene_scores$accum_p                  <- NA_real_
  gene_scores$accum_fdr                <- NA_real_
  gene_scores$accum_significant        <- NA
  gene_scores$accum_p_top              <- NA_real_
  gene_scores$accum_fdr_top            <- NA_real_
  gene_scores$accum_significant_top    <- NA
  gene_scores$accum_p_bottom           <- NA_real_
  gene_scores$accum_fdr_bottom         <- NA_real_
  gene_scores$accum_significant_bottom <- NA
}

if (nrow(gene_rer) > 0) {
  gene_scores <- gene_scores %>% full_join(gene_rer, by = "Gene")
} else {
  gene_scores$rer_min_pval      <- NA_real_
  gene_scores$rer_P             <- NA_real_
  gene_scores$rer_significant   <- NA
  gene_scores$rer_perm_padj     <- NA_real_
  gene_scores$rer_significant_fdr <- NA
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

# ── 5. Correlation analysis (gene-level) ──────────────────────────────────────

cat("\n─── Correlation analysis ──────────────────────────────────────\n")

# Correlations use gene_caas_score as the primary axis. Accumulation, RER and FADE
# are represented by their native significance (accum_p / RER p / FADE BF),
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

# ── 6. Write outputs ──────────────────────────────────────────────────────────

cat("\n─── Writing outputs ───────────────────────────────────────────\n")

# Position scores. The per-row factors (asr_score, derived_agreement) and asr_path_score are not
# written: CAAS_score is the position-level number (the asr_path_score of a position equals it), and
# the breakdown lives in filtered_discovery.tsv. top_species_residues, bottom_species_residues
# and n_conserved_pairs (alignment-based, over the contrast pairs of the hypotheses that call the position)
# are the raw-residue descriptors.
pos_out <- pos_scores %>%
  select(Gene, Position,
         biochem_groupings, is_biochem,
         any_of(c("n_hypotheses", "participating_hypotheses",
                  "top_species_residues", "bottom_species_residues",
                  "n_top_species", "n_bottom_species", "n_conserved_pairs")), CAAS_score,
         side,
         any_of("caas"), any_of(c("ancestral_aa", "derived_aa")),
         # p.emp: the pooled "detects AND exceeds" position p; p.adj_bh: its BH adjustment (§2h).
         any_of(c("p.emp", "p.adj_bh")),
         # nd: cycles of the null that score the position; null_unseen: nd == 0. Null-adequacy diagnostics, no demotion.
         any_of(c("nd", "null_unseen"))) %>%
  arrange(desc(CAAS_score))

write_tsv(pos_out, "position_scores.tsv")
cat(sprintf("  position_scores.tsv: %d rows\n", nrow(pos_out)))

# Gene scores
gene_out <- gene_scores %>%
  select(
    Gene,
    n_positions, n_positions_top, n_positions_bottom,
    gene_caas_score, gene_caas_score_top, gene_caas_score_bottom,
    any_of(c("accum_p", "accum_fdr", "accum_significant",
             "accum_p_top", "accum_fdr_top", "accum_significant_top",
             "accum_p_bottom", "accum_fdr_bottom", "accum_significant_bottom")),
    any_of(c("rer_min_pval", "rer_P", "rer_significant", "rer_perm_padj", "rer_significant_fdr", "rer_rho", "rer_acceleration")),
    any_of(c("fade_max_bf_top", "fade_significant_top",
             "fade_max_bf_bottom", "fade_significant_bottom"))
  ) %>%
  arrange(desc(gene_caas_score))

write_tsv(gene_out, "gene_scores.tsv")
cat(sprintf("  gene_scores.tsv: %d rows\n", nrow(gene_out)))

# ── FCS stats table ───────────────────────────────────────────────────────────
# Read by 12.FCS_general_report.Rmd. Columns: gene; score_<ranking>, zero-floored downstream over the
# cleaned_background universe; flag_<name>, per-gene booleans. score_global, score_top and score_bottom
# are the *_all gene scores (magnitude); the top and bottom ones use only the positions of that side.
# The flags annotate the leading edge and never gate the FCS input.
.has  <- function(d, c) c %in% names(d)
.col  <- function(d, c, default = NA) if (.has(d, c)) d[[c]] else rep(default, nrow(d))
.istrue <- function(x) !is.na(x) & x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
.rer_dir <- tolower(as.character(.col(gene_scores, "rer_acceleration")))

# The CAAS rankings are the gene scores: gene_caas_score, gene_caas_score_top and gene_caas_score_bottom.
fcs_stats <- tibble(
  gene             = gene_scores$Gene,
  score_global     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score"))),
  score_top        = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_top"))),
  score_bottom     = suppressWarnings(as.numeric(.col(gene_scores, "gene_caas_score_bottom"))),
  flag_fade_top    = .istrue(.col(gene_scores, "fade_significant_top")),
  flag_fade_bottom = .istrue(.col(gene_scores, "fade_significant_bottom")),
  flag_rer_acc     = .istrue(.col(gene_scores, "rer_significant")) & grepl("acc", .rer_dir),
  flag_rer_decc    = .istrue(.col(gene_scores, "rer_significant")) & grepl("dec", .rer_dir),
  flag_rer_acc_adjusted = .istrue(.col(gene_scores, "rer_significant_fdr")) & grepl("acc", .rer_dir),
  flag_rer_decc_adjusted = .istrue(.col(gene_scores, "rer_significant_fdr")) & grepl("dec", .rer_dir),
  flag_accum        = .istrue(.col(gene_scores, "accum_significant")),
  flag_accum_top    = .istrue(.col(gene_scores, "accum_significant_top")),
  flag_accum_bottom = .istrue(.col(gene_scores, "accum_significant_bottom"))
) %>%
  mutate(flag_fade = flag_fade_top | flag_fade_bottom)
write_tsv(fcs_stats, "fcs_stats.tsv")
cat(sprintf("  fcs_stats.tsv: %d genes (%d top, %d bottom)\n",
            nrow(fcs_stats), sum(!is.na(fcs_stats$score_top)),
            sum(!is.na(fcs_stats$score_bottom))))

# ── Per-module FCS ranking files ──────────────────────────────────────────────
# Each carries only the score_<ranking> columns of its module. The cross-module annotation (CAAS scores,
# FADE, accumulation flags) reaches the FCS report separately, as annot_file = fcs_stats.tsv. SCORING can
# thus render one FCS report per module, ranked by that module's own statistic and annotated with the
# evidence of every module.
.nonempty <- function(col) .has(gene_scores, col) && any(!is.na(gene_scores[[col]]))

if (has_rer) {
  # The RER ranking is the signed parametric statistic sign(Rho) * -log10(P), the gene statistic of RERconverge (getStat) and the
  # one its permulations (corStat) carry for the null.
  .rp <- suppressWarnings(as.numeric(.col(gene_scores, "rer_P")))
  .rr <- suppressWarnings(as.numeric(.col(gene_scores, "rer_rho")))
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_global = sign(.rr) * -log10(.rp)
  ), "fcs_stats_rer.tsv")
  cat("  fcs_stats_rer.tsv written (RER global ranking)\n")
}

if (.nonempty("fade_max_bf_top") || .nonempty("fade_max_bf_bottom")) {
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_top    = suppressWarnings(as.numeric(.col(gene_scores, "fade_max_bf_top"))),
    score_bottom = suppressWarnings(as.numeric(.col(gene_scores, "fade_max_bf_bottom")))
  ), "fcs_stats_fade.tsv")
  cat("  fcs_stats_fade.tsv written (FADE BF rankings)\n")
}

if (.nonempty("accum_p")) {
  # FCS ranks in descending order (higher = more significant) while empirical p runs the other way, so
  # the ranking axis is -log10(accum_p), derived here and not stored as a column of gene_scores.
  .afp <- suppressWarnings(as.numeric(.col(gene_scores, "accum_p")))
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_global = ifelse(is.na(.afp), NA_real_, -log10(pmax(.afp, 1e-300)))
  ), "fcs_stats_accum.tsv")
  cat("  fcs_stats_accum.tsv written (accumulation rankings, -log10 p)\n")
}

# Correlations
write_tsv(corr_results, "gene_correlations.tsv")


cat("Base scoring complete; calibrated evidence and activity lists are materialized next.\n")
