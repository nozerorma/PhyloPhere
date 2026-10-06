#!/usr/bin/env Rscript
# scoring_compute.R — Position-level and gene-level CAAS scores, integrated with FADE, RER and accumulation.
# PhyloPhere | subworkflows/SCORING/local/src/
# =============================================================================
# Called by:  SCORING_COMPUTE Nextflow process (scoring_compute.nf → Rscript scoring_compute.R ...),
#             after observed_core_scores.py, whose two tables it reads
#
# The position score (CAAS_score) and the gene score (size_adj_max) come from core.scores via
# observed_core_scores.py. This script adds the permulation p-values of positions (p.emp and its
# BH and permutation-FDR adjustments) and genes (gene_caas_pperm), joins the per-gene evidence of
# FADE, RERConverge and accumulation, and writes the tables, the ranked slices and the enrichment
# curves that the ENRICHMENT subworkflow and the reports read. It runs once on the full pool of
# filtered_discovery.tsv; direction is carried by the `side` column (top, bottom, none).
#
# Args (named flags, from task.script). An input whose file name starts with NO_ counts as absent.
#   --postproc            filtered_discovery.tsv (mandatory)
#   --core_positions      observed_core_scores.py positions table (mandatory)
#   --core_genes          observed_core_scores.py genes table (mandatory)
#   --fade_top, --fade_bottom         gene-level FADE summaries (fade_summary_{top,bottom}.tsv)
#   --fade_site_top, --fade_site_bot  per-site FADE tables (comma-delimited)
#   --rer                 rerconverge_summary_<trait>.tsv
#   --accum_dir           directory with accumulation_<direction>_<scheme>_aggregated_results.csv
#   --caas_perms          caas_perms.rds, the gene-level null (scoring_caas_perms.R)
#   --caas_pos_cycle_caas perm_pos_cycle_caas.tsv.gz, the position-level null for p.emp
#   --gene_perm_pooled    true|false: also write the n-stratified pooled gene p (default false)
#   --p_emp_thr           threshold of flag_caas_significant (default 0.05)
#   --hypotheses_pairs, --top_pct, --gene_top_pct   parsed but not used by the computation
#
# Outputs (working directory):
#   position_scores.tsv            one row per Gene, Position and side
#   gene_scores.tsv                one row per Gene
#   gene_correlations.tsv          pairwise correlations between gene scores (header only while one score exists)
#   fcs_stats.tsv                  gene scores and flags for the FCS reports
#   fcs_stats_{rer,fade,accum}.tsv per-module FCS rankings, written when that evidence exists
#   gene_lists/slice_*.tsv         12 ranked gene slices (top, bottom, global x 25, 10, 5, 1%)
#   position_lists/slice_*.tsv     12 ranked position slices, same layout
#   gene_threshold_enrichment.tsv  odds ratio and Fisher p of FADE, RER and accumulation across CAAS thresholds
#   pos_threshold_enrichment.tsv   the same curve for position-level FADE
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

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
caas_perms_file      <- parse_arg("--caas_perms")  # caas_perms.rds (CAAS permulation-excess null); NO_FILE otherwise
core_positions_file  <- parse_arg("--core_positions")  # observed_core_scores.py: CAAS_score per (Gene, Position, side)
core_genes_file      <- parse_arg("--core_genes")      # observed_core_scores.py: size-adjusted gene CAAS scores
caas_pos_cycle_caas_file <- parse_arg("--caas_pos_cycle_caas")  # perm_pos_cycle_caas.tsv.gz (p.emp numerator/denominator); NO_FILE otherwise
gene_perm_pooled_raw <- parse_arg("--gene_perm_pooled", "false")
# Threshold of flag_caas_significant (fcs_stats.tsv): gene_caas_pperm_adj <= p_emp_thr. It is the value of
# the scoring_p_emp_thr param of 11.Scoring_report.Rmd (conf/scoring.config), so the "% significant"
# figures of the enrichment reports agree with the Scoring report.
p_emp_thr         <- as.numeric(parse_arg("--p_emp_thr", "0.05"))
# Rows reach this script already pooled over hypotheses by CT_DISAMBIGUATION (one row per Gene,
# Position, scheme and side, hypothesis NA), so scoring never pools hypotheses itself.
top_pct           <- as.numeric(parse_arg("--top_pct",  "0.10"))
top25_pct         <- 0.25
top5_pct          <- 0.05
top1_pct          <- 0.01
gene_top_pct      <- as.numeric(parse_arg("--gene_top_pct",  "0.10"))
gene_top25_pct    <- 0.25
gene_top5_pct     <- 0.05
gene_top1_pct     <- 0.01
gene_perm_pooled      <- tolower(as.character(gene_perm_pooled_raw)) %in% c("true", "1", "yes")
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
# (caap_group, asr_path_score, side, ...), used as they are; position_scores.tsv keeps the same schema.
cat(sprintf("  %d rows, %d unique Gene×Position pairs\n",
            nrow(df), n_distinct(paste(df$Gene, df$Position))))

# ── 2. Position-level scoring ─────────────────────────────────────────────────

# ── 2a. Scheme scope ──────────────────────────────────────────────────────────
# The five scoring schemes.
#
# No per-scheme weight: how many of the five schemes detect a substitution is a
# deterministic property of which amino acids are involved (a discretised
# biochemical distance, see the report's Biochemistry tab), not evidence
# strength. Section 2g aggregates schemes with a mean for this reason.
scoring_schemes <- c("US", "GS4", "GS3", "GS2", "GS1")

# Priority only picks the representative scheme whose display columns (side, caap_group, ...) the
# Gene×Position aggregation of section 2g carries. It is separate from the scoring itself, which treats
# all five schemes symmetrically.
scheme_priority_int <- c(US = 5, GS4 = 4, GS3 = 3, GS2 = 2, GS1 = 1)

df <- df %>%
  mutate(
    scheme_priority = scheme_priority_int[caap_group]
  ) %>%
  filter(caap_group %in% scoring_schemes)
cat(sprintf("  %d rows across %d scoring schemes after dropping non-scoring schemes\n",
            nrow(df), n_distinct(df$caap_group)))

# ── 2b. Hypothesis pooling (already done upstream) ────────────────────────────
# Rows arrive pooled over hypotheses by CT_DISAMBIGUATION (fop_pool.pool_domains). H1..Hn are
# overlapping K-pair designs over the same Voronoi domains, not independent replicates: averaging over
# them in section 2g would dilute a strong canonical signal and let a position with many harvested
# hypotheses distort every genome-wide rank. The pooling averages the per-hypothesis domain scores over
# the K fixed domains, weighted by the PSS of contrast_hypotheses_pairs.tsv; a single contrast reduces
# to the plain PSS-weighted domain mean.
cat("  FOP pooling: done in-tree per (Gene, Position, scheme, side) [core v3]\n")
# Backfill the descriptor columns that section 2g and the reports expect, so that a missing column is
# never hit.
if (!"n_hypotheses" %in% names(df))          df$n_hypotheses <- 1L
if (!"participating_hypotheses" %in% names(df)) df$participating_hypotheses <- ""
# Descriptors of the species that carry the change; empty text when the input has none.
for (.c in c("top_species_residues", "bottom_species_residues",
             "n_top_species", "n_bottom_species")) {
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
# Sort descending by scheme_priority (US > GS4 > GS3 > GS2 > GS1) so first()
# deterministically picks the US scheme (falling back to GS4..GS1) for
# display/gating-only columns (asr_is_conserved, etc.).
# Priority is display-only and never enters a scored quantity.
df <- df %>% arrange(desc(scheme_priority))

# `side` belongs to the aggregation key: a position detected on both sides has two rows, and
# CAAS_score (the mean over the five schemes) is taken per side.
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
    n_schemes          = dplyr::n(),
    scheme_set         = paste(sort(unique(as.character(caap_group))), collapse = "+"),
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
    # The per-scheme factors (asr_score / caas_row) and the ASR diagnostic columns (asr_path_score,
    # derived_agreement) are not carried to the position level: CAAS_score is the mean of asr_path_score
    # over the schemes, and a position-level mean of each sub-factor would hide scheme disagreement (a
    # split V->{I,L} shows derived_agreement ~ 0.9 when US strongly disagrees). They stay per
    # (Gene, Position, caap_group) in `df` for anything that needs the breakdown.
    caap_group         = first(caap_group),
    .groups = "drop"
  )

# CAAS_score = mean of asr_path_score over the schemes that scored the (Gene, Position, side),
# computed by core.scores.
core_pos <- read_tsv(core_positions_file, show_col_types = FALSE,
                     col_types = cols(Gene = col_character(), Position = col_integer(),
                                      side = col_character(), CAAS_score = col_double()))
.pk  <- function(g, p, sd) paste(g, p, sd, sep = "\r")
.hit <- match(.pk(pos_scores$Gene, pos_scores$Position, pos_scores$side),
              .pk(core_pos$Gene, core_pos$Position, core_pos$side))
if (anyNA(.hit) || nrow(core_pos) != nrow(pos_scores)) {
  stop(sprintf("core positions (%d) and scored positions (%d) disagree: observed_core_scores.py and this script read different rows",
               nrow(core_pos), nrow(pos_scores)))
}
pos_scores$CAAS_score <- core_pos$CAAS_score[.hit]
rm(core_pos, .hit)

# The asr_path_score of a position row is the mean over its schemes for that side, which is CAAS_score
# itself (caas_row is the row's asr_path_score).
pos_scores <- pos_scores %>% mutate(asr_path_score = CAAS_score)

# Ancestral and derived residues of each (Gene, Position, side), as the ASR inferred them. They are read from
# the K fixed domains of the US scheme row: domain_<d>_anc_aa is the ancestral residue of domain d and
# domain_<d>_top_aa / domain_<d>_bot_aa the derived one on the side that carries the change. Residues are pooled
# over the domains, most frequent first, comma-separated. Empty text when the position has no US row, no
# resolved domain, or side "none" for the derived set. The `caas` pattern (top/bottom residues) has no direction.
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
.rep_us <- df %>%
  filter(caap_group == "US", !duplicated(paste(Gene, Position, side, sep = "\r"))) %>%
  select(Gene, Position, side, any_of(c(.anc_cols, sub("_anc_aa$", "_top_aa", .anc_cols), sub("_anc_aa$", "_bot_aa", .anc_cols))))
if (length(.anc_cols) && nrow(.rep_us)) {
  .mat <- function(cols) as.matrix(.rep_us[, cols, drop = FALSE])
  .der <- .mat(sub("_anc_aa$", "_top_aa", .anc_cols))
  .bot <- .mat(sub("_anc_aa$", "_bot_aa", .anc_cols))
  .der[.rep_us$side == "bottom", ] <- .bot[.rep_us$side == "bottom", ]
  .der[!.rep_us$side %in% c("top", "bottom"), ] <- NA_character_
  .rep_us <- .rep_us %>% transmute(Gene, Position, side, ancestral_aa = .pool_aa(.mat(.anc_cols)), derived_aa = .pool_aa(.der))
  pos_scores <- pos_scores %>% left_join(.rep_us, by = c("Gene", "Position", "side"))
  pos_scores$ancestral_aa[is.na(pos_scores$ancestral_aa)] <- ""
  pos_scores$derived_aa[is.na(pos_scores$derived_aa)] <- ""
} else {
  pos_scores$ancestral_aa <- ""
  pos_scores$derived_aa <- ""
}
rm(.rep_us, .anc_cols)

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
# The pooled statistic is the max-over-sides "all" axis -- identical to
# .pos_undirected (§4a) on the observed side and _build_cycle_score_pools'
# pc["all"] on the null. The null's per-cycle score is the caas_score column, the same core.scores value the
# observed position gets. Values within TIE_TOL of the observed score count as ties (>=).
pos_scores$p.emp <- NA_real_
has_caas_pos_cycle_caas <- file_exists(caas_pos_cycle_caas_file)
# A null table with a header and no row (N = 0, or no permuted cycle re-detected any position) has no cycle to count:
# (k + 1) / (N + 1) would read 1 for every position, a value rather than the absence of one. p.emp, p.adj_bh and
# p.adj_sam stay NA, as when no null is given.
if (has_caas_pos_cycle_caas && length(read_lines(caas_pos_cycle_caas_file, n_max = 2)) < 2) {
  cat("  perm_pos_cycle_caas.tsv.gz has no row: no null cycle, so p.emp, p.adj_bh and p.adj_sam stay NA\n")
  has_caas_pos_cycle_caas <- FALSE
}
# caas_perms.rds is loaded here when present (its columns are the cycle roster for N below) and
# reused by the gene-level p in §4f.
caas_perms <- NULL
if (has_caas_pos_cycle_caas) {
  cat("Loading per-cycle CAAS null (p.emp):", caas_pos_cycle_caas_file, "\n")
  # caas_score is read as text and converted with as.numeric: readr's own double parser
  # is off by an ulp for ~13% of 17-digit values.
  cyc_caas <- read_tsv(caas_pos_cycle_caas_file, show_col_types = FALSE,
                       col_types = cols(.default = col_guess(), caas_score = col_character()))
  if (!"caas_score" %in% names(cyc_caas)) {
    stop("perm_pos_cycle_caas.tsv.gz has no caas_score column: it predates the shared position score. ",
         "Regenerate the CAAS permulation null.")
  }
  cyc_caas <- cyc_caas %>%
    mutate(Position = as.integer(Position),
           caas_score = suppressWarnings(as.numeric(caas_score))) %>%
    filter(!is.na(caas_score))

  # per-cycle pooled null statistic = max over detected sides of the scheme-mean
  cyc_pooled <- cyc_caas %>%
    group_by(Gene, Position, cycle) %>%
    summarise(caas_max = max(caas_score), .groups = "drop")
  rm(cyc_caas)

  # observed pooled statistic = best side per (Gene, Position) (== .pos_undirected)
  obs_max <- pos_scores %>%
    filter(!is.na(CAAS_score)) %>%
    group_by(Gene, Position) %>%
    summarise(.obs = max(CAAS_score), .groups = "drop")

  # N = cycles replayed, not cycles present in the file: a cycle that re-detects
  # no position emits no rows. The replay roster is the column set of
  # caas_perms.rds (one column per cycle, zero-detection cycles included). It is
  # used only if it contains every cycle present in the file; otherwise N falls
  # back to the cycles present, which can only make p.emp larger.
  N_present <- dplyr::n_distinct(cyc_pooled$cycle)
  N_emp     <- N_present
  if (file_exists(caas_perms_file)) {
    caas_perms <- tryCatch(readRDS(caas_perms_file), error = function(e) NULL)
    .roster <- colnames(caas_perms[["caas_corStat_byrank"]][["global"]])
    if (!is.null(.roster) && all(unique(cyc_pooled$cycle) %in% .roster)) {
      N_emp <- length(.roster)
    } else {
      cat("  WARNING: caas_perms.rds cycle roster missing or inconsistent with ",
          "perm_pos_cycle_caas.tsv.gz; N = cycles present in the null file\n",
          file = stderr())
    }
  }

  # An observed position that no null cycle re-detects has k_emp = 0 (its
  # null statistic is -Inf in every cycle), hence the left join from obs_max.
  .k_emp <- obs_max %>%
    left_join(cyc_pooled, by = c("Gene", "Position")) %>%
    group_by(Gene, Position) %>%
    summarise(k_emp    = sum(caas_max >= .obs - TIE_TOL, na.rm = TRUE),
              null_hit = any(!is.na(caas_max)),
              .groups  = "drop") %>%
    mutate(p.emp = (k_emp + 1) / (N_emp + 1))

  .n_obs_pos_e <- nrow(.k_emp)
  .n_matched_e <- sum(.k_emp$null_hit)
  .rate_e <- if (.n_obs_pos_e > 0) .n_matched_e / .n_obs_pos_e else 0
  cat(sprintf("  p.emp: %d/%d observed positions re-detected by the null (%.1f%%), N=%d cycles (%d with detections)\n",
              .n_matched_e, .n_obs_pos_e, 100 * .rate_e, N_emp, N_present))
  # The overlap rate only diagnoses a coordinate mismatch when there are enough
  # observed positions to estimate it. Below P_EMP_GUARD_MIN_POS an unmatched
  # position is taken at face value: no null cycle re-detects it, so k_emp = 0.
  P_EMP_GUARD_MIN_POS <- 10L
  if (.n_obs_pos_e < P_EMP_GUARD_MIN_POS && .n_matched_e < .n_obs_pos_e) {
    cat(sprintf(paste0("  p.emp: %d of %d observed position(s) never re-detected by the ",
                       "null; scored at k_emp = 0 (too few positions for the ",
                       "coordinate-mismatch check)\n"),
                .n_obs_pos_e - .n_matched_e, .n_obs_pos_e))
  }
  if (.n_obs_pos_e >= P_EMP_GUARD_MIN_POS && .rate_e < 0.5) {
    # Low overlap means the two tables are on different coordinate systems, not
    # that the unmatched positions are strong: leave them untested (NA).
    .k_emp$p.emp[!.k_emp$null_hit] <- NA_real_
    cat(sprintf(paste0("  WARNING: p.emp join rate %.1f%% < 50%% -- the observed ",
                       "positions (filtered_discovery.tsv) and the null's ",
                       "perm_pos_cycle_caas.tsv.gz positions are likely on different ",
                       "coordinate systems. Unmatched positions left NA; treat ",
                       "p.emp/p.adj_bh/p.adj_sam as unreliable.\n"),
                100 * .rate_e), file = stderr())
  }
  pos_scores <- pos_scores %>%
    select(-any_of("p.emp")) %>%
    left_join(.k_emp %>% select(Gene, Position, p.emp), by = c("Gene", "Position"))
} else if (!file_exists(caas_pos_cycle_caas_file)) {
  cat("  no --caas_pos_cycle_caas provided, skipping p.emp\n")
}

# ── 2h. Position-level multiple testing: p.adj_bh and p.adj_sam ───────────────
# p.adj_bh: BH over the permutation family, one test per (Gene, Position) (the
# side rows of a position share one pooled p.emp and enter BH once). The
# family is every position the null detects in >= 1 cycle plus every observed
# position with a p.emp. p.emp's statistic is "max-side CAAS_score if detected,
# -Inf otherwise", so a null-detectable position the observed data did not
# detect is a tested position with p = 1; restricting BH to observed-detected
# positions would select on the statistic itself. Positions detected neither
# by the null nor by the observed data are left out, so m counts only columns
# the finite null sample happened to reach.
#
# p.adj_sam: permutation FDR (SAM-style, Tusher et al. 2001) on the pooled
# max-over-sides score. For a threshold t,
#   FDR(t) = [#(cycle, position) null pairs with score >= t] / N
#            / #{observed positions with score >= t},
# with pi0 = 1; a position's value is the minimum FDR(t) over thresholds t at or
# below its own score. The expected null count is taken over every position
# each cycle detects, so it needs no family definition and counts columns the
# observed data did not detect.
pos_scores$p.adj_bh  <- NA_real_
pos_scores$p.adj_sam <- NA_real_
if (has_caas_pos_cycle_caas) {
  .fam_e <- cyc_pooled %>%
    distinct(Gene, Position) %>%
    full_join(.k_emp %>% filter(!is.na(p.emp)) %>% select(Gene, Position, p.emp),
              by = c("Gene", "Position")) %>%
    mutate(observed = !is.na(p.emp),
           p.emp    = coalesce(p.emp, 1))
  if (nrow(.fam_e) > 0) {
    .fam_e$p.adj_bh <- p.adjust(.fam_e$p.emp, method = "BH")
    pos_scores <- pos_scores %>%
      select(-p.adj_bh) %>%
      left_join(.fam_e %>% filter(observed) %>% select(Gene, Position, p.adj_bh),
                by = c("Gene", "Position"))
  }
  cat(sprintf("  p.adj_bh: BH over %d positions (%d observed, %d null-only at p = 1)\n",
              nrow(.fam_e), sum(.fam_e$observed), sum(!.fam_e$observed)))

  .sam <- obs_max %>%
    semi_join(.k_emp %>% filter(!is.na(p.emp)), by = c("Gene", "Position"))
  if (nrow(.sam) > 0) {
    .null_sorted <- sort(cyc_pooled$caas_max)
    .s   <- .sam$.obs
    .ord <- order(.s, decreasing = TRUE)
    .s_d <- .s[.ord]
    # null pairs with score >= t, per cycle; observed positions with score >= t
    .exp_null <- (length(.null_sorted) -
                    findInterval(.s_d - TIE_TOL, .null_sorted, left.open = TRUE)) / N_emp
    .n_obs    <- vapply(.s_d, function(t) sum(.s >= t - TIE_TOL), numeric(1))
    .fdr      <- pmin(1, .exp_null / .n_obs)
    .q        <- rev(cummin(rev(.fdr)))
    .sam$p.adj_sam <- NA_real_
    .sam$p.adj_sam[.ord] <- .q
    pos_scores <- pos_scores %>%
      select(-p.adj_sam) %>%
      left_join(.sam %>% select(Gene, Position, p.adj_sam), by = c("Gene", "Position"))
  }
  cat(sprintf("  p.adj_sam: permutation FDR over %d observed positions, %d null (cycle, position) pairs, N=%d\n",
              nrow(.sam), nrow(cyc_pooled), N_emp))
  rm(cyc_pooled)
}

# ── 2i. FADE (gene-level - see section 4d) ────────────────────────────────────
# FADE is read at gene level here (max Bayes Factor per gene across sites); significance is BF >= 100.
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

# All positions are retained; the directional splits (top, bottom) are applied per analysis after scoring.

# ── 4. Gene-level scoring ─────────────────────────────────────────────────────

cat("\n─── Gene-level scoring ────────────────────────────────────────\n")

# ── 4a. Gene CAAS Scores: size-adjusted max of CAAS_score per gene ────────────
# Three scores computed from different position subsets:
#   gene_caas_score - all positions (full pool)
#   gene_caas_score_top - positions with side == "top"
#   gene_caas_score_bottom - positions with side == "bottom"
#
# The scores (size_adj_max, direction-matched reference pools, a "both" position counted
# once in the undirected score) come from core.scores via observed_core_scores.py; NA when
# the gene has no scored position in the direction. n_positions* are descriptors.
core_gene <- read_tsv(core_genes_file, show_col_types = FALSE,
                      col_types = cols(Gene = col_character(), .default = col_double()))
.core_gene <- function(col, g) core_gene[[col]][match(g, core_gene$Gene)]
.pos_undirected <- pos_scores %>%
  filter(!is.na(CAAS_score)) %>%
  group_by(Gene, Position) %>%
  summarise(CAAS_score = max(CAAS_score), .groups = "drop")
cat(sprintf("  size-adjust reference pools: all=%d, top=%d, bottom=%d positions\n",
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
    gene_caas_score_top_all    = .core_gene("gene_caas_score_top_all",    dplyr::first(Gene)),
    gene_caas_score_bottom_all = .core_gene("gene_caas_score_bottom_all", dplyr::first(Gene)),
    n_positions_top    = sum(side == "top",    na.rm = TRUE),
    n_positions_bottom = sum(side == "bottom", na.rm = TRUE),
    max_hypotheses     = if ("n_hypotheses" %in% names(pos_scores)) max(n_hypotheses, na.rm = TRUE) else NA_integer_,
    mean_hypotheses    = if ("n_hypotheses" %in% names(pos_scores)) round(mean(n_hypotheses, na.rm = TRUE), 1) else NA_real_,
    .groups = "drop"
  ) %>%
  left_join(.gene_undirected, by = "Gene") %>%
  mutate(
    gene_caas_score_top    = gene_caas_score_top_all,
    gene_caas_score_bottom = gene_caas_score_bottom_all
  )

cat(sprintf("  gene_caas_score: %d genes (%d with top positions, %d with bottom)\n",
            nrow(gene_caas),
            sum(!is.na(gene_caas$gene_caas_score_top)),
            sum(!is.na(gene_caas$gene_caas_score_bottom))))

# ── 4b. Gene Accumulation Score (optional) ────────────────────────────────────
# Reads accumulation_<direction>_<scheme>_aggregated_results.csv for direction in {all, top, bottom};
# ct_accumulation.nf runs the three directions and stages all of them in accum_dir. "all" pools every
# position with side != "none"; "top" / "bottom" restrict to that side, so the flag is
# direction-aware like FADE and RER. Per direction, the per-scheme PValueEmpirical columns are combined
# into one p per gene (Cauchy combination below), BH-adjusted over the genes with at least one CAAS,
# and flagged at FDR < 0.05.
#
# Returns a tibble with Gene + accum_cct_p<suffix> / accum_fdr<suffix> / accum_significant<suffix> and
# accum_pval_<scheme><suffix> (suffix = "" for "all", "_top" / "_bottom" otherwise); the suffixes let
# the three directions full_join onto gene_scores without colliding.
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
  # The five schemes (US, GS4, GS3, GS2, GS1) partition the amino acids but test the same positions: a
  # position counted under one scheme is frequently counted under others, so the per-scheme p-values
  # are positively correlated. CCT keeps its null valid under arbitrary dependence, which a combiner
  # that assumes independence does not. Weights are equal (1 / number of schemes available). The same
  # combiner is applied in accum_gene_lists.nf and 10.Accumulation_report.Rmd.
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

# ── 4d. Gene FADE Significance (optional) ─────────────────────────────────────
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

# ── 4f. Per-gene CAAS permulation p ───────────────────────────────────────────
# gene_caas_score{,_top_all,_bottom_all} against the matching direction's matrix of
# caas_corStat_byrank in caas_perms.rds (genes × N cycles; zeros are kept, so N is always the full cycle
# count). Per gene, p is right-tailed, (k + 1) / (N + 1), and BH-adjusted within direction. It floors
# near 1/(N+1) per gene, so it ranks genes with FDR context and is not a genome-wide-significant call
# by itself; the optional pooled variant below trades that floor for validity only within an
# n-stratum. A gene absent from the null's universe gets NA, never the implicit-zero-row value 1/(N+1).
cat("\n─── CAAS permulation gene p (Tier 1A) ─────────────────────────\n")
for (col in c("gene_caas_pperm", "gene_caas_pperm_top", "gene_caas_pperm_bottom",
              "gene_caas_pperm_adj", "gene_caas_pperm_adj_top", "gene_caas_pperm_adj_bottom")) {
  gene_scores[[col]] <- NA_real_
}
if (gene_perm_pooled) {
  for (col in c("gene_caas_pperm_pooled", "gene_caas_pperm_pooled_top", "gene_caas_pperm_pooled_bottom")) {
    gene_scores[[col]] <- NA_real_
  }
}

if (!file_exists(caas_perms_file)) {
  cat("  no --caas_perms provided, skipping gene p.perm\n")
} else {
  if (is.null(caas_perms)) caas_perms <- tryCatch(readRDS(caas_perms_file), error = function(e) NULL)
  if (is.null(caas_perms) || is.null(caas_perms[["caas_corStat_byrank"]])) {
    cat("  CAAS permulation null unreadable or missing caas_corStat_byrank, skipping gene p.perm\n")
  } else {
    .stat <- caas_perms[["gene_stat"]]
    if (!is.null(.stat) && !identical(as.character(.stat), "size_adj_max")) {
      stop(sprintf(paste0(
        "gene_caas_pperm: CAAS permulation null is STALE (gene_stat='%s'), but the observed ",
        "gene_caas_score uses 'size_adj_max'. Rebuild it (no ASR replay needed, minutes):\n",
        "  python3 subworkflows/CT_DISAMBIGUATION/local/reaggregate_perm_scores.py --detail ",
        "<run>/caas_permulation/perm_pos_detail --output-dir <run>/caas_permulation\n",
        "  Rscript subworkflows/SCORING/local/src/scoring_caas_perms.R --gene-cycle-scores ",
        "<run>/caas_permulation/gene_cycle_scores.tsv --output <run>/caas_permulation/caas_perms.rds"),
        as.character(.stat)))
    }
    byrank <- caas_perms[["caas_corStat_byrank"]]

    # Right-tailed nominal p against a single gene's own null row.
    .perm_p_row <- function(obs, null_row) {
      if (is.na(obs) || length(null_row) == 0) return(NA_real_)
      nr <- null_row[!is.na(null_row)]
      if (length(nr) == 0) return(NA_real_)
      (sum(nr >= obs) + 1) / (length(nr) + 1)
    }
    .apply_direction <- function(obs_col, mat) {
      out <- rep(NA_real_, nrow(gene_scores))
      if (is.null(mat) || nrow(mat) == 0) return(out)
      idx <- match(gene_scores$Gene, rownames(mat))
      for (i in seq_along(idx)) {
        if (is.na(idx[i])) next  # gene outside the null's universe -> NA, never 1/(N+1)
        out[i] <- .perm_p_row(gene_scores[[obs_col]][i], mat[idx[i], ])
      }
      out
    }

    gene_scores$gene_caas_pperm        <- .apply_direction("gene_caas_score",            byrank[["global"]])
    gene_scores$gene_caas_pperm_top    <- .apply_direction("gene_caas_score_top_all",    byrank[["top"]])
    gene_scores$gene_caas_pperm_bottom <- .apply_direction("gene_caas_score_bottom_all", byrank[["bottom"]])

    for (pair in list(c("gene_caas_pperm",        "gene_caas_pperm_adj"),
                       c("gene_caas_pperm_top",    "gene_caas_pperm_adj_top"),
                       c("gene_caas_pperm_bottom", "gene_caas_pperm_adj_bottom"))) {
      p <- gene_scores[[pair[1]]]
      tested <- !is.na(p)
      if (any(tested)) gene_scores[[pair[2]]][tested] <- p.adjust(p[tested], method = "BH")
    }
    .n_cycles <- if (!is.null(byrank[["global"]])) ncol(byrank[["global"]]) else NA_integer_
    cat(sprintf("  gene_caas_pperm: %d/%d genes scored (N=%s cycles, floor~%.4g)\n",
                sum(!is.na(gene_scores$gene_caas_pperm)), nrow(gene_scores),
                ifelse(is.na(.n_cycles), "?", as.character(.n_cycles)),
                if (!is.na(.n_cycles) && .n_cycles > 0) 1 / (.n_cycles + 1) else NA_real_))

    if (gene_perm_pooled) {
      # Optional n-stratified pooled null: all genes' null rows within the same n_positions decile as the
      # tested gene are pooled, for higher resolution than the per-gene-row floor above. It is valid only
      # within an n-stratum: the size_adj_max null depends on n by construction (F(max)^n), so pooling
      # across strata would compare a gene with a null it was never drawn from. Off by default
      # (params.scoring_gene_perm_pooled).
      .pooled_direction <- function(obs_col, n_col, mat) {
        out <- rep(NA_real_, nrow(gene_scores))
        if (is.null(mat) || nrow(mat) == 0 || !(n_col %in% names(gene_scores))) return(out)
        n_vals <- gene_scores[[n_col]]
        strata <- suppressWarnings(dplyr::ntile(n_vals, 10))
        idx <- match(gene_scores$Gene, rownames(mat))
        for (s in sort(unique(strata[!is.na(strata)]))) {
          members <- which(strata == s & !is.na(idx))
          if (length(members) == 0) next
          pool <- as.vector(mat[idx[members], , drop = FALSE])
          pool <- pool[!is.na(pool)]
          if (length(pool) == 0) next
          for (i in members) {
            obs <- gene_scores[[obs_col]][i]
            if (!is.na(obs)) out[i] <- (sum(pool >= obs) + 1) / (length(pool) + 1)
          }
        }
        out
      }
      gene_scores$gene_caas_pperm_pooled        <- .pooled_direction("gene_caas_score",            "n_positions",        byrank[["global"]])
      gene_scores$gene_caas_pperm_pooled_top    <- .pooled_direction("gene_caas_score_top_all",    "n_positions_top",    byrank[["top"]])
      gene_scores$gene_caas_pperm_pooled_bottom <- .pooled_direction("gene_caas_score_bottom_all", "n_positions_bottom", byrank[["bottom"]])
      cat(sprintf("  gene_caas_pperm_pooled: %d/%d genes scored (n-stratified pooled null)\n",
                  sum(!is.na(gene_scores$gene_caas_pperm_pooled)), nrow(gene_scores)))
    }
  }
}

# ── 5. Correlation analysis (gene-level) ──────────────────────────────────────

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

# ── 6. Write outputs ──────────────────────────────────────────────────────────

cat("\n─── Writing outputs ───────────────────────────────────────────\n")

# Position scores. The per-scheme factors (asr_score, derived_agreement) and asr_path_score are not
# written: CAAS_score is the position-level number (the asr_path_score of a position equals it), and
# the per-scheme breakdown lives in filtered_discovery.tsv. top_species_residues and
# bottom_species_residues (alignment-based, independent of the hypothesis) are the raw-residue
# descriptors.
pos_out <- pos_scores %>%
  select(Gene, Position,
         n_schemes, any_of("scheme_set"),
         any_of(c("n_hypotheses", "participating_hypotheses",
                  "top_species_residues", "bottom_species_residues",
                  "n_top_species", "n_bottom_species")), CAAS_score,
         side,
         any_of("caas"), any_of(c("ancestral_aa", "derived_aa")),
         # p.emp: the pooled "detects AND exceeds" position p; p.adj_bh and
         # p.adj_sam: its BH and permutation-FDR adjustments (§2h).
         any_of(c("p.emp", "p.adj_bh", "p.adj_sam"))) %>%
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
    any_of(c("gene_caas_pperm", "gene_caas_pperm_top", "gene_caas_pperm_bottom",
             "gene_caas_pperm_adj", "gene_caas_pperm_adj_top", "gene_caas_pperm_adj_bottom",
             "gene_caas_pperm_pooled", "gene_caas_pperm_pooled_top", "gene_caas_pperm_pooled_bottom")),
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

# ── FCS stats table ───────────────────────────────────────────────────────────
# Read by 12.FCS_general_report.Rmd. Columns: gene; score_<ranking>, zero-floored downstream over the
# cleaned_background universe; flag_<name>, per-gene booleans. score_global, score_top and score_bottom
# are the *_all gene scores (magnitude); the top and bottom ones use only the positions of that side.
# The flags annotate the leading edge and never gate the FCS input.
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
  flag_accum_bottom = .istrue(.col(gene_scores, "accum_significant_bottom")),
  # CAAS permulation-null significance (gene_caas_pperm_adj <= p_emp_thr). It is a separate evidence
  # axis from flag_fade, flag_rer and flag_accum, which record whether the other tools flagged the gene;
  # keeping it apart lets a reader ask whether a gene is CAAS-significant independently of them.
  flag_caas_significant = .istrue(.col(gene_scores, "gene_caas_pperm_adj") <= p_emp_thr)
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
  # FCS ranks in descending order (higher = more significant) while the CCT p runs the other way, so
  # the ranking axis is -log10(accum_cct_p), derived here and not stored as a column of gene_scores.
  .afp <- suppressWarnings(as.numeric(.col(gene_scores, "accum_cct_p")))
  write_tsv(tibble(
    gene         = gene_scores$Gene,
    score_global = ifelse(is.na(.afp), NA_real_, -log10(pmax(.afp, 1e-300)))
  ), "fcs_stats_accum.tsv")
  cat("  fcs_stats_accum.tsv written (accumulation rankings, -log10 CCT p)\n")
}

# Correlations
write_tsv(corr_results, "gene_correlations.tsv")

# ── 7. Ranked gene slices ─────────────────────────────────────────────────────

dir.create("gene_lists", showWarnings = FALSE)

# 12 percentile slices (top, bottom, global × 25, 10, 5, 1%). The genes with a defined score are
# ranked in descending order and the top frac of that ranked set is kept (not of the full
# background), the same selection axis as the position slices of section 7b. They are percentile
# slices, not a significance gate: a foreground close to the whole gene universe is a poor input for
# interaction-density tests. The global slices rank on gene_caas_score, the undirected size-adjusted
# score. fade_sig is NA for them because FADE has no undirected flag, and is_fade in the export loop
# treats a missing or NA fade_sig as FALSE.
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
    # Write an empty slice so that the process always finds its declared output
    write_tsv(tibble(Gene = character(), score = numeric(), is_fade = logical(), is_rer = logical(), is_accum = logical()), file_name)
    cat(sprintf("  %s: empty (0 genes) exported\n", file_name))
    next
  }

  n_keep <- max(1, round(slice$frac * n_total))
  slice_df <- ranked %>% slice_head(n = n_keep)

  # Evidence columns: whether the other tools also call the gene significant, matched to the direction
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
      # Direction-matched like is_fade: a "top" slice checks accum_significant_top (accumulation_top_*),
      # a "bottom" slice accum_significant_bottom, a "global" slice accum_significant (accumulation_all_*).
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

# ── 7b. Ranked position slices ────────────────────────────────────────────────
# Position-level analog of the gene slices above. 14.Position_enrichment_report.Rmd and posenrich.nf read
# them (position_lists_dir), so the position percentiles are defined in one place. The direction filters
# on side (top: side == "top", bottom: side == "bottom", global: every row of position_scores.tsv); only
# positions with CAAS_score > 0 are ranked; the cut is round(frac * n_scored) after a stable sort by
# (desc score, Gene, Position).
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
    pos_out %>% dplyr::filter(side == "top")
  } else if (slice$direction == "bottom") {
    pos_out %>% dplyr::filter(side == "bottom")
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

# ── 8. Threshold enrichment ───────────────────────────────────────────────────
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

# ── 8b. Position-level FADE enrichment ────────────────────────────────────────
has_fade_site_data <- (has_fade_site_top && nrow(fade_site_top_df) > 0) ||
                      (has_fade_site_bot && nrow(fade_site_bot_df) > 0)

if (has_fade_site_data) {
  cat("\n─── Position-level FADE threshold enrichment ─────────────────\n")
  pos_enrich_rows <- list()
  k_pe            <- 1

  for (dir in c("top", "bottom")) {
    pos_subset <- pos_scores %>% filter(side == dir)
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
