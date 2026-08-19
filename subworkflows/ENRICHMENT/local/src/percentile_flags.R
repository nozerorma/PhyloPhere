#!/usr/bin/env Rscript
# =============================================================================
# percentile_flags.R - shared gene percentile-flag helper
# =============================================================================
# Defines fcs_percentile_flags() only. No library() calls: sourced by both
# fcs_enrich.R (which needs RERconverge/dplyr/etc for its own machinery) and
# 15.Comparison_report.Rmd (which doesn't), so neither picks up a dependency
# it doesn't already have.
# =============================================================================

# ── Per-gene percentile flags (directional) ──────────────────────────────────
# From score_top / score_bottom columns, derive cumulative top-X% membership
# among the NON-ZERO (signal) genes of each directional ranking:
#   pct_top25  = gene is in the top 25% of the top ranking (⊇ top10 ⊇ top5 ⊇ top1)
#   pct_top10  = gene is in the top 10% of the top ranking      (⊇ top5 ⊇ top1)
#   pct_bottom10 = gene is in the top 10% of the bottom ranking (strongly bottom)
# Returns a per-gene data.frame of logical columns (gene + pct_*). Quartet, not
# a trio, 25% is a real, tested tier elsewhere in this pipeline (posenrich's
# own Fisher test runs at char_fracs 0.25/0.10/0.05/0.01, AMI's DOMINO module
# search runs on top25/bottom25 gene lists), so leaving it out of the
# percentile-membership criteria here (and downstream in the Interesting
# Genes/Positions tables) would be inconsistent with what's actually tested.
#
# gene_lists_dir (optional): SCORING's published gene_lists/slice_<dir><pct>.tsv
# (scoring_compute.R), when supplied AND all 4 files for a direction exist,
# membership comes from there instead of re-ranking score_top/score_bottom
# here, so this and every other report/module (AMI already does) agree on
# gene percentile membership by construction. Only valid when score_top/
# score_bottom actually ARE CAAS's gene_caas_score_top_all/bottom_all (true
# for CAAS's own report and for 15.Comparison_report.Rmd's gene_attr_df,
# NEVER for RER's report, which reuses this same function with its own
# ranking under the same column names), callers must only pass
# gene_lists_dir when that holds. Falls back to re-ranking (all 4 files for a
# direction, not a mix) when gene_lists_dir is NULL or a file is missing.
fcs_percentile_flags <- function(stats, breaks = c(0.25, 0.10, 0.05, 0.01), gene_lists_dir = NULL) {
  pct_of <- function(v) {
    ok <- !is.na(v) & v > 0
    p <- rep(NA_real_, length(v))
    if (any(ok)) p[ok] <- 1 - (rank(v[ok], ties.method = "max") - 1) / sum(ok)
    p
  }
  .slice_genes <- function(direction, pct) {
    if (is.null(gene_lists_dir)) return(NULL)
    f <- file.path(gene_lists_dir, sprintf("slice_%s%d.tsv", direction, pct))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !"Gene" %in% names(d)) return(NULL)
    unique(as.character(d$Gene))
  }
  .published_quartet <- function(direction) {
    g25 <- .slice_genes(direction, 25); g10 <- .slice_genes(direction, 10)
    g5  <- .slice_genes(direction, 5);  g1  <- .slice_genes(direction, 1)
    if (is.null(g25) || is.null(g10) || is.null(g5) || is.null(g1)) return(NULL)
    list(g25 = g25, g10 = g10, g5 = g5, g1 = g1)
  }

  out <- data.frame(gene = stats$gene, stringsAsFactors = FALSE)
  if ("score_top" %in% names(stats)) {
    pub <- .published_quartet("top")
    if (!is.null(pub)) {
      out$pct_top25 <- out$gene %in% pub$g25
      out$pct_top10 <- out$gene %in% pub$g10
      out$pct_top5  <- out$gene %in% pub$g5
      out$pct_top1  <- out$gene %in% pub$g1
    } else {
      pt <- pct_of(stats$score_top)
      out$pct_top25 <- !is.na(pt) & pt <= breaks[1]
      out$pct_top10 <- !is.na(pt) & pt <= breaks[2]
      out$pct_top5  <- !is.na(pt) & pt <= breaks[3]
      out$pct_top1  <- !is.na(pt) & pt <= breaks[4]
    }
  }
  if ("score_bottom" %in% names(stats)) {
    pub <- .published_quartet("bottom")
    if (!is.null(pub)) {
      out$pct_bottom25 <- out$gene %in% pub$g25
      out$pct_bottom10 <- out$gene %in% pub$g10
      out$pct_bottom5  <- out$gene %in% pub$g5
      out$pct_bottom1  <- out$gene %in% pub$g1
    } else {
      pb <- pct_of(stats$score_bottom)
      out$pct_bottom25 <- !is.na(pb) & pb <= breaks[1]
      out$pct_bottom10 <- !is.na(pb) & pb <= breaks[2]
      out$pct_bottom5  <- !is.na(pb) & pb <= breaks[3]
      out$pct_bottom1  <- !is.na(pb) & pb <= breaks[4]
    }
  }
  out
}
