#!/usr/bin/env Rscript
# fcs_enrich.R — Functional class scoring (FCS) core: three gene-set tests per ranking.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: fcs_compute.R (FCS_COMPUTE_BATCHED process, fcs.nf) and 12.FCS_general_report.Rmd.
# Defines functions only.
#
# Tests the gene sets of GMT files against a gene ranking with three complementary tests:
#   1. Wilcoxon-AUC (RERconverge::fastwilcoxGMTall): rank shift of the set.
#   2. Lachenbruch two-part: prevalence of nonzero scorers (Fisher) plus magnitude among
#      them (Wilcoxon), combined as chi-square with 2 df. Safe for zero-inflated scores.
#   3. Path-sum permulation: observed sum of the scores of the set against a null; NES is
#      the z-score of the observed sum. Safe for zero-inflated scores.
# Each test casts one vote (FDR below its threshold plus a direction check, and the
# permulation p-value gate where a null exists). A term is "Hard evidence" with 3 votes,
# "Supported" with 2, and "Exploratory" with 1: "(relative)" for a lone Wilcoxon or
# Lachenbruch pass, "(phylogenetic)" for a lone permulation pass, which makes no claim
# relative to the other genes. Tests 2 and 3 run only on non-negative (zero-floored)
# rankings; signed rankings (RER, two-sided) use the Wilcoxon test alone.
#
# Design:
#   * The universe is the full tested background (cleaned_background). Genes without
#     signal are floored to 0, never dropped or penalized, as in RERconverge's
#     accelerating and decelerating rankings.
#   * Direction by membership: a "top" ranking keeps the score of the top and both
#     genes and floors the rest to 0; "bottom" keeps bottom and both. A gene of both
#     directions contributes its full score to each.
#   * Significance is not used to filter the input. Callers attach it to the leading-edge
#     genes as annotation (fcs_annotate_leading_edge).
#   * Multiple testing is BH within each GMT (each database is its own family of
#     hypotheses); fastwilcoxGMTall already adjusts per GMT.
# =============================================================================

suppressPackageStartupMessages({
  library(RERconverge)   # fastwilcoxGMT, read.gmt
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(parallel)
  library(Matrix)        # sparse membership matrix + %*% for the vectorized null
})

# fcs_percentile_flags() lives in percentile_flags.R, which 15.Comparison_report.Rmd also
# sources directly to avoid the library() calls above. It is sourced here unless a caller
# already did. The lookup uses paths relative to the working directory, as the callers
# stage this file (see src_candidates in 12.FCS_general_report.Rmd); sys.frame()$ofile is
# not reliable when source() runs inside a chunk of rmarkdown::render().
if (!exists("fcs_percentile_flags", mode = "function")) {
  .pflags_candidates <- c("src/percentile_flags.R", "percentile_flags.R")
  .pflags_path <- .pflags_candidates[file.exists(.pflags_candidates)][1]
  if (is.na(.pflags_path)) stop("percentile_flags.R not found under src/ or next to fcs_enrich.R")
  source(.pflags_path)
}


# ── GMT loading ───────────────────────────────────────────────────────────────

# Database label of a GMT file, without extension or version:
#   MSigDB     c6.all.v2026.1.Hs.symbols.gmt  -> c6
#   WebGestalt pathway_KEGG.symbols.gmt        -> pathway_KEGG
fcs_db_name <- function(path) {
  x <- basename(path)
  x <- sub("\\.gmt$", "", x)
  x <- sub("\\.v[0-9.]+\\.Hs\\.symbols$", "", x)
  x <- sub("\\.symbols$", "", x)
  x <- sub("\\.all$", "", x)
  x
}

# Table (database, pathway, description) from the GMT files.
# GMT line: <setID>\t<description>\t<gene1>\t<gene2>... WebGestalt GMTs carry the term name
# in column 2 and MSigDB a URL, so a URL or empty description falls back to the set ID.
fcs_load_descriptions <- function(gmt_dir) {
  files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
  rows <- list()
  for (f in files) {
    db <- fcs_db_name(f)
    sp <- strsplit(readLines(f, warn = FALSE), "\t", fixed = TRUE)
    sp <- sp[lengths(sp) >= 1]
    id   <- vapply(sp, function(x) x[1], character(1))
    desc <- vapply(sp, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    desc <- ifelse(is.na(desc) | !nzchar(desc) | grepl("^https?://", desc), id, desc)
    # TRANSFAC-style motif IDs (e.g. V$HNF4_01, TGTTTGY_V$HNF3_Q6) have no separate
    # description: the transcription factor symbol is shown instead, only when desc repeats id.
    tf <- ifelse(grepl("V\\$", id), sub(".*V\\$([A-Za-z0-9]+).*", "\\1", id), NA_character_)
    desc <- ifelse(!is.na(tf) & nzchar(tf) & desc == id,
                   paste0(tf, " - TF target motif"), desc)
    rows[[db]] <- data.frame(database = db, pathway = id, description = desc,
                             stringsAsFactors = FALSE)
  }
  dplyr::distinct(dplyr::bind_rows(rows))
}

# Every *.gmt of gmt_dir as a named list of RERconverge gmt objects ($genesets and
# $geneset.names, the structure fastwilcoxGMT expects). Files that fail to parse are skipped.
fcs_load_gmts <- function(gmt_dir) {
  files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
  if (length(files) == 0) stop(sprintf("No GMT files found in: %s", gmt_dir))
  gmts <- list()
  for (f in files) {
    obj <- tryCatch(RERconverge::read.gmt(f), error = function(e) {
      message(sprintf("  read.gmt failed for %s: %s", basename(f), e$message)); NULL
    })
    if (!is.null(obj) && length(obj$genesets) > 0) gmts[[fcs_db_name(f)]] <- obj
  }
  if (length(gmts) == 0) stop("No GMT files could be parsed.")
  gmts
}

# ── Ranking construction ──────────────────────────────────────────────────────

# scores: named numeric vector (gene -> score) of the signal genes only.
# universe: all tested genes (cleaned_background).
# Returns a named numeric vector over union(universe, names(scores)); genes without a
# score, or with NA, get `floor` (0). The caller handles direction by passing only the
# scores of the directional subset (top and both, or bottom and both).
fcs_build_vals <- function(scores, universe, floor = 0) {
  scores <- scores[!is.na(scores)]
  genes  <- union(universe, names(scores))
  vals   <- setNames(rep(floor, length(genes)), genes)
  vals[names(scores)] <- as.numeric(scores)
  vals
}

# ── Wilcoxon-AUC test ─────────────────────────────────────────────────────────

# Wilcoxon-AUC test of one ranking over all GMTs.
# vals: named numeric vector (full universe, zero-floored).
# gmts: named list of RERconverge gmt objects.
# num_g, max_g: minimum and maximum genes per set (max_g 0 or Inf = no limit).
# alternative: "greater" or "two.sided" (see fcs_alternative).
# Returns a tibble: database, pathway, stat (AUC - 0.5), pval, p.adj (BH per GMT),
#   num.genes, gene.vals (leading-edge members as "gene:rank").
# fastwilcoxGMT is internal to RERconverge; the exported fastwilcoxGMTall loops it over a
# named list of GMTs and adjusts per GMT, so it is called once and the databases stitched.
fcs_run_ranking <- function(vals, gmts, num_g = 10, max_g = 500, alternative = "two.sided") {
  if (!is.null(max_g) && is.finite(max_g) && max_g > 0) {
    gmts_to_run <- lapply(gmts, function(gmt) {
      gs <- gmt$genesets
      if (is.null(names(gs)) && !is.null(gmt$geneset.names)) names(gs) <- gmt$geneset.names
      keep <- vapply(gs, function(set) length(intersect(set, names(vals))) <= max_g, logical(1))
      list(
        genesets = gs[keep],
        geneset.names = if (!is.null(gmt$geneset.names)) gmt$geneset.names[keep] else names(gs)[keep]
      )
    })
  } else {
    gmts_to_run <- gmts
  }

  reslist <- tryCatch(
    RERconverge::fastwilcoxGMTall(vals, gmts_to_run, outputGeneVals = TRUE, num.g = num_g),
    error = function(e) { message(sprintf("  fastwilcoxGMTall failed: %s", e$message)); NULL }
  )
  if (is.null(reslist) || length(reslist) == 0) return(tibble::tibble())
  vals <- vals[!is.na(vals)]
  out <- list()
  for (db in names(reslist)) {
    res <- reslist[[db]]
    if (is.null(res) || nrow(res) == 0) next
    if (!is.null(max_g) && is.finite(max_g) && max_g > 0 && "num.genes" %in% names(res)) {
      res <- res[res$num.genes <= max_g, , drop = FALSE]
      if (nrow(res) == 0) next
    }
    # The p-value of fastwilcoxGMT (simpleAUCgenesRanks) is always two-sided. For
    # one-sided magnitude rankings (non-negative, zero-floored: CAAS, FADE and
    # accumulation scores, RER accelerating and decelerating) the analytic p-value is
    # recomputed one-sided from the AUC, so a depleted set (stat < 0) is not flagged,
    # and BH is applied again per GMT. The background is the annotated genes of the
    # GMT, as in fastwilcoxGMT.
    if (alternative == "greater") {
      gmt  <- gmts_to_run[[db]]
      n_db <- length(intersect(unique(unlist(gmt$genesets)), names(vals)))
      n1   <- res$num.genes
      n2   <- n_db - n1
      U    <- (res$stat + 0.5) * n1 * n2
      mu   <- n1 * n2 / 2
      sdv  <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
      res$pval  <- pnorm(U, mu, sdv, lower.tail = FALSE)
      res$p.adj <- p.adjust(res$pval, method = "BH")
    }
    out[[db]] <- tibble::as_tibble(res, rownames = "pathway") %>%
      dplyr::mutate(database = db)
  }
  if (length(out) == 0) return(tibble::tibble())
  dplyr::bind_rows(out) %>%
    dplyr::relocate(database, pathway) %>%
    dplyr::arrange(p.adj, dplyr::desc(abs(stat)))
}

# ── Test sidedness ────────────────────────────────────────────────────────────

# A ranking with negative values (a signed statistic such as sign(Rho) * -log10(P)) is
# two-sided; a non-negative magnitude ranking is one-sided "greater". It is read from
# the values, so callers need not pass it.
fcs_alternative <- function(vals) if (any(vals < 0, na.rm = TRUE)) "two.sided" else "greater"

# ── Lachenbruch two-part test ─────────────────────────────────────────────────

# Enrichment in a zero-inflated distribution, as the combination of two parts:
#   Part 1: one-sided Fisher exact test of the 2x2 prevalence table (score > 0 or = 0).
#   Part 2: one-sided Wilcoxon test on the positive-scoring genes only (magnitude).
# Each p-value (floored at 1e-15) becomes chi-square with 1 df; their sum is chi-square
# with 2 df, which gives the combined p-value. BH per GMT. Only for non-negative,
# zero-floored vals.
# Returns per pathway: lach_pval, lach_p.adj, lach_chi_binary (Part 1), lach_chi_nonzero
#   (Part 2; 0 when Part 2 cannot run), lach_chi_total and lach_frac_magnitude
#   (lach_chi_nonzero / lach_chi_total; high when magnitude and not prevalence drives it).
fcs_run_lachenbruch <- function(vals, gmts, num_g = 10, max_g = 500) {
  vals[is.na(vals)] <- 0
  universe   <- names(vals)
  hits       <- names(vals[vals > 0])
  pos_scores <- vals[vals > 0]

  out <- list()
  for (db in names(gmts)) {
    gmt <- gmts[[db]]
    gs  <- gmt$genesets
    if (is.null(names(gs))) names(gs) <- gmt$geneset.names

    rows <- list()
    for (pname in names(gs)) {
      p_genes <- intersect(gs[[pname]], universe)
      n1 <- length(p_genes)
      if (n1 < num_g || (!is.null(max_g) && is.finite(max_g) && max_g > 0 && n1 > max_g)) next

      k1 <- length(intersect(p_genes, hits))
      k0 <- n1 - k1
      b1 <- length(hits) - k1
      b0 <- (length(universe) - n1) - b1

      mat  <- matrix(c(k1, k0, b1, b0), nrow = 2, byrow = TRUE)
      ft   <- tryCatch(fisher.test(mat, alternative = "greater"), error = function(e) NULL)
      if (is.null(ft)) next
      p1   <- max(ft$p.value, 1e-15)
      chi1 <- qchisq(p1, df = 1, lower.tail = FALSE)

      p_pos  <- intersect(p_genes, names(pos_scores))
      n1_pos <- length(p_pos)
      chi2   <- 0
      if (n1_pos >= 2 && (length(pos_scores) - n1_pos) >= 2) {
        bg_pos <- setdiff(names(pos_scores), p_pos)
        wt     <- tryCatch(
          wilcox.test(pos_scores[p_pos], pos_scores[bg_pos], alternative = "greater"),
          error = function(e) NULL)
        if (!is.null(wt)) {
          p2   <- max(wt$p.value, 1e-15)
          chi2 <- qchisq(p2, df = 1, lower.tail = FALSE)
        }
      }

      chi_total <- chi1 + chi2
      rows[[pname]] <- tibble::tibble(
        database            = db,
        pathway             = pname,
        lach_pval           = pchisq(chi_total, df = 2, lower.tail = FALSE),
        lach_chi_binary     = chi1,
        lach_chi_nonzero    = chi2,
        lach_chi_total      = chi_total,
        lach_frac_magnitude = if (chi_total > 0) chi2 / chi_total else NA_real_
      )
    }
    if (length(rows) == 0) next
    db_df <- dplyr::bind_rows(rows)
    db_df$lach_p.adj <- p.adjust(db_df$lach_pval, method = "BH")
    out[[db]] <- db_df
  }
  if (length(out) == 0) return(tibble::tibble())
  dplyr::bind_rows(out)
}

# ── Path-sum permulation ──────────────────────────────────────────────────────

# Score accumulation: the observed sum of the scores of a set against a null of sums.
# Zero-floored genes add 0 and do not distort the sums. BH per GMT.
# NES = (obs_sum - null_mean) / null_sd; NES > 0 means enrichment.
# p = (1 + number of null sums >= obs_sum) / (n_perms + 1), so its floor is 1 / (n_perms + 1).
# Only for non-negative vals.
# null_mat: genes x N permulation null (CAAS or RER). Without it, vals are shuffled
# n_perms times. fcs_run_all always passes a null_mat.
# Returns per pathway: perm_pval, perm_p.adj, perm_nes.
fcs_run_permulation <- function(vals, gmts, num_g = 10, max_g = 500, n_perms = 2000, seed = 1998,
                                null_mat = NULL) {
  vals[is.na(vals)] <- 0
  genes_all <- names(vals)
  N_genes   <- length(genes_all)

  if (!is.null(null_mat)) {
    # The shared CAAS or RER permulation null, the one behind the Wilcoxon p.perm, is
    # preferred over a label shuffle. Its rows are aligned to the genes of this ranking:
    # a gene of `vals` missing from null_mat gets an all-0 row (an unscored gene adds
    # nothing to the observed sum either), and a gene of null_mat missing from `vals`
    # is dropped.
    common   <- intersect(genes_all, rownames(null_mat))
    perm_mat <- matrix(0, nrow = N_genes, ncol = ncol(null_mat),
                       dimnames = list(genes_all, NULL))
    perm_mat[common, ] <- null_mat[common, , drop = FALSE]
    n_perms <- ncol(perm_mat)   # the size of the null itself, not the n_perms argument
  } else {
    set.seed(seed)
    perm_mat <- vapply(seq_len(n_perms), function(i) sample(vals), numeric(N_genes))
    rownames(perm_mat) <- genes_all
  }

  out <- list()
  for (db in names(gmts)) {
    gmt <- gmts[[db]]
    gs  <- gmt$genesets
    if (is.null(names(gs))) names(gs) <- gmt$geneset.names

    valid <- sapply(gs, function(g) {
      n <- length(intersect(g, genes_all))
      n >= num_g && (is.null(max_g) || !is.finite(max_g) || max_g <= 0 || n <= max_g)
    })
    gs_v  <- gs[valid]
    if (length(gs_v) == 0) next
    pnames <- names(gs_v)

    gene_idx <- setNames(seq_along(genes_all), genes_all)
    ri <- integer(0); ci <- integer(0)
    for (i in seq_along(gs_v)) {
      g_in <- intersect(gs_v[[i]], genes_all)
      if (length(g_in)) { ri <- c(ri, rep(i, length(g_in))); ci <- c(ci, gene_idx[g_in]) }
    }
    M <- Matrix::sparseMatrix(i = ri, j = ci, x = 1,
                              dims = c(length(gs_v), N_genes),
                              dimnames = list(pnames, genes_all))

    obs_sums  <- as.numeric(M %*% vals)
    null_sums <- as.matrix(M %*% perm_mat)   # pathways × n_perms
    null_mu   <- rowMeans(null_sums)
    null_sd   <- apply(null_sums, 1, sd)
    nes       <- (obs_sums - null_mu) / ifelse(null_sd == 0, 1, null_sd)
    pvals     <- (rowSums(null_sums >= obs_sums) + 1) / (n_perms + 1)

    db_df <- tibble::tibble(
      database   = db,
      pathway    = pnames,
      perm_pval  = pvals,
      perm_p.adj = p.adjust(pvals, method = "BH"),
      perm_nes   = nes
    )
    out[[db]] <- db_df
  }
  if (length(out) == 0) return(tibble::tibble())
  dplyr::bind_rows(out)
}

# ── Progress logging ──────────────────────────────────────────────────────────

# Writes a timestamped line to stderr and appends it to fcs_progress.log in the working
# directory. knitr buffers the output of a chunk, so the file is the channel to follow
# a run with `tail -f`.
fcs_progress <- function(msg, file = "fcs_progress.log") {
  line <- sprintf("[FCS %s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  message(line)
  try(cat(line, "\n", sep = "", file = file, append = TRUE), silent = TRUE)
}

# ── Vectorized permulation null ───────────────────────────────────────────────

# Per-column ranks (average ties); NA gets 0 so absent genes drop out of the rank sum, as
# fastwilcoxGMT does by removing NA values before rank().
fcs_colranks <- function(m) {
  if (!anyNA(m)) {
    r <- apply(m, 2, rank, ties.method = "average")
  } else {
    r <- apply(m, 2, function(v) { x <- numeric(length(v)); o <- !is.na(v)
                                   x[o] <- rank(v[o], ties.method = "average"); x })
  }
  if (is.null(dim(r))) r <- matrix(r, nrow = nrow(m))
  rownames(r) <- rownames(m)
  r
}

# Sparse membership matrix (sets x genes) over a fixed gene space.
fcs_membership_matrix <- function(genesets_named, set_names, genes) {
  gi <- setNames(seq_along(genes), genes)
  ii <- integer(0); jj <- integer(0)
  for (s in seq_along(set_names)) {
    g <- genesets_named[[ set_names[s] ]]
    idx <- gi[g]; idx <- idx[!is.na(idx)]
    if (length(idx)) { ii <- c(ii, rep.int(s, length(idx))); jj <- c(jj, as.integer(idx)) }
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = 1,
                       dims = c(length(set_names), length(genes)),
                       dimnames = list(set_names, genes))
}

# Null of the Wilcoxon-AUC statistic: fastwilcoxGMTall(corStat[, j], gmts)$stat for every
# permulation column j, in one sparse matrix product per GMT instead of N calls.
# It follows fastwilcoxGMT: the background is the annotated genes of the GMT given (the
# observed side drops the sets over max_g first, and so does this), ranks are computed
# per GMT and per column with average ties, and a set that fails num.g (or has <= 2
# background genes) in a column is NA for that column.
# Returns a named list db -> (sets x N) matrix of AUC - 0.5, rows aligned to
# rownames(realenrich[[db]]).
fcs_null_enrichstat_vectorized <- function(corStat, gmts, realenrich, num_g = 10, max_g = 500) {
  enrichStat <- list()
  genes_all  <- rownames(corStat)
  for (db in names(realenrich)) {
    gmt <- gmts[[db]]
    set_names <- rownames(realenrich[[db]])
    if (is.null(gmt) || length(set_names) == 0) {
      enrichStat[[db]] <- matrix(NA_real_, nrow = length(set_names), ncol = ncol(corStat),
                                 dimnames = list(set_names, NULL))
      next
    }
    gs <- gmt$genesets; names(gs) <- gmt$geneset.names
    # The observed side (fcs_run_ranking) gives fastwilcoxGMT the GMT without its sets
    # over max_g, and fastwilcoxGMT takes the genes of the sets it receives as the
    # background, so genes held only by dropped sets are not in it. The null uses the same
    # background; a dropped set has no member left here and is NA below.
    if (!is.null(max_g) && is.finite(max_g) && max_g > 0) {
      gs <- gs[vapply(gs, function(set) length(intersect(set, genes_all)) <= max_g, logical(1))]
    }
    genes_db <- intersect(unique(unlist(gs)), genes_all)
    if (length(genes_db) < 3) {
      enrichStat[[db]] <- matrix(NA_real_, nrow = length(set_names), ncol = ncol(corStat),
                                 dimnames = list(set_names, NULL))
      next
    }
    M   <- fcs_membership_matrix(gs, set_names, genes_db)  # sets x genes_db
    sub <- corStat[genes_db, , drop = FALSE]               # genes_db x N
    Rk  <- fcs_colranks(sub)                               # NA -> 0
    ranksum <- as.matrix(M %*% Rk)                         # sets x N
    if (!anyNA(sub)) {
      n1   <- as.numeric(Matrix::rowSums(M))               # per set, constant over cols
      n2   <- length(genes_db) - n1
      U    <- ranksum - (n1 * (n1 + 1) / 2)                # n1/U recycle down columns
      stat <- (U / (n1 * n2)) - 0.5
      stat[(n1 < num_g) | (!is.null(max_g) & is.finite(max_g) & max_g > 0 & n1 > max_g) | (n2 <= 2), ] <- NA_real_
    } else {
      notNA <- !is.na(sub)
      n1    <- as.matrix(M %*% (notNA * 1.0))              # sets x N
      ntot  <- matrix(colSums(notNA), nrow = nrow(M), ncol = ncol(sub), byrow = TRUE)
      n2    <- ntot - n1
      U     <- ranksum - (n1 * (n1 + 1) / 2)
      stat  <- (U / (n1 * n2)) - 0.5
      stat[(n1 < num_g) | (!is.null(max_g) & is.finite(max_g) & max_g > 0 & n1 > max_g) | (n2 <= 2)] <- NA_real_
    }
    rownames(stat) <- set_names
    enrichStat[[db]] <- stat
  }
  enrichStat
}

# Empirical permulation p-value per set, with a pseudo-count, by ranking sidedness:
#   greater   -> (#{null >= obs} + 1) / (N_valid + 1)       [magnitude rankings]
#   two.sided -> (#{|null| >= |obs|} + 1) / (N_valid + 1)   [signed rankings]
# N_valid counts the non-NA null columns of the set.
fcs_permpvalenrich_vectorized <- function(realenrich, enrichStat, alternative = "two.sided") {
  out <- list()
  for (db in names(realenrich)) {
    null <- enrichStat[[db]]
    if (is.null(null) || nrow(null) == 0) next
    obs <- realenrich[[db]]$stat[match(rownames(null), rownames(realenrich[[db]]))]
    count <- if (alternative == "greater") rowSums(null >= obs, na.rm = TRUE)
             else                          rowSums(abs(null) >= abs(obs), na.rm = TRUE)
    denom <- rowSums(!is.na(null))
    p <- (count + 1) / (denom + 1)
    p[is.na(obs)] <- NA_real_
    names(p) <- rownames(null)
    out[[db]] <- p
  }
  out
}

# ── Lachenbruch null, Part 1 ──────────────────────────────────────────────────

# Null of the prevalence part (Fisher). A one-sided "greater" Fisher exact test of a 2x2
# table equals the upper tail of the hypergeometric distribution,
# phyper(k1 - 1, m = k1 + b1, n = k0 + b0, k = k1 + k0, lower.tail = FALSE), so one
# phyper() call per database vectorizes it exactly over every (set, null column) cell.
#
# The universe is the whole ranking (rownames(corStat_rk)), as in the observed side of
# fcs_run_lachenbruch (universe = names(vals)), and not the GMT-annotated genes that
# fcs_null_enrichstat_vectorized uses for the Wilcoxon background. Each test's null must
# share the background of its own observed side.
#
# Returns a named list db -> (sets x N) matrix of chi-square values (1 df), NA where the
# set fails num_g or max_g, as in fcs_null_enrichstat_vectorized.
fcs_null_lachenbruch_binary_vectorized <- function(corStat_rk, gmts, realenrich, num_g = 10, max_g = 500) {
  chi1 <- list()
  genes_all  <- rownames(corStat_rk)
  N          <- ncol(corStat_rk)
  n_universe <- length(genes_all)

  hitmat <- (corStat_rk > 0) * 1.0
  hitmat[is.na(hitmat)] <- 0        # NA scores treated as non-hit (zero-floor convention)
  m_col <- Matrix::colSums(hitmat)  # total hits in universe, per column

  for (db in names(realenrich)) {
    gmt <- gmts[[db]]
    set_names <- rownames(realenrich[[db]])
    if (is.null(gmt) || length(set_names) == 0) {
      chi1[[db]] <- matrix(NA_real_, nrow = length(set_names), ncol = N,
                           dimnames = list(set_names, NULL))
      next
    }
    gs <- gmt$genesets; names(gs) <- gmt$geneset.names
    M  <- fcs_membership_matrix(gs, set_names, genes_all)   # sets x genes_all

    n1 <- as.numeric(Matrix::rowSums(M))          # geneset size within universe, constant per set
    k1 <- as.matrix(M %*% hitmat)                 # sets x N -- hits within set, per column
    m  <- matrix(m_col, nrow = nrow(M), ncol = N, byrow = TRUE)
    n  <- n_universe - m
    k  <- matrix(n1, nrow = nrow(M), ncol = N)

    p1 <- matrix(
      phyper(as.vector(k1) - 1, m = as.vector(m), n = as.vector(n), k = as.vector(k),
             lower.tail = FALSE),
      nrow = nrow(M), ncol = N, dimnames = list(set_names, NULL))
    p1  <- pmax(p1, 1e-15)                         # same floor as Part 1 of fcs_run_lachenbruch
    chi <- qchisq(p1, df = 1, lower.tail = FALSE)
    fail <- (n1 < num_g) | (!is.null(max_g) & is.finite(max_g) & max_g > 0 & n1 > max_g) | (n1 == 0)
    chi[fail, ] <- NA_real_
    chi1[[db]] <- chi
  }
  chi1
}

# ── Lachenbruch null, Part 2 ──────────────────────────────────────────────────

# Null of the magnitude part (Wilcoxon on the positive-scoring genes). Which genes are
# positive varies by null column (a gene can be 0 in one permuted labeling and positive
# in another), the case that the NA branch of fcs_null_enrichstat_vectorized handles by
# masking the other genes to NA before ranking. That function is called with a copy of
# corStat_rk restricted to its positive entries.
#
# Its num_g is fixed at 2, not the caller's, because Part 2 of fcs_run_lachenbruch
# requires only 2 positive members (n1_pos >= 2), independently of and below the num_g
# of Part 1. The caller's num_g would set to NA sets with 2 to num_g - 1 positive
# scorers that the observed side accepts.
#
# Returns a named list db -> (sets x N) matrix of chi-square values (1 df) from the
# normal approximation of the Wilcoxon U statistic. It is 0, not NA, where Part 2 cannot
# run, as in fcs_run_lachenbruch, so that chi_total is defined whenever Part 1 is.
fcs_null_lachenbruch_magnitude_vectorized <- function(corStat_rk, gmts, realenrich, max_g = 500) {
  corStat_pos <- corStat_rk
  corStat_pos[corStat_rk <= 0] <- NA

  auc_pos <- fcs_null_enrichstat_vectorized(corStat_pos, gmts, realenrich, num_g = 2, max_g = max_g)

  chi2 <- list()
  for (db in names(realenrich)) {
    auc <- auc_pos[[db]]
    if (is.null(auc) || nrow(auc) == 0) { chi2[[db]] <- auc; next }

    gmt <- gmts[[db]]
    set_names <- rownames(realenrich[[db]])
    gs <- gmt$genesets; names(gs) <- gmt$geneset.names
    genes_db <- intersect(unique(unlist(gs)), rownames(corStat_pos))
    M_pos    <- fcs_membership_matrix(gs, set_names, genes_db)
    sub_pos  <- corStat_pos[genes_db, , drop = FALSE]

    # n1 and n2 per column, as fcs_null_enrichstat_vectorized derives them internally from
    # the NA pattern of sub_pos (that function returns only the statistic).
    notNA <- !is.na(sub_pos)
    n1n   <- as.matrix(M_pos %*% (notNA * 1.0))
    ntot  <- matrix(colSums(notNA), nrow = nrow(M_pos), ncol = ncol(sub_pos), byrow = TRUE)
    n2n   <- ntot - n1n

    U    <- (auc + 0.5) * n1n * n2n
    mu   <- n1n * n2n / 2
    sdv  <- sqrt(n1n * n2n * (n1n + n2n + 1) / 12)
    # Continuity correction of 0.5, as in the normal approximation of wilcox.test() for
    # alternative = "greater"; without it the observed and null values diverge, even at
    # large n1 and n2, because qchisq is steep at small p.
    z    <- (U - mu - 0.5) / sdv
    p2   <- pmax(pnorm(z, lower.tail = FALSE), 1e-15)   # same floor as Part 2 of fcs_run_lachenbruch
    chi  <- qchisq(p2, df = 1, lower.tail = FALSE)
    chi[is.na(auc)] <- 0   # Part 2 unavailable: 0, never NA, as on the observed side
    chi2[[db]] <- chi
  }
  chi2
}

# ── Lachenbruch empirical p-value ─────────────────────────────────────────────

# Empirical lach_p.perm of every row of lach_rk: the null chi-square of both parts is
# summed per column and compared with the observed lach_chi_total through
# fcs_permpvalenrich_vectorized. Returns a vector aligned to the rows of lach_rk, next
# to the analytic lach_pval and lach_p.adj.
fcs_compute_lach_p_perm <- function(lach_rk, corStat_rk, gmts, num_g = 10, max_g = 500) {
  realenrich <- list()
  for (db in unique(lach_rk$database)) {
    db_df <- lach_rk[lach_rk$database == db, , drop = FALSE]
    realenrich[[db]] <- data.frame(stat = db_df$lach_chi_total, row.names = db_df$pathway)
  }

  chi1_null <- fcs_null_lachenbruch_binary_vectorized(corStat_rk, gmts, realenrich, num_g = num_g, max_g = max_g)
  chi2_null <- fcs_null_lachenbruch_magnitude_vectorized(corStat_rk, gmts, realenrich, max_g = max_g)
  chi_total_null <- setNames(
    lapply(names(realenrich), function(db) chi1_null[[db]] + chi2_null[[db]]),
    names(realenrich))

  ppv <- fcs_permpvalenrich_vectorized(realenrich, chi_total_null, alternative = "greater")

  out <- rep(NA_real_, nrow(lach_rk))
  for (db in names(ppv)) {
    idx <- which(lach_rk$database == db)
    out[idx] <- ppv[[db]][lach_rk$pathway[idx]]
  }
  out
}

# ── Null per ranking ──────────────────────────────────────────────────────────

# Null matrix (genes x N) of ranking `rk`, or NULL when there is none. Both loops of
# fcs_run_all use the same null for a ranking, so it is resolved once here (corStat_byrk).
# CAAS supplies one matrix per direction (corStat_byrank). RER shares one base matrix
# (base_corStat) and derives the accelerating and decelerating rankings from the
# sign of base_corRho.
fcs_resolve_corstat_rk <- function(rk, corStat_byrank, base_corStat, base_corRho) {
  if (!is.null(corStat_byrank) && !is.null(corStat_byrank[[rk]])) {
    return(as.matrix(corStat_byrank[[rk]]))
  }
  if (!is.null(base_corStat)) {
    corStat_rk <- base_corStat
    if (rk == "accelerating" && !is.null(base_corRho)) {
      corStat_rk <- ifelse(base_corRho > 0,  base_corStat, 0)
    } else if (rk == "decelerating" && !is.null(base_corRho)) {
      corStat_rk <- ifelse(base_corRho < 0, -base_corStat, 0)
    }
    rownames(corStat_rk) <- rownames(base_corStat)
    return(corStat_rk)
  }
  NULL
}

# ── Result table ──────────────────────────────────────────────────────────────

# Column types of the FCS result table (fcs_enrich_merged.tsv), for the readr::read_tsv
# of the row-concatenated batch outputs (the enrich_file path of 12.FCS_general_report.Rmd).
# A header-only TSV, which a small gene universe can produce, gives readr nothing to
# infer from and every column would be read as logical; this spec keeps the types stable
# whatever the number of rows.
fcs_enrich_col_types <- function() {
  readr::cols(
    ranking = readr::col_character(), database = readr::col_character(),
    pathway = readr::col_character(), stat = readr::col_double(),
    pval = readr::col_double(), p.adj = readr::col_double(), p.perm = readr::col_double(),
    num.genes = readr::col_double(), gene.vals = readr::col_character(),
    lach_pval = readr::col_double(), lach_p.adj = readr::col_double(),
    lach_chi_binary = readr::col_double(), lach_chi_nonzero = readr::col_double(),
    lach_chi_total = readr::col_double(), lach_frac_magnitude = readr::col_double(),
    lach_p.perm = readr::col_double(),
    perm_pval = readr::col_double(), perm_p.adj = readr::col_double(), perm_nes = readr::col_double(),
    sig_wilcoxon = readr::col_logical(), sig_lachenbruch = readr::col_logical(),
    sig_permulation = readr::col_logical(),
    evidence_count = readr::col_integer(), evidence_label = readr::col_character(),
    .default = readr::col_guess()
  )
}

# ── Evidence classification ───────────────────────────────────────────────────

# Evidence gates and labels. The FDR threshold of each test is set in conf/enrichment.config
# (fdr_wilcoxon, fdr_lachenbruch, fdr_permsum).
#   sig_wilcoxon: FDR gate, direction (stat > 0) and, when p.perm exists, the permulation
#     gate p.perm < p_perm_thr. p.perm is NA without a perms file or null, and the gate is skipped.
#   sig_lachenbruch: FDR gate and, when lach_p.perm exists, the same permulation gate on the shared null.
#   sig_permulation: FDR gate and direction (NES > 0) against the shared null; without one
#     its columns are NA (no label shuffle stands in for it).
# A ranking without a permulation null (no perms file, an empty null, a stale one) has no
# phylogenetic gate: its sig_permulation cannot pass, and a row that passes the gates that
# exist is "Exploratory (relative)", never "Supported", "Hard evidence" or "(phylogenetic)".
# no_null_rankings lists those rankings.
fcs_classify_evidence <- function(enrich_df, fdr_wilcoxon, fdr_lachenbruch, fdr_permsum, p_perm_thr, no_null_rankings = character(0)) {
  enrich_df %>%
    dplyr::mutate(
      sig_wilcoxon    = !is.na(p.adj)       & p.adj       < fdr_wilcoxon &
                        (is.na(p.perm) | p.perm < p_perm_thr) &
                        !is.na(stat)  & stat > 0,
      sig_lachenbruch = !is.na(lach_p.adj)  & lach_p.adj  < fdr_lachenbruch &
                        (is.na(lach_p.perm) | lach_p.perm < p_perm_thr),
      lacks_null      = ranking %in% no_null_rankings,
      sig_permulation = !lacks_null & !is.na(perm_p.adj)  & perm_p.adj  < fdr_permsum &
                        !is.na(perm_nes) & perm_nes > 0,
      evidence_count  = as.integer(sig_wilcoxon) +
                        as.integer(sig_lachenbruch) +
                        as.integer(sig_permulation),
      evidence_label  = dplyr::case_when(
        # no null for this ranking: no phylogenetic gate was applied, whatever else passed
        lacks_null & evidence_count >= 1L ~ "Exploratory (relative)",
        lacks_null                        ~ "Not significant",
        evidence_count == 3L ~ "Hard evidence",
        evidence_count == 2L ~ "Supported",
        # Permulation has no relative-to-other-genes requirement (Wilcoxon and Lachenbruch
        # have it besides their own phylogenetic gate): a lone permulation pass is a
        # different claim and gets its own label.
        evidence_count == 1L & sig_permulation ~ "Exploratory (phylogenetic)",
        evidence_count == 1L                   ~ "Exploratory (relative)",
        TRUE                 ~ "Not significant"
      )
    ) %>% dplyr::select(-lacks_null)
}

# ── RER null on the scale of the observed gene scores ─────────────────────────

# The observed RER scores are sign(Rho) * -log10(p.perm), p.perm being the empirical permulation p of the gene
# (RERconverge::permpvalcor(): the observed correlation against the gene's own null correlations, on its side of
# the null median, pseudo-count (num + 1) / (denom + 1)). The null of the pathway tests must be on that scale, so
# every null correlation gets the same p against the other null correlations of its gene (it is left out of its
# own row): p = #(null >= x) / #(null >= median) on the upper side, #(null <= x) / #(null <= median) on the
# lower one. The sign is the side of the median. Genes x N matrix like corStat, with the same dimnames.
fcs_empirical_corstat <- function(corRho) {
  Z   <- as.matrix(corRho)
  out <- matrix(0, nrow(Z), ncol(Z), dimnames = dimnames(Z))
  for (i in seq_len(nrow(Z))) {
    ok <- !is.na(Z[i, ])
    if (sum(ok) < 2L) next
    x    <- Z[i, ok]
    m    <- stats::median(x)
    up   <- x >= m
    ge   <- length(x) - rank(x, ties.method = "min") + 1   # null values >= x, itself included
    le   <- rank(x, ties.method = "max")                   # null values <= x, itself included
    p    <- ifelse(up, ge / sum(x >= m), le / sum(x <= m))
    out[i, ok] <- ifelse(up, 1, -1) * -log10(p)
  }
  out
}


# ── Full run ──────────────────────────────────────────────────────────────────

# Runs the three tests on every ranking and classifies the evidence.
# rankings: named list of named numeric vectors (zero-floored, see fcs_build_vals).
# gmts: named list of RERconverge gmt objects (fcs_load_gmts).
# perms_file: RDS of the permulation null, or "NO_FILE". Two shapes are read: RER
#   (corStat, genes x N, with optional corRho; one matrix shared by all rankings) and
#   CAAS (caas_corStat_byrank, one genes x N matrix per ranking global, top and bottom).
# Returns one row per (ranking, database, pathway): the columns of fcs_enrich_col_types.
fcs_run_all <- function(rankings, gmts, num_g = 10, max_g = 500, perms_file = "NO_FILE",
                        fdr_thr = 0.15, p_perm_thr = 0.025, n_perms_sum = 10000,
                        fdr_wilcoxon = fdr_thr, fdr_lachenbruch = fdr_thr,
                        fdr_permsum = fdr_thr, seed = 1998) {
  # Defined even without a perms file, because the Lachenbruch and path-sum loop below
  # runs regardless and reads corStat_byrk (all NULL here). It is replaced by the
  # per-ranking resolution once a perms file loads.
  corStat_byrank <- NULL; base_corStat <- NULL; base_corRho <- NULL
  corStat_byrk   <- setNames(vector("list", length(rankings)), names(rankings))
  res <- list(); alts <- list()
  for (rk in names(rankings)) {
    alts[[rk]] <- fcs_alternative(rankings[[rk]])
    r <- fcs_run_ranking(rankings[[rk]], gmts, num_g = num_g, max_g = max_g, alternative = alts[[rk]])
    if (nrow(r) > 0) res[[rk]] <- dplyr::mutate(r, ranking = rk)
  }
  if (length(res) == 0) {
    return(tibble::tibble(
      ranking = character(), database = character(), pathway = character(),
      stat = numeric(), pval = numeric(), p.adj = numeric(), p.perm = numeric(),
      num.genes = numeric(), gene.vals = character(),
      lach_pval = numeric(), lach_p.adj = numeric(), lach_chi_binary = numeric(),
      lach_chi_nonzero = numeric(), lach_chi_total = numeric(), lach_frac_magnitude = numeric(),
      lach_p.perm = numeric(),
      perm_pval = numeric(), perm_p.adj = numeric(), perm_nes = numeric(),
      sig_wilcoxon = logical(), sig_lachenbruch = logical(), sig_permulation = logical(),
      evidence_count = integer(), evidence_label = character()
    ))
  }
  enrich_df <- dplyr::bind_rows(res)
  enrich_df$p.perm <- NA_real_

  if (!is.null(perms_file) && perms_file != "NO_FILE" && file.exists(perms_file)) {
    fcs_progress(paste0("Loading null permulations from: ", perms_file))
    corperms <- tryCatch(readRDS(perms_file), error = function(e) NULL)

    # Two shapes of the perms RDS:
    #   * RER: corperms$corStat (genes x N) and optionally corRho, one matrix shared by
    #          the rankings; accelerating and decelerating are derived from corRho.
    #   * CAAS: corperms$caas_corStat_byrank, a list keyed by ranking name
    #          (global, top, bottom) of genes x N null matrices. The direction is already
    #          in each matrix (the side of each labeling partitions it), so no corRho
    #          transform is needed.
    corStat_byrank <- c(corperms[["corStat_byrank"]], corperms[["caas_corStat_byrank"]])
    # The corStat_byrank entries named *_asr (global_asr, top_asr, bottom_asr) do not match
    # the ranking names (global, top, bottom), so they are never selected here; the
    # ancestral-state-axis null is used by 11.Scoring_report.Rmd.

    if (!is.null(corperms[["caas_corStat_byrank"]])) {
      .stat_have <- corperms[["gene_stat"]]
      if (is.null(.stat_have)) {
        # An RDS without a gene_stat stamp is taken as built with size_adj_max, the
        # aggregator scoring_caas_perms.R uses; it is reported and kept.
        fcs_progress(paste0(
          "CAAS permulation null carries no gene_stat stamp (pre-stamp build); ",
          "assuming 'size_adj_max' (the only historical aggregator) and using it."))
      } else if (!identical(as.character(.stat_have), "size_adj_max")) {
        warning(sprintf(paste0(
          "CAAS permulation null is STALE: built with gene statistic '%s', but the ",
          "observed gene_caas_score uses 'size_adj_max'. p.perm left NA for the CAAS ",
          "rankings. Rebuild it (no ASR replay needed, minutes):\n",
          "  python3 subworkflows/CT_DISAMBIGUATION/local/reaggregate_perm_scores.py ",
          "--detail <run>/caas_permulation/perm_pos_detail --output-dir <run>/caas_permulation\n",
          "  Rscript subworkflows/SCORING/local/src/scoring_caas_perms.R ",
          "--gene-cycle-scores <run>/caas_permulation/gene_cycle_scores.tsv ",
          "--output <run>/caas_permulation/caas_perms.rds"),
          as.character(.stat_have)))
        corStat_byrank <- corStat_byrank[
          !names(corStat_byrank) %in% names(corperms[["caas_corStat_byrank"]])]
      }
    }
    if (!is.null(corperms) && (!is.null(corperms[["corStat"]]) || !is.null(corStat_byrank))) {
      base_corStat <- if (!is.null(corperms[["corStat"]])) as.matrix(corperms[["corStat"]]) else NULL
      base_corRho  <- if (!is.null(corperms[["corRho"]])) as.matrix(corperms[["corRho"]]) else NULL
      # RER: the null is the signed -log10 empirical p of each null correlation (the scale of the observed scores,
      # see fcs_empirical_corstat); its sign, the side of the gene's null median, gives the direction of the
      # accelerating and decelerating rankings. Without corRho the parametric corStat is kept.
      if (!is.null(base_corStat) && !is.null(base_corRho)) {
        base_corStat <- fcs_empirical_corstat(base_corRho)
        base_corRho  <- base_corStat
      }
      n_perms      <- if (!is.null(base_corStat)) ncol(base_corStat)
                      else if (length(corStat_byrank)) ncol(corStat_byrank[[1]]) else 0L
      fcs_progress(sprintf("Vectorized pathway permulations: N=%d, %d GMTs, %d rankings",
                           n_perms, length(gmts), length(rankings)))

      # Null matrix of every ranking, resolved once for this loop and for the
      # Lachenbruch and path-sum loop below.
      corStat_byrk <- setNames(
        lapply(names(rankings), fcs_resolve_corstat_rk,
               corStat_byrank = corStat_byrank, base_corStat = base_corStat, base_corRho = base_corRho),
        names(rankings))

      for (rk in names(rankings)) {
        obs_rk <- enrich_df %>% dplyr::filter(ranking == rk)
        if (nrow(obs_rk) == 0) next
        t0  <- Sys.time()
        alt <- alts[[rk]]

        # Null of this ranking (resolved above).
        corStat_rk <- corStat_byrk[[rk]]
        if (is.null(corStat_rk)) next  # no null for this ranking → p.perm stays NA

        # Observed statistic per set and database (rows = pathways).
        realenrich <- list()
        for (db in unique(obs_rk$database)) {
          db_df <- obs_rk %>% dplyr::filter(database == db) %>% as.data.frame()
          rownames(db_df) <- db_df$pathway
          realenrich[[db]] <- db_df[, c("pval", "stat"), drop = FALSE]
        }

        enrichStat <- tryCatch(
          fcs_null_enrichstat_vectorized(corStat_rk, gmts, realenrich, num_g = num_g, max_g = max_g),
          error = function(e) { fcs_progress(sprintf("null stats failed [%s]: %s", rk, e$message)); NULL })
        if (is.null(enrichStat)) next

        ppv <- fcs_permpvalenrich_vectorized(realenrich, enrichStat, alternative = alt)
        for (db in names(ppv)) {
          idx <- which(enrich_df$ranking == rk & enrich_df$database == db)
          if (length(idx)) enrich_df$p.perm[idx] <- ppv[[db]][enrich_df$pathway[idx]]
        }
        fcs_progress(sprintf("ranking %-13s done (%s) | %d sets across %d GMTs | %.1fs",
                             rk, alt, nrow(obs_rk), length(realenrich),
                             as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      }
    }
  }

  # ── Lachenbruch two-part and path-sum permulation (non-negative rankings) ──
  # Skipped for two-sided rankings (signed RER values), where score > 0 does not mean
  # "signal present"; those get the Wilcoxon test only.
  lach_res <- list()
  perm_res <- list()
  for (rk in names(rankings)) {
    if (alts[[rk]] != "greater") next
    vals_rk <- rankings[[rk]]
    # The null the Wilcoxon loop uses (corStat_byrk). NULL without a perms file or a null
    # for this ranking; the empirical p-values of both tests below then stay NA.
    corStat_rk <- corStat_byrk[[rk]]

    fcs_progress(sprintf("Lachenbruch two-part test: ranking %s", rk))
    lach_rk <- tryCatch(
      fcs_run_lachenbruch(vals_rk, gmts, num_g = num_g, max_g = max_g),
      error = function(e) {
        fcs_progress(sprintf("  Lachenbruch failed [%s]: %s", rk, e$message))
        tibble::tibble()
      })
    if (nrow(lach_rk) > 0) {
      lach_rk$lach_p.perm <- if (!is.null(corStat_rk)) {
        tryCatch(
          fcs_compute_lach_p_perm(lach_rk, corStat_rk, gmts, num_g = num_g, max_g = max_g),
          error = function(e) {
            fcs_progress(sprintf("  Lachenbruch null failed [%s]: %s", rk, e$message))
            NA_real_
          })
      } else NA_real_
      lach_res[[rk]] <- dplyr::mutate(lach_rk, ranking = rk)
    }

    if (is.null(corStat_rk)) {
      # A label shuffle ignores the phylogeny: it would fill the permulation columns with a
      # weaker null under the same names, so they stay NA.
      fcs_progress(sprintf("Path sum permulation: ranking %s skipped (no permulation null for this ranking)", rk))
      next
    }
    fcs_progress(sprintf("Path sum permulation: ranking %s (CAAS/RER null, %d perms)", rk, ncol(corStat_rk)))
    perm_rk <- tryCatch(
      fcs_run_permulation(vals_rk, gmts, num_g = num_g, max_g = max_g, n_perms = n_perms_sum,
                          seed = seed, null_mat = corStat_rk),
      error = function(e) {
        fcs_progress(sprintf("  Permulation failed [%s]: %s", rk, e$message))
        tibble::tibble()
      })
    if (nrow(perm_rk) > 0) perm_res[[rk]] <- dplyr::mutate(perm_rk, ranking = rk)
  }

  lach_df <- if (length(lach_res) > 0) dplyr::bind_rows(lach_res) else tibble::tibble()
  perm_df <- if (length(perm_res) > 0) dplyr::bind_rows(perm_res) else tibble::tibble()

  if (nrow(lach_df) > 0) {
    enrich_df <- dplyr::left_join(enrich_df, lach_df, by = c("ranking", "database", "pathway"))
  } else {
    enrich_df <- dplyr::mutate(enrich_df,
      lach_pval = NA_real_, lach_p.adj = NA_real_, lach_chi_binary = NA_real_,
      lach_chi_nonzero = NA_real_, lach_chi_total = NA_real_, lach_frac_magnitude = NA_real_,
      lach_p.perm = NA_real_)
  }
  if (nrow(perm_df) > 0) {
    enrich_df <- dplyr::left_join(enrich_df, perm_df, by = c("ranking", "database", "pathway"))
  } else {
    enrich_df <- dplyr::mutate(enrich_df,
      perm_pval = NA_real_, perm_p.adj = NA_real_, perm_nes = NA_real_)
  }

  # ── Evidence gates and classification ──
  enrich_df <- fcs_classify_evidence(enrich_df, fdr_wilcoxon = fdr_wilcoxon, fdr_lachenbruch = fdr_lachenbruch,
                                     fdr_permsum = fdr_permsum, p_perm_thr = p_perm_thr,
                                     no_null_rankings = names(corStat_byrk)[vapply(corStat_byrk, is.null, logical(1))])

  enrich_df %>% dplyr::relocate(ranking, database, pathway,
                                 evidence_count, evidence_label)
}

# ── Leading edge ──────────────────────────────────────────────────────────────

# One row per (ranking, database, pathway, gene) from the gene.vals ("gene:rank" list) of
# enrich_df, joined by gene to the annotation columns of attr_df when given.
fcs_annotate_leading_edge <- function(enrich_df, attr_df = NULL) {
  if (nrow(enrich_df) == 0) {
    # matrix(nrow = 0) makes every column logical, which would fail a later left_join
    # against a character `gene` column. The columns are typed explicitly, and those of
    # attr_df keep their types through a 0-row slice.
    base_df <- tibble::tibble(
      ranking = character(0), database = character(0), pathway = character(0),
      stat = numeric(0), p.adj = numeric(0),
      gene = character(0), gene_rank = numeric(0)
    )
    if (!is.null(attr_df)) {
      extra_cols <- setdiff(names(attr_df), names(base_df))
      base_df <- dplyr::bind_cols(base_df, attr_df[0, extra_cols, drop = FALSE])
    }
    return(base_df)
  }
  le <- enrich_df %>%
    dplyr::filter(!is.na(gene.vals) & nzchar(gene.vals)) %>%
    dplyr::select(dplyr::any_of(c("ranking", "database", "pathway", "stat", "p.adj", "gene.vals"))) %>%
    tidyr::separate_rows(gene.vals, sep = ",\\s*") %>%
    tidyr::separate(gene.vals, into = c("gene", "gene_rank"), sep = ":", fill = "right", extra = "merge") %>%
    dplyr::mutate(gene = trimws(gene), gene_rank = suppressWarnings(as.numeric(gene_rank)))
  if (!is.null(attr_df) && "gene" %in% names(attr_df)) {
    le <- le %>% dplyr::left_join(attr_df, by = "gene")
  }
  le
}

# Per-pathway summary of a leading-edge table: n_le (number of leading-edge genes) and,
# for each column of attr_cols present, the percentage of those genes with a true value.
fcs_leading_edge_summary <- function(le, attr_cols) {
  present <- intersect(attr_cols, names(le))
  if (nrow(le) == 0) {
    # Typed empty table, for the same reason as in fcs_annotate_leading_edge: database and
    # pathway must stay character for the later joins.
    base_df <- tibble::tibble(
      ranking = character(0), database = character(0), pathway = character(0),
      stat = numeric(0), p.adj = numeric(0), n_le = integer(0)
    )
    if (length(present) > 0) {
      base_df <- dplyr::bind_cols(base_df, le[0, present, drop = FALSE])
    }
    return(base_df)
  }
  
  res_df <- le %>%
    dplyr::group_by(ranking, database, pathway, stat, p.adj)
    
  if (length(present) > 0) {
    res_df <- res_df %>%
      dplyr::summarise(
        n_le = dplyr::n(),
        dplyr::across(dplyr::all_of(present), ~ round(100 * mean(.x, na.rm = TRUE), 1)),
        .groups = "drop"
      )
  } else {
    res_df <- res_df %>%
      dplyr::summarise(
        n_le = dplyr::n(),
        .groups = "drop"
      )
  }
  
  res_df %>% dplyr::arrange(p.adj, dplyr::desc(stat))
}
