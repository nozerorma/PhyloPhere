#!/usr/bin/env Rscript
# fcs_enrich.R — Functional class scoring (FCS) core: two gene-set tests per ranking.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: fcs_compute.R (FCS_COMPUTE_BATCHED process, fcs.nf) and 12.FCS_general_report.Rmd.
# Defines functions only.
#
# Tests the gene sets of GMT files against a gene ranking with two complementary tests:
#   1. Wilcoxon-AUC (RERconverge::fastwilcoxGMTall): rank shift of the set. On a zero-floored
#      ranking its one-sided p carries the tie term of the rank sum, and its permulation p
#      compares that p (not the AUC) with the same p of every null column, so that columns
#      with a different number of tied genes share one scale.
#   2. Lachenbruch two-part: prevalence of nonzero scorers (hypergeometric) plus magnitude among
#      them (rank-sum), combined as chi-square with 2 df. Safe for zero-inflated scores. The
#      observed ranking and every null column go through one code (fcs_lach_prepare,
#      fcs_lach_stats), which POSENRICH (posenrich_enrich.py, lachenbruch_columns) repeats.
# Each test casts one vote (FDR below its threshold, a direction check for the Wilcoxon, and
# the permulation p-value gate where a null exists). A term is "Supported" with 2 votes and
# "Exploratory" with 1. Lachenbruch runs only on non-negative (zero-floored) rankings; signed
# rankings (RER, two-sided) use the Wilcoxon test alone.
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

# ── Ties ──────────────────────────────────────────────────────────────────────

# Tie term of the variance of a rank-sum statistic: sum over the groups of equal values of t^3 - t. A zero-floored
# ranking has most of its genes tied at 0, and the variance of the statistic shrinks with that term.
fcs_tie_term <- function(x) {
  l <- rle(sort(x))$lengths
  sum(as.numeric(l)^3 - l)
}

# Standard deviation of the rank-sum U of n1 genes against n2 others, N = n1 + n2, corrected for ties. Without ties
# (tie = 0) it is sqrt(n1 * n2 * (N + 1) / 12), the one of RERconverge::simpleAUCgenesRanks.
fcs_rank_sum_sd <- function(n1, n2, tie) {
  N <- n1 + n2
  sqrt(n1 * n2 / 12 * ((N + 1) - tie / pmax(N * (N - 1), 1)))
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
    # GMT, as in fastwilcoxGMT. The variance of the rank sum carries the tie term of
    # that background: a zero-floored ranking is mostly ties, and the variance without
    # the term is several times too large.
    if (alternative == "greater") {
      gmt  <- gmts_to_run[[db]]
      bg   <- vals[intersect(unique(unlist(gmt$genesets)), names(vals))]
      n_db <- length(bg)
      n1   <- res$num.genes
      n2   <- n_db - n1
      U    <- (res$stat + 0.5) * n1 * n2
      mu   <- n1 * n2 / 2
      sdv  <- fcs_rank_sum_sd(n1, n2, fcs_tie_term(bg))
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
#   Part 1 (prevalence): one-sided Fisher exact test of the 2x2 table (score > 0 or = 0, in the set or not), the upper tail
#     of a hypergeometric distribution over the whole ranking.
#   Part 2 (magnitude): one-sided Wilcoxon rank-sum of the positive scores of the set against the positive scores of the
#     rest of the ranking, by the normal approximation with continuity correction and with the tie term of the positives.
#     It needs at least 2 positives in the set and 2 outside it.
# Each p-value (floored at 1e-15) becomes chi-square with 1 df; their sum is chi-square with 2 df, which gives the combined
# p-value. BH per GMT. Only for non-negative, zero-floored vals.
# The observed ranking and every permulation column go through the same code (fcs_lach_prepare, fcs_lach_stats), so a null
# value is the value the observed ranking would have if it were that column.

# Positives of every column of X (genes x columns; NA and 0 are no signal), their ranks among the positives of the column
# (average ties) and the tie term of those ranks.
fcs_lach_prepare <- function(X) {
  X <- as.matrix(X)
  X[is.na(X)] <- 0
  ii <- vector("list", ncol(X)); xx <- vector("list", ncol(X)); tie <- numeric(ncol(X))
  for (j in seq_len(ncol(X))) {
    p <- which(X[, j] > 0)
    ii[[j]] <- p
    if (length(p)) { xx[[j]] <- rank(X[p, j]); tie[j] <- fcs_tie_term(X[p, j]) } else xx[[j]] <- numeric(0)
  }
  jj <- rep(seq_len(ncol(X)), lengths(ii)); iv <- as.integer(unlist(ii)); rv <- as.numeric(unlist(xx))
  dims <- dim(X)
  list(pos  = Matrix::sparseMatrix(i = iv, j = jj, x = rep(1, length(iv)), dims = dims),
       rank = Matrix::sparseMatrix(i = iv, j = jj, x = rv, dims = dims),
       m = lengths(ii), tie = tie, n_genes = nrow(X))
}

# Sets of one GMT inside a universe: the sets with num_g to max_g genes of the universe, as a sparse membership matrix
# (sets x universe) and the size of each. NULL when no set passes.
fcs_lach_index <- function(gmt, universe, num_g = 10, max_g = 500) {
  gs <- gmt$genesets
  if (is.null(names(gs))) names(gs) <- gmt$geneset.names
  n1 <- vapply(gs, function(set) length(intersect(set, universe)), 1L)
  has_max <- !is.null(max_g) && is.finite(max_g) && max_g > 0
  keep <- n1 >= num_g & (!has_max | n1 <= max_g)
  if (!any(keep)) return(NULL)
  gs <- gs[keep]
  list(names = names(gs), n1 = as.numeric(n1[keep]), M = fcs_membership_matrix(gs, names(gs), universe))
}

# Lachenbruch chi-squares of the sets of one GMT in every column prepared by fcs_lach_prepare. Returns matrices sets x
# columns: chi1 (prevalence), chi2 (magnitude, 0 where it cannot run) and chi = chi1 + chi2.
fcs_lach_stats <- function(prep, M, n1) {
  nc <- length(prep$m)
  k1 <- as.matrix(M %*% prep$pos)
  mm <- matrix(prep$m, nrow(M), nc, byrow = TRUE)
  p1 <- phyper(as.vector(k1) - 1, as.vector(mm), prep$n_genes - as.vector(mm), rep(n1, nc), lower.tail = FALSE)
  chi1 <- matrix(qchisq(pmax(p1, 1e-15), df = 1, lower.tail = FALSE), nrow(M), nc)
  n2p <- mm - k1
  U2 <- as.matrix(M %*% prep$rank) - k1 * (k1 + 1) / 2
  z <- (U2 - k1 * n2p / 2 - 0.5) / fcs_rank_sum_sd(k1, n2p, matrix(prep$tie, nrow(M), nc, byrow = TRUE))
  chi2 <- matrix(qchisq(pmax(pnorm(z, lower.tail = FALSE), 1e-15), df = 1, lower.tail = FALSE), nrow(M), nc)
  chi2[!(k1 >= 2 & n2p >= 2) | is.na(chi2)] <- 0
  list(chi1 = chi1, chi2 = chi2, chi = chi1 + chi2)
}

# Observed side. Returns per pathway: lach_pval, lach_p.adj, lach_chi_binary (Part 1), lach_chi_nonzero (Part 2; 0 when
#   Part 2 cannot run), lach_chi_total and lach_frac_magnitude (lach_chi_nonzero / lach_chi_total; high when magnitude and
#   not prevalence drives it).
fcs_run_lachenbruch <- function(vals, gmts, num_g = 10, max_g = 500) {
  vals[is.na(vals)] <- 0
  prep <- fcs_lach_prepare(matrix(vals, ncol = 1, dimnames = list(names(vals), NULL)))
  out <- list()
  for (db in names(gmts)) {
    idx <- fcs_lach_index(gmts[[db]], names(vals), num_g, max_g)
    if (is.null(idx)) next
    st <- fcs_lach_stats(prep, idx$M, idx$n1)
    chi <- st$chi[, 1]
    db_df <- tibble::tibble(
      database            = db,
      pathway             = idx$names,
      lach_pval           = pchisq(chi, df = 2, lower.tail = FALSE),
      lach_chi_binary     = st$chi1[, 1],
      lach_chi_nonzero    = st$chi2[, 1],
      lach_chi_total      = chi,
      lach_frac_magnitude = ifelse(chi > 0, st$chi2[, 1] / chi, NA_real_)
    )
    db_df$lach_p.adj <- p.adjust(db_df$lach_pval, method = "BH")
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

# ── Wilcoxon: null of the tie-corrected p ─────────────────────────────────────

# The null of the one-sided p of the Wilcoxon-AUC test: for every permulation column, the p that fcs_run_ranking gives a set
# if that column were the ranking (normal approximation with the tie term of the column over the annotated genes of the GMT).
# A column of a zero-floored ranking has its own number of tied genes, so the AUC is not on one scale across the columns of
# the null and across the observed ranking; the standardized p is. Built from the null of the statistic (enrichStat, from
# fcs_null_enrichstat_vectorized) and the tie term of every column.
# Returns a named list db -> (sets x N) matrix, NA where the statistic is NA.
fcs_null_wilcoxon_p_vectorized <- function(corStat, gmts, realenrich, enrichStat, num_g = 10, max_g = 500) {
  out <- list()
  genes_all <- rownames(corStat)
  for (db in names(realenrich)) {
    stat <- enrichStat[[db]]
    gmt <- gmts[[db]]
    set_names <- rownames(realenrich[[db]])
    if (is.null(gmt) || is.null(stat) || length(set_names) == 0) { out[[db]] <- stat; next }
    gs <- gmt$genesets; names(gs) <- gmt$geneset.names
    if (!is.null(max_g) && is.finite(max_g) && max_g > 0) {
      gs <- gs[vapply(gs, function(set) length(intersect(set, genes_all)) <= max_g, logical(1))]
    }
    genes_db <- intersect(unique(unlist(gs)), genes_all)
    if (length(genes_db) < 3) { out[[db]] <- stat * NA_real_; next }
    M     <- fcs_membership_matrix(gs, set_names, genes_db)
    sub   <- corStat[genes_db, , drop = FALSE]
    notNA <- !is.na(sub)
    n1    <- as.matrix(M %*% (notNA * 1.0))
    n2    <- matrix(colSums(notNA), nrow(M), ncol(sub), byrow = TRUE) - n1
    tie   <- matrix(vapply(seq_len(ncol(sub)), function(j) fcs_tie_term(sub[notNA[, j], j]), 0), nrow(M), ncol(sub), byrow = TRUE)
    out[[db]] <- pnorm((stat + 0.5) * n1 * n2, n1 * n2 / 2, fcs_rank_sum_sd(n1, n2, tie), lower.tail = FALSE)
  }
  out
}

# Empirical permulation p-value of the one-sided Wilcoxon p, per set: (1 + #{null p <= observed p}) / (N_valid + 1); the
# smaller p is the more extreme.
fcs_permpval_from_p_vectorized <- function(realenrich, nullP) {
  out <- list()
  for (db in names(realenrich)) {
    null <- nullP[[db]]
    if (is.null(null) || nrow(null) == 0) next
    obs   <- realenrich[[db]]$pval[match(rownames(null), rownames(realenrich[[db]]))]
    count <- rowSums(null <= obs * (1 + 1e-9), na.rm = TRUE)
    p <- (count + 1) / (rowSums(!is.na(null)) + 1)
    p[is.na(obs)] <- NA_real_
    names(p) <- rownames(null)
    out[[db]] <- p
  }
  out
}

# The null of a ranking on the genes of that ranking: rows in the order of `genes`, and a gene without a row in the null has
# no signal (0) in every permulation. Observed and null then share the universe, so a set has the same genes on both sides.
fcs_align_null <- function(m, genes) {
  m <- as.matrix(m)
  out <- matrix(0, length(genes), ncol(m), dimnames = list(genes, colnames(m)))
  common <- intersect(genes, rownames(m))
  out[common, ] <- m[common, , drop = FALSE]
  out
}

# ── Lachenbruch empirical p-value ─────────────────────────────────────────────

# Empirical lach_p.perm of every row of lach_rk: the chi-square of the set in every column of the null, computed as for the
# observed ranking, against the observed lach_chi_total: (1 + #{null >= observed}) / (N + 1). The null has the genes of
# the ranking as rows (fcs_align_null). Returns a vector aligned to the rows of lach_rk, next to the analytic lach_pval and
# lach_p.adj.
fcs_compute_lach_p_perm <- function(lach_rk, corStat_rk, gmts, num_g = 10, max_g = 500) {
  prep <- fcs_lach_prepare(corStat_rk)
  out  <- rep(NA_real_, nrow(lach_rk))
  for (db in unique(lach_rk$database)) {
    idx <- fcs_lach_index(gmts[[db]], rownames(corStat_rk), num_g, max_g)
    if (is.null(idx)) next
    nullchi <- fcs_lach_stats(prep, idx$M, idx$n1)$chi
    rows <- which(lach_rk$database == db)
    k    <- match(lach_rk$pathway[rows], idx$names)
    cnt  <- vapply(seq_along(rows), function(i) if (is.na(k[i])) NA_real_ else sum(nullchi[k[i], ] >= lach_rk$lach_chi_total[rows[i]] - 1e-12), 0)
    out[rows] <- (cnt + 1) / (ncol(nullchi) + 1)
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
    sig_wilcoxon = readr::col_logical(), sig_lachenbruch = readr::col_logical(),
    evidence_count = readr::col_integer(), evidence_label = readr::col_character(),
    .default = readr::col_guess()
  )
}

# ── Evidence classification ───────────────────────────────────────────────────

# Evidence gates and labels. The FDR threshold of each test is set in conf/enrichment.config
# (fdr_wilcoxon, fdr_lachenbruch).
#   sig_wilcoxon: FDR gate, direction (stat > 0) and, when p.perm exists, the permulation
#     gate p.perm < p_perm_thr. p.perm is NA without a perms file or null, and the gate is skipped.
#   sig_lachenbruch: FDR gate and, when lach_p.perm exists, the same permulation gate on the shared null.
# A ranking without a permulation null (no perms file, an empty null, a stale one) has no
# phylogenetic gate: a row that passes the gates that exist is "Exploratory", never "Supported".
# no_null_rankings lists those rankings.
fcs_classify_evidence <- function(enrich_df, fdr_wilcoxon, fdr_lachenbruch, p_perm_thr, no_null_rankings = character(0)) {
  enrich_df %>%
    dplyr::mutate(
      sig_wilcoxon    = !is.na(p.adj)       & p.adj       < fdr_wilcoxon &
                        (is.na(p.perm) | p.perm < p_perm_thr) &
                        !is.na(stat)  & stat > 0,
      sig_lachenbruch = !is.na(lach_p.adj)  & lach_p.adj  < fdr_lachenbruch &
                        (is.na(lach_p.perm) | lach_p.perm < p_perm_thr),
      lacks_null      = ranking %in% no_null_rankings,
      evidence_count  = as.integer(sig_wilcoxon) + as.integer(sig_lachenbruch),
      evidence_label  = dplyr::case_when(
        # no null for this ranking: no phylogenetic gate was applied, whatever else passed
        lacks_null & evidence_count >= 1L ~ "Exploratory",
        lacks_null                        ~ "Not significant",
        evidence_count == 2L ~ "Supported",
        evidence_count == 1L ~ "Exploratory",
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

# Runs the tests on every ranking and classifies the evidence.
# rankings: named list of named numeric vectors (zero-floored, see fcs_build_vals).
# gmts: named list of RERconverge gmt objects (fcs_load_gmts).
# perms_file: RDS of the permulation null, or "NO_FILE". Two shapes are read: RER
#   (corStat, genes x N, with optional corRho; one matrix shared by all rankings) and
#   CAAS (caas_corStat_byrank, one genes x N matrix per ranking global, top and bottom).
# Returns one row per (ranking, database, pathway): the columns of fcs_enrich_col_types.
fcs_run_all <- function(rankings, gmts, num_g = 10, max_g = 500, perms_file = "NO_FILE",
                        fdr_thr = 0.15, p_perm_thr = 0.025,
                        fdr_wilcoxon = fdr_thr, fdr_lachenbruch = fdr_thr) {
  # Defined even without a perms file, because the Lachenbruch loop below
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
      sig_wilcoxon = logical(), sig_lachenbruch = logical(),
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
      # Lachenbruch loop below.
      corStat_byrk <- setNames(
        lapply(names(rankings), fcs_resolve_corstat_rk,
               corStat_byrank = corStat_byrank, base_corStat = base_corStat, base_corRho = base_corRho),
        names(rankings))
      # The null of a magnitude ranking on the genes of that ranking: observed and null share the universe, so a set has
      # the same genes on both sides and the ties of the null are counted over the same genes.
      for (rk in names(corStat_byrk)) {
        if (!is.null(corStat_byrk[[rk]]) && alts[[rk]] == "greater") {
          corStat_byrk[[rk]] <- fcs_align_null(corStat_byrk[[rk]], names(rankings[[rk]]))
        }
      }

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

        # Magnitude rankings: the permulated statistic is the tie-corrected p of the AUC, on one scale across the columns
        # of the null whatever their number of tied genes. Signed rankings have no ties: the AUC itself, as RERconverge.
        if (alt == "greater") {
          nullP <- fcs_null_wilcoxon_p_vectorized(corStat_rk, gmts, realenrich, enrichStat, num_g = num_g, max_g = max_g)
          ppv   <- fcs_permpval_from_p_vectorized(realenrich, nullP)
        } else {
          ppv <- fcs_permpvalenrich_vectorized(realenrich, enrichStat, alternative = alt)
        }
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

  # ── Lachenbruch two-part (non-negative rankings) ──
  # Skipped for two-sided rankings (signed RER values), where score > 0 does not mean
  # "signal present"; those get the Wilcoxon test only.
  lach_res <- list()
  for (rk in names(rankings)) {
    if (alts[[rk]] != "greater") next
    vals_rk <- rankings[[rk]]
    # The null the Wilcoxon loop uses (corStat_byrk). NULL without a perms file or a null
    # for this ranking; the empirical p-value below then stays NA.
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
  }

  lach_df <- if (length(lach_res) > 0) dplyr::bind_rows(lach_res) else tibble::tibble()

  if (nrow(lach_df) > 0) {
    enrich_df <- dplyr::left_join(enrich_df, lach_df, by = c("ranking", "database", "pathway"))
  } else {
    enrich_df <- dplyr::mutate(enrich_df,
      lach_pval = NA_real_, lach_p.adj = NA_real_, lach_chi_binary = NA_real_,
      lach_chi_nonzero = NA_real_, lach_chi_total = NA_real_, lach_frac_magnitude = NA_real_,
      lach_p.perm = NA_real_)
  }

  # ── Evidence gates and classification ──
  enrich_df <- fcs_classify_evidence(enrich_df, fdr_wilcoxon = fdr_wilcoxon, fdr_lachenbruch = fdr_lachenbruch,
                                     p_perm_thr = p_perm_thr,
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
