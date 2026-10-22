#!/usr/bin/env Rscript
# fcs_enrich.R — Functional class scoring (FCS) core: Wilcoxon-AUC test with permulation p.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: fcs_compute.R (FCS_COMPUTE_BATCHED process, fcs.nf) and 12.FCS_general_report.Rmd.
# Defines functions only.
#
# Tests the gene sets of GMT files against a gene ranking with the Wilcoxon-AUC test
# (RERconverge::fastwilcoxGMTall): rank shift of the set. On a zero-floored ranking its
# one-sided p carries the tie term of the rank sum. Permulation p-values compare the
# enrichment statistic (AUC - 0.5) with the null statistics, matching RERconverge.
#
# Significance (fcs_classify_significance) rests on the permulation p of each set: BH over
# p.perm within each ranking and database gives q.perm, a set is significant when
# q.perm < fdr_thr. A ranking without a permulation null has no p.perm and so no significant sets.
#
# The rule of Redlich et al. 2024 is reported next to it (redlich_level): the analytic p.adj (BH of the
# test within each database) below a level, and p.perm below p_perm_thr. redlich_level is the smallest of
# the levels (default 0.05, 0.10, 0.15) that the set meets, NA when it meets none. It does not enter sig_perm.
#
# Design:
#   * CAAS uses the full tested background (cleaned_background), with absent signal
#     floored to 0. RER uses only genes with a non-missing association statistic.
#   * Direction by membership: a "top" ranking keeps the score of the top and both
#     genes and floors the rest to 0; "bottom" keeps bottom and both. A gene of both
#     directions contributes its full score to each.
#   * Significance is not used to filter the input. Callers attach it to the leading-edge
#     genes as annotation (fcs_annotate_leading_edge).
#   * Multiple testing is BH within each GMT (each database is its own family of
#     hypotheses).
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
# CAAS returns a named vector over union(universe, names(scores)); genes without a
# score, or with NA, get `floor` (0). RER excludes missing scores and adds no genes
# from the background, matching getStat. Measured zero scores remain in both cases.
# The caller handles CAAS direction by passing the directional subset's scores.
# Gene IDs are opaque: keep suffixes; never silently collapse distinct IDs or
# select one of several scores for an exact duplicate. Missing scores do not
# exempt an input from this contract.
fcs_check_gene_ids <- function(ids, n, context) {
  if (n == 0L) return(invisible(NULL))
  if (length(ids) != n || anyNA(ids) || any(!nzchar(trimws(ids))))
    stop(context, ": missing or empty gene IDs")
  if (anyDuplicated(ids))
    stop(context, ": duplicate gene IDs: ", paste(head(unique(ids[duplicated(ids)]), 5), collapse = ", "))
  invisible(NULL)
}

# getStat treats a shared prefix before a dot as one gene. Reject that ambiguity
# rather than reproducing its inconsistent choice of which score to retain.
fcs_check_rer_gene_ids <- function(ids, n, context) {
  fcs_check_gene_ids(ids, n, context)
  base_ids <- sub("\\..*", "", ids)
  if (anyDuplicated(base_ids))
    stop(context, ": ambiguous RER gene IDs sharing a prefix before '.': ",
         paste(head(unique(base_ids[duplicated(base_ids)]), 5), collapse = ", "))
  invisible(NULL)
}

fcs_build_vals <- function(scores, universe, floor = 0, score_source = c("caas", "rer")) {
  fcs_check_gene_ids(names(scores), length(scores), "gene scores")
  score_source <- match.arg(score_source)
  if (score_source == "rer")
    fcs_check_rer_gene_ids(names(scores), length(scores), "RER gene scores")
  scores <- scores[!is.na(scores)]
  if (score_source == "rer") return(scores)
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
# vals: named numeric vector (CAAS background with zero floor, or scored RER genes).
# gmts: named list of RERconverge gmt objects.
# num_g, max_g: minimum and maximum genes per set (max_g 0 or Inf = no limit).
# alternative: "greater" or "two.sided" (see fcs_alternative).
# Returns a tibble: database, pathway, stat (AUC - 0.5), pval, p.adj (BH per GMT),
#   num.genes, gene.vals (leading-edge members as "gene:rank").
# fastwilcoxGMT is internal to RERconverge; the exported fastwilcoxGMTall loops it over a
# named list of GMTs and adjusts per GMT, so it is called once and the databases stitched.
fcs_run_ranking <- function(vals, gmts, num_g = 10, max_g = 0, alternative = "two.sided") {
  fcs_check_gene_ids(names(vals), length(vals), "observed ranking")
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
    # Use the same evaluable domain as the vectorized null: at least three
    # annotated background genes outside the set. Native fastwilcoxGMT falls
    # back outside the GMT when its complement is empty, using incompatible
    # ranks; that degenerate result is not an annotated-background AUC.
    gmt <- gmts_to_run[[db]]
    bg <- vals[intersect(unique(unlist(gmt$genesets)), names(vals))]
    res <- res[length(bg) - res$num.genes > 2, , drop = FALSE]
    if (nrow(res) == 0) next
    res$p.adj <- p.adjust(res$pval, method = "BH")
    # The p-value of fastwilcoxGMT (simpleAUCgenesRanks) is always two-sided. For
    # one-sided magnitude rankings (non-negative, zero-floored: CAAS, FADE and
    # accumulation scores) the analytic p-value is
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

# RER binary/continuous rankings are two-sided, even if this sample has no negative
# scores. Other rankings infer direction from their values unless the caller sets RER.
fcs_alternative <- function(vals, score_source = c("auto", "caas", "rer")) {
  score_source <- match.arg(score_source)
  if (score_source == "rer" || any(vals < 0, na.rm = TRUE)) "two.sided" else "greater"
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
fcs_null_enrichstat_vectorized <- function(corStat, gmts, realenrich, num_g = 10, max_g = 0) {
  fcs_check_gene_ids(rownames(corStat), nrow(corStat), "null matrix")
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

# The null of a ranking on the genes of that ranking: rows in the order of `genes`, and a gene without a row in the null has
# no signal (0) in every permulation. Observed and null then share the universe, so a set has the same genes on both sides.
fcs_align_null <- function(m, genes) {
  m <- as.matrix(m)
  fcs_check_gene_ids(rownames(m), nrow(m), "null matrix")
  fcs_check_gene_ids(genes, length(genes), "observed universe")
  out <- matrix(0, length(genes), ncol(m), dimnames = list(genes, colnames(m)))
  common <- intersect(genes, rownames(m))
  out[common, ] <- m[common, , drop = FALSE]
  out
}

# ── Observed statistic ───────────────────────────────────────────────────────

# The observed enrichment of every ranking: one row per (ranking, database, pathway) from fcs_run_ranking, the test
# alternative of each ranking (`alts`), and `enrich` = NULL when no ranking has a set. Needs no null.
fcs_observed <- function(rankings, gmts, num_g = 10, max_g = 0, score_source = "auto") {
  res <- list(); alts <- list()
  for (rk in names(rankings)) {
    alts[[rk]] <- fcs_alternative(rankings[[rk]], score_source = score_source)
    r <- fcs_run_ranking(rankings[[rk]], gmts, num_g = num_g, max_g = max_g, alternative = alts[[rk]])
    if (nrow(r) > 0) res[[rk]] <- dplyr::mutate(r, ranking = rk)
  }
  list(enrich = if (length(res) > 0) dplyr::bind_rows(res) else NULL, alts = alts)
}

# ── Permulation counts by block of null columns ──────────────────────────────

# Counts of the permulation p of every set of one ranking over the columns of `corStat` (genes x columns, any subset of
# the null): magnitude rankings compare AUC - 0.5 with `>=`; signed rankings compare `abs(AUC - 0.5)` with `>=`, matching
# RERconverge's enrichment-statistic comparison. `valid` counts evaluable null columns. Both are sums over columns, so
# counts of disjoint blocks add up to the counts of the whole null.
fcs_null_count <- function(corStat, gmts, realenrich, alt, num_g = 10, max_g = 0) {
  out <- list()
  null_all <- fcs_null_enrichstat_vectorized(corStat, gmts, realenrich, num_g = num_g, max_g = max_g)
  for (db in names(realenrich)) {
    null <- null_all[[db]]
    if (is.null(null) || nrow(null) == 0) next
    obs <- realenrich[[db]]$stat[match(rownames(null), rownames(realenrich[[db]]))]
    count <- if (alt == "greater") rowSums(null >= obs, na.rm = TRUE)
             else                  rowSums(abs(null) >= abs(obs), na.rm = TRUE)
    count[is.na(obs)] <- NA_real_
    out[[db]] <- tibble::tibble(database = db, pathway = rownames(null),
                                count = as.numeric(count), valid = as.numeric(rowSums(!is.na(null))))
  }
  dplyr::bind_rows(out)
}

# Counts of the enrich rows `obs_rk` over a list of null matrices (genes x columns, one per block of cycles). The counts
# of the matrices are summed; a failed matrix gives NULL for the whole ranking, as the unblocked code did.
fcs_null_counts_by_matrices <- function(rk, obs_rk, mats, gmts, alt, num_g = 10, max_g = 0) {
  realenrich <- list()
  for (db in unique(obs_rk$database)) {
    db_df <- obs_rk %>% dplyr::filter(database == db) %>% as.data.frame()
    rownames(db_df) <- db_df$pathway
    realenrich[[db]] <- db_df[, c("pval", "stat"), drop = FALSE]
  }
  for (m in mats) fcs_check_gene_ids(rownames(m), nrow(m), "null matrix")
  parts <- lapply(mats, function(m) {
    tryCatch(fcs_null_count(m, gmts, realenrich, alt, num_g = num_g, max_g = max_g),
             error = function(e) { fcs_progress(sprintf("null stats failed [%s]: %s", rk, e$message)); NULL })
  })
  if (any(vapply(parts, is.null, logical(1)))) return(NULL)
  counts <- dplyr::bind_rows(parts) %>%
    dplyr::group_by(database, pathway) %>%
    dplyr::summarise(count = sum(count), valid = sum(valid), .groups = "drop")
  dplyr::mutate(counts, ranking = rk)
}

# Counts over the null in blocks of `block_cols` columns of one matrix (NULL: one block).
fcs_null_counts_by_block <- function(rk, obs_rk, corStat, gmts, alt, num_g = 10, max_g = 0, block_cols = NULL) {
  n <- ncol(corStat)
  if (is.null(block_cols) || block_cols >= n) block_cols <- n
  mats <- lapply(seq(1L, n, by = block_cols), function(s) corStat[, s:min(n, s + block_cols - 1L), drop = FALSE])
  fcs_null_counts_by_matrices(rk, obs_rk, mats, gmts, alt, num_g = num_g, max_g = max_g)
}

# Sets the permulation p of the enrich rows from the counts: p.perm = (count + 1) / (valid + 1), NA where there is no count.
fcs_close_counts <- function(enrich_df, counts) {
  if (is.null(counts) || nrow(counts) == 0) return(enrich_df)
  enrich_df$.row <- seq_len(nrow(enrich_df))
  m <- dplyr::left_join(enrich_df[, c(".row", "ranking", "database", "pathway")], counts,
                        by = c("ranking", "database", "pathway"))
  p <- (m$count + 1) / (m$valid + 1)
  hit <- !is.na(p)
  enrich_df$p.perm[m$.row[hit]] <- p[hit]
  enrich_df$.row <- NULL
  enrich_df
}

# ── Null per ranking ──────────────────────────────────────────────────────────

# Null matrix (genes x N) of ranking `rk`, or NULL when there is none.
# CAAS supplies one matrix per direction (corStat_byrank). RER supplies one matrix (base_corStat,
# the corStat of the permulations: sign(Rho) * -log10(P) of each null phenotype), used for its
# only ranking, "global".
fcs_resolve_corstat_rk <- function(rk, corStat_byrank, base_corStat) {
  if (!is.null(corStat_byrank) && !is.null(corStat_byrank[[rk]])) {
    return(as.matrix(corStat_byrank[[rk]]))
  }
  if (!is.null(base_corStat) && rk == "global") return(as.matrix(base_corStat))
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
    q.perm = readr::col_double(),
    num.genes = readr::col_double(), gene.vals = readr::col_character(),
    sig_perm = readr::col_logical(),
    direction = readr::col_character(), redlich_level = readr::col_double(),
    .default = readr::col_guess()
  )
}

# ── Significance classification ───────────────────────────────────────────────

# BH over p.perm within each (ranking, database) gives q.perm. A set is significant (sig_perm) when q.perm < fdr_thr.
# direction is the sign of stat for a
# significant set (a signed ranking is tested two-sided, so a set can be significant on either side). Rows without p.perm
# (no null for the ranking, or a set the null could not evaluate) have q.perm NA and are not significant.
# redlich_level: the smallest of analytic_levels with p.adj below it, among the sets with p.perm below p_perm_thr;
# NA otherwise. Neither classification adds a sign filter.
fcs_classify_significance <- function(enrich_df, fdr_thr, p_perm_thr, analytic_levels = c(0.05, 0.10, 0.15)) {
  enrich_df <- dplyr::select(enrich_df, -dplyr::any_of("sig_top"))
  levels <- sort(unique(analytic_levels), decreasing = TRUE)
  out <- enrich_df %>%
    dplyr::group_by(ranking, database) %>%
    dplyr::mutate(q.perm = p.adjust(p.perm, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sig_perm  = !is.na(q.perm) & q.perm < fdr_thr,
      direction = dplyr::case_when(!sig_perm ~ NA_character_, stat > 0 ~ "positive", stat < 0 ~ "negative", TRUE ~ NA_character_)
    )
  eligible <- !is.na(out$p.adj) & !is.na(out$p.perm) & out$p.perm < p_perm_thr
  out$redlich_level <- NA_real_
  for (L in levels) out$redlich_level[eligible & out$p.adj < L] <- L
  out
}


# ── Full run ──────────────────────────────────────────────────────────────────

# Runs the Wilcoxon-AUC test on every ranking and classifies the significance.
# rankings: named list of named numeric vectors (see fcs_build_vals).
# gmts: named list of RERconverge gmt objects (fcs_load_gmts).
# perms_file: RDS of the permulation null, or "NO_FILE". Two shapes are read: RER
#   (corStat, genes x N, the sign(Rho) * -log10(P) of every null phenotype; one matrix for its
#   "global" ranking) and CAAS (caas_corStat_byrank, one genes x N matrix per ranking global,
#   top and bottom).
# gene_stat: the gene statistic of the observed gene scores (size_adj_max, fmax or max); the CAAS null is used only when
#   it was built with the same one.
# score_source: "rer" fixes a two-sided test; "auto"/"caas" infer it from the scores.
# Returns one row per (ranking, database, pathway): the columns of fcs_enrich_col_types.
fcs_run_all <- function(rankings, gmts, num_g = 10, max_g = 0, perms_file = "NO_FILE",
                        fdr_thr = 0.15, p_perm_thr = 0.025, analytic_levels = c(0.05, 0.10, 0.15),
                        gene_stat = "size_adj_max", null_block_cols = NULL, score_source = "auto") {
  corStat_byrank <- NULL; base_corStat <- NULL
  observed <- fcs_observed(rankings, gmts, num_g = num_g, max_g = max_g, score_source = score_source)
  alts <- observed$alts
  if (is.null(observed$enrich)) {
    return(tibble::tibble(
      ranking = character(), database = character(), pathway = character(),
      stat = numeric(), pval = numeric(), p.adj = numeric(), p.perm = numeric(), q.perm = numeric(),
      num.genes = numeric(), gene.vals = character(),
      sig_perm = logical(), direction = character(), redlich_level = numeric()
    ))
  }
  enrich_df <- observed$enrich
  enrich_df$p.perm <- NA_real_

  if (!is.null(perms_file) && perms_file != "NO_FILE" && file.exists(perms_file)) {
    fcs_progress(paste0("Loading null permulations from: ", perms_file))
    corperms <- tryCatch(readRDS(perms_file), error = function(e) NULL)

    # Two shapes of the perms RDS:
    #   * RER: corperms$corStat (genes x N), one matrix for the "global" ranking.
    #   * CAAS: corperms$caas_corStat_byrank, a list keyed by ranking name
    #          (global, top, bottom) of genes x N null matrices. The direction is already
    #          in each matrix (the side of each labeling partitions it).
    corStat_byrank <- c(corperms[["corStat_byrank"]], corperms[["caas_corStat_byrank"]])
    # The RDS holds the CAAS axis only (caas_corStat_byrank: one matrix per ranking global, top and bottom).

    if (!is.null(corperms[["caas_corStat_byrank"]])) {
      .stat_have <- corperms[["gene_stat"]]
      if (is.null(.stat_have)) {
        # An RDS without a gene_stat stamp is taken as built with size_adj_max, the
        # aggregator of the default; it is reported and kept when that is the observed one.
        fcs_progress(paste0(
          "CAAS permulation null carries no gene_stat stamp (pre-stamp build); ",
          "assuming 'size_adj_max' (the only historical aggregator) and using it."))
        .stat_have <- "size_adj_max"
      }
      if (!identical(as.character(.stat_have), as.character(gene_stat))) {
        warning(sprintf(paste0(
          "CAAS permulation null is STALE: built with gene statistic '%s', but the ",
          "observed gene_caas_score uses '%s'. p.perm left NA for the CAAS ",
          "rankings. Rerun the null (CAAS_CORE_MERGE) to rebuild it."),
          as.character(.stat_have), as.character(gene_stat)))
        corStat_byrank <- corStat_byrank[
          !names(corStat_byrank) %in% names(corperms[["caas_corStat_byrank"]])]
      }
    }
    if (!is.null(corperms) && (!is.null(corperms[["corStat"]]) || !is.null(corStat_byrank))) {
      base_corStat <- if (!is.null(corperms[["corStat"]])) as.matrix(corperms[["corStat"]]) else NULL
      if (!is.null(base_corStat))
        fcs_check_rer_gene_ids(rownames(base_corStat), nrow(base_corStat), "RER null matrix")
      n_perms      <- if (!is.null(base_corStat)) ncol(base_corStat)
                      else if (length(corStat_byrank)) ncol(corStat_byrank[[1]]) else 0L
      fcs_progress(sprintf("Vectorized pathway permulations: N=%d, %d GMTs, %d rankings",
                           n_perms, length(gmts), length(rankings)))

      # Null matrix of every ranking.
      corStat_byrk <- setNames(
        lapply(names(rankings), fcs_resolve_corstat_rk,
               corStat_byrank = corStat_byrank, base_corStat = base_corStat),
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
        corStat_rk <- corStat_byrk[[rk]]
        if (is.null(corStat_rk)) next  # no null for this ranking: p.perm stays NA
        t0 <- Sys.time()
        counts <- fcs_null_counts_by_block(rk, obs_rk, corStat_rk, gmts, alts[[rk]], num_g = num_g, max_g = max_g,
                                           block_cols = null_block_cols)
        enrich_df <- fcs_close_counts(enrich_df, counts)
        fcs_progress(sprintf("ranking %-13s done (%s) | %d sets across %d GMTs | %.1fs",
                             rk, alts[[rk]], nrow(obs_rk), length(unique(obs_rk$database)),
                             as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      }
    }
  }

  # ── Significance ──
  enrich_df <- fcs_classify_significance(enrich_df, fdr_thr = fdr_thr, p_perm_thr = p_perm_thr, analytic_levels = analytic_levels)

  enrich_df %>% dplyr::relocate(ranking, database, pathway, stat, pval, p.adj, p.perm, q.perm,
                                num.genes, gene.vals, sig_perm, direction, redlich_level)
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
