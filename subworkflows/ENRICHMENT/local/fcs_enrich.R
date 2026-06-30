#!/usr/bin/env Rscript
# =============================================================================
# fcs_enrich.R — shared Functional Class Scoring (FCS) enrichment core
# =============================================================================
# Rank-based, threshold-free gene-set enrichment via the Wilcoxon rank-sum /
# Mann-Whitney AUC test (RERconverge::fastwilcoxGMT). This is the single
# enrichment engine for the ranked question across CAAS / RER / FADE /
# accumulation, replacing GSEA. STRING handles the set/network question.
#
# Design (locked in design discussion):
#   * Universe = the full tested background (cleaned_background). Genes with no
#     signal are floored to 0 — never dropped, never penalised (zero-floor),
#     mirroring RERconverge's accel/decel rankings.
#   * Direction via membership: a "top" ranking keeps the score for {top, both}
#     genes and floors everyone else to 0; "bottom" keeps {bottom, both}. A
#     both-direction gene contributes its full score to BOTH — never penalised.
#   * Significance is NEVER gated into the input. It rides along as annotation on
#     the leading-edge genes (gate_sig / gate_fdr / top-1/5/10% / cross-module).
#   * Multiple testing: BH per GMT (each database is its own hypothesis family);
#     fastwilcoxGMT already BH-adjusts within a single GMT call.
#
# Statistic returned by fastwilcoxGMT: stat = AUC - 0.5 in [-0.5, 0.5].
#   stat > 0 → pathway genes rank systematically HIGH on the score.
#
# This file defines functions only; it is sourced by FCS_general.Rmd.
# =============================================================================

suppressPackageStartupMessages({
  library(RERconverge)   # fastwilcoxGMT, read.gmt
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(parallel)
  library(Matrix)        # sparse membership matrix + %*% for the vectorized null
})

# ── GMT loading ──────────────────────────────────────────────────────────────
# Clean, version-agnostic database label from a GMT file name:
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

# Build an ID -> human-readable description map from the GMT files themselves.
# GMT line format: <setID>\t<description>\t<gene1>\t<gene2>...
# WebGestalt GMTs carry the term name in column 2; MSigDB carries a URL there, so
# for URL/empty descriptions we fall back to the (already readable) set ID.
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
    # TRANSFAC-style TF motif IDs (e.g. V$HNF4_01, TGTTTGY_V$HNF3_Q6) carry no
    # separate description — surface the transcription factor symbol so the
    # Description column is informative. Only applied when desc just repeats id.
    tf <- ifelse(grepl("V\\$", id), sub(".*V\\$([A-Za-z0-9]+).*", "\\1", id), NA_character_)
    desc <- ifelse(!is.na(tf) & nzchar(tf) & desc == id,
                   paste0(tf, " — TF target motif"), desc)
    rows[[db]] <- data.frame(database = db, pathway = id, description = desc,
                             stringsAsFactors = FALSE)
  }
  dplyr::distinct(dplyr::bind_rows(rows))
}

# Load every *.gmt in gmt_dir as a named list of RERconverge gmt objects
# (each has $genesets and $geneset.names — the structure fastwilcoxGMT expects).
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

# ── Ranking construction (zero-floor over the universe) ──────────────────────
# scores: named numeric vector (gene -> score) for the SIGNAL genes only.
# universe: character vector of all tested genes (cleaned_background).
# Returns a named numeric vector over union(universe, names(scores)) where genes
# without a score (or NA) are set to floor (0). Direction is handled by the
# caller passing only the directional subset's scores (top/both or bottom/both).
fcs_build_vals <- function(scores, universe, floor = 0) {
  scores <- scores[!is.na(scores)]
  genes  <- union(universe, names(scores))
  vals   <- setNames(rep(floor, length(genes)), genes)
  vals[names(scores)] <- as.numeric(scores)
  vals
}

# ── Run FCS for one ranking over all GMTs ────────────────────────────────────
# vals: named numeric vector (full universe, zero-floored).
# gmts: named list of RERconverge gmt objects.
# num_g: minimum genes per set (fastwilcoxGMT num.g).
# Returns a tibble: database, pathway, stat (AUC-0.5), pval, p.adj (BH per GMT),
#                   num.genes, gene.vals (leading-edge members as "gene:rank").
# NOTE: fastwilcoxGMT is internal to RERconverge; fastwilcoxGMTall is the exported
# wrapper (loops fastwilcoxGMT over a named annotation list and BH-adjusts per
# GMT). We call it once over all GMTs and stitch the per-database results.
fcs_run_ranking <- function(vals, gmts, num_g = 10, alternative = "two.sided") {
  reslist <- tryCatch(
    RERconverge::fastwilcoxGMTall(vals, gmts, outputGeneVals = TRUE, num.g = num_g),
    error = function(e) { message(sprintf("  fastwilcoxGMTall failed: %s", e$message)); NULL }
  )
  if (is.null(reslist) || length(reslist) == 0) return(tibble::tibble())
  vals <- vals[!is.na(vals)]
  out <- list()
  for (db in names(reslist)) {
    res <- reslist[[db]]
    if (is.null(res) || nrow(res) == 0) next
    # fastwilcoxGMT's simple-mode p-value (simpleAUCgenesRanks) is ALWAYS two-sided.
    # For one-sided ("greater") magnitude rankings — CAAS/FADE/accum scores and RER
    # accel/decel (all non-negative, zero-floored) — recompute the analytic p one-
    # sided from the AUC so a DEPLETED set (stat < 0) is not falsely flagged, then
    # re-BH per GMT. Matches RERconverge's alternative="greater" omnibus convention.
    # (Background = the GMT's own annotated genes, exactly as fastwilcoxGMT.)
    if (alternative == "greater") {
      gmt  <- gmts[[db]]
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

# ── Per-ranking test sidedness ────────────────────────────────────────────────
# A signed ranking (negatives present, e.g. RER global getStat = sign(Rho)*-log10P)
# is two-sided; a non-negative magnitude ranking (CAAS/FADE/accum scores, RER
# accel/decel zero-floored) is one-sided "greater". Detected from the values so no
# per-module plumbing is needed.
fcs_alternative <- function(vals) if (any(vals < 0, na.rm = TRUE)) "two.sided" else "greater"

# ── Progress logging (flushed; survives knitr chunk buffering) ────────────────
# Writes a timestamped line to stderr AND appends to fcs_progress.log in the work
# dir, so a run can be followed with `tail -f` even while an Rmd chunk is mid-flight
# (knitr buffers chunk stdout/stderr, so the file is the reliable channel).
fcs_progress <- function(msg, file = "fcs_progress.log") {
  line <- sprintf("[FCS %s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  message(line)
  try(cat(line, "\n", sep = "", file = file, append = TRUE), silent = TRUE)
}

# Per-column ranks (average ties); NA -> 0 so absent genes drop out of the rank-sum
# (mirrors fastwilcoxGMT's `vals = vals[!is.na(vals)]` then `rank()`).
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

# Sparse set x gene membership matrix over a fixed gene space.
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

# ── Vectorized permulation null statistics ────────────────────────────────────
# Reproduces fastwilcoxGMTall(corStat[,j], gmts)$stat for EVERY permulation column
# j in ONE sparse matrix multiply per GMT, replacing N x fastwilcoxGMTall calls.
# Faithful to fastwilcoxGMT: background = the GMT's own annotated genes; ranks are
# per-GMT, per-column, average-tie; sets failing num.g (or bkgenes<=2) in a column
# become NA for that column. Output: named list db -> (sets x N) AUC-0.5 matrix,
# rows aligned to rownames(realenrich[[db]]). Equivalence to fastwilcoxGMT is
# proven numerically by fcs_enrich_equivtest.R.
fcs_null_enrichstat_vectorized <- function(corStat, gmts, realenrich, num_g = 10) {
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
      stat[(n1 < num_g) | (n2 <= 2), ] <- NA_real_
    } else {
      notNA <- !is.na(sub)
      n1    <- as.matrix(M %*% (notNA * 1.0))              # sets x N
      ntot  <- matrix(colSums(notNA), nrow = nrow(M), ncol = ncol(sub), byrow = TRUE)
      n2    <- ntot - n1
      U     <- ranksum - (n1 * (n1 + 1) / 2)
      stat  <- (U / (n1 * n2)) - 0.5
      stat[(n1 < num_g) | (n2 <= 2)] <- NA_real_
    }
    rownames(stat) <- set_names
    enrichStat[[db]] <- stat
  }
  enrichStat
}

# Empirical permulation p-value per set (pseudo-count built in), per-ranking sidedness:
#   greater   -> (#{null >= obs} + 1) / (N_valid + 1)   [magnitude rankings]
#   two.sided -> (#{|null| >= |obs|} + 1) / (N_valid + 1) [signed rankings]
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

# rankings: named list of named-numeric vectors (already zero-floored).
fcs_run_all <- function(rankings, gmts, num_g = 10, perms_file = "NO_FILE") {
  res <- list(); alts <- list()
  for (rk in names(rankings)) {
    alts[[rk]] <- fcs_alternative(rankings[[rk]])
    r <- fcs_run_ranking(rankings[[rk]], gmts, num_g = num_g, alternative = alts[[rk]])
    if (nrow(r) > 0) res[[rk]] <- dplyr::mutate(r, ranking = rk)
  }
  if (length(res) == 0) {
    return(tibble::tibble(ranking = character(), database = character(), pathway = character(),
                          stat = numeric(), pval = numeric(), p.adj = numeric(), p.perm = numeric(),
                          num.genes = numeric(), gene.vals = character()))
  }
  enrich_df <- dplyr::bind_rows(res)
  enrich_df$p.perm <- NA_real_

  if (!is.null(perms_file) && perms_file != "NO_FILE" && file.exists(perms_file)) {
    fcs_progress(paste0("Loading null permulations from: ", perms_file))
    corperms <- tryCatch(readRDS(perms_file), error = function(e) NULL)

    # Two perms-RDS shapes are supported:
    #   • RER  : corperms$corStat (genes×N) [+ optional corRho] — ONE matrix shared
    #            across rankings; accel/decel derived from corRho per ranking.
    #   • CAAS : corperms$corStat_byrank, a named list keyed by ranking name
    #            (global/top/bottom), each a genes×N null matrix. The matching
    #            matrix is selected per ranking — directionality is precomputed
    #            upstream (each labeling's change_side partitions the directions),
    #            so no corRho transform is needed.
    corStat_byrank <- corperms$corStat_byrank
    if (!is.null(corperms) && (!is.null(corperms$corStat) || !is.null(corStat_byrank))) {
      base_corStat <- if (!is.null(corperms$corStat)) as.matrix(corperms$corStat) else NULL
      base_corRho  <- if (!is.null(corperms$corRho)) as.matrix(corperms$corRho) else NULL
      n_perms      <- if (!is.null(base_corStat)) ncol(base_corStat)
                      else if (length(corStat_byrank)) ncol(corStat_byrank[[1]]) else 0L
      fcs_progress(sprintf("Vectorized pathway permulations: N=%d, %d GMTs, %d rankings",
                           n_perms, length(gmts), length(rankings)))

      for (rk in names(rankings)) {
        obs_rk <- enrich_df %>% dplyr::filter(ranking == rk)
        if (nrow(obs_rk) == 0) next
        t0  <- Sys.time()
        alt <- alts[[rk]]

        # Per-ranking null matrix. CAAS supplies a precomputed per-direction matrix;
        # RER shares one base matrix and derives the accel/decel transform from corRho.
        corStat_rk <- NULL
        if (!is.null(corStat_byrank) && !is.null(corStat_byrank[[rk]])) {
          corStat_rk <- as.matrix(corStat_byrank[[rk]])
        } else if (!is.null(base_corStat)) {
          corStat_rk <- base_corStat
          if (rk == "accelerating" && !is.null(base_corRho)) {
            corStat_rk <- ifelse(base_corRho > 0,  base_corStat, 0)
          } else if (rk == "decelerating" && !is.null(base_corRho)) {
            corStat_rk <- ifelse(base_corRho < 0, -base_corStat, 0)
          }
          rownames(corStat_rk) <- rownames(base_corStat)
        }
        if (is.null(corStat_rk)) next  # no null for this ranking → p.perm stays NA

        # Observed stat per set, per db (rows = pathways).
        realenrich <- list()
        for (db in unique(obs_rk$database)) {
          db_df <- obs_rk %>% dplyr::filter(database == db) %>% as.data.frame()
          rownames(db_df) <- db_df$pathway
          realenrich[[db]] <- db_df[, c("pval", "stat"), drop = FALSE]
        }

        enrichStat <- tryCatch(
          fcs_null_enrichstat_vectorized(corStat_rk, gmts, realenrich, num_g = num_g),
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

  enrich_df %>% dplyr::relocate(ranking, database, pathway)
}

# ── Per-gene percentile flags (directional) ──────────────────────────────────
# From score_top / score_bottom columns, derive cumulative top-X% membership
# among the NON-ZERO (signal) genes of each directional ranking:
#   pct_top10  = gene is in the top 10% of the top ranking      (⊇ top5 ⊇ top1)
#   pct_bottom10 = gene is in the top 10% of the bottom ranking (strongly bottom)
# Returns a per-gene data.frame of logical columns (gene + pct_*).
fcs_percentile_flags <- function(stats, breaks = c(0.10, 0.05, 0.01)) {
  pct_of <- function(v) {
    ok <- !is.na(v) & v > 0
    p <- rep(NA_real_, length(v))
    if (any(ok)) p[ok] <- 1 - (rank(v[ok], ties.method = "max") - 1) / sum(ok)
    p
  }
  out <- data.frame(gene = stats$gene, stringsAsFactors = FALSE)
  if ("score_top" %in% names(stats)) {
    pt <- pct_of(stats$score_top)
    out$pct_top10 <- !is.na(pt) & pt <= breaks[1]
    out$pct_top5  <- !is.na(pt) & pt <= breaks[2]
    out$pct_top1  <- !is.na(pt) & pt <= breaks[3]
  }
  if ("score_bottom" %in% names(stats)) {
    pb <- pct_of(stats$score_bottom)
    out$pct_bottom10 <- !is.na(pb) & pb <= breaks[1]
    out$pct_bottom5  <- !is.na(pb) & pb <= breaks[2]
    out$pct_bottom1  <- !is.na(pb) & pb <= breaks[3]
  }
  out
}

# ── Leading-edge annotation ──────────────────────────────────────────────────
# Parse fastwilcoxGMT's "gene:rank, gene:rank, ..." gene.vals into one row per
# (ranking, database, pathway, gene) and left-join the per-gene attribute table
# (flag_* + pct_* logical columns).
fcs_annotate_leading_edge <- function(enrich_df, attr_df = NULL) {
  if (nrow(enrich_df) == 0) {
    cols <- c("ranking", "database", "pathway", "stat", "p.adj", "gene", "gene_rank")
    if (!is.null(attr_df)) {
      cols <- unique(c(cols, names(attr_df)))
    }
    empty_df <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
    colnames(empty_df) <- cols
    return(tibble::as_tibble(empty_df))
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

# ── Leading-edge SUMMARY ──────────────────────────────────────────────────────
fcs_leading_edge_summary <- function(le, attr_cols) {
  present <- intersect(attr_cols, names(le))
  if (nrow(le) == 0) {
    cols <- c("ranking", "database", "pathway", "stat", "p.adj", "n_le", present)
    empty_df <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
    colnames(empty_df) <- cols
    return(tibble::as_tibble(empty_df))
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
