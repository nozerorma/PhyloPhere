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
fcs_run_ranking <- function(vals, gmts, num_g = 10, enrich = FALSE) {
  reslist <- tryCatch(
    RERconverge::fastwilcoxGMTall(vals, gmts, outputGeneVals = TRUE, num.g = num_g),
    error = function(e) { message(sprintf("  fastwilcoxGMTall failed: %s", e$message)); NULL }
  )
  if (is.null(reslist) || length(reslist) == 0) return(tibble::tibble())
  out <- list()
  for (db in names(reslist)) {
    res <- reslist[[db]]
    if (is.null(res) || nrow(res) == 0) next
    out[[db]] <- tibble::as_tibble(res, rownames = "pathway") %>%
      dplyr::mutate(database = db)
  }
  if (length(out) == 0) return(tibble::tibble())
  dplyr::bind_rows(out) %>%
    dplyr::relocate(database, pathway) %>%
    dplyr::arrange(p.adj, dplyr::desc(abs(stat)))
}

# Convenience: run several named rankings (e.g. global/top/bottom) and tag each.
# Parallelized getEnrichPerms and vectorized permpvalenrich
fcs_get_enrich_perms_parallel <- function(corperms, realenrich, annotlist, ncores = 1) {
  numperms <- ncol(corperms$corP)
  groups <- length(realenrich)
  
  c <- 1
  while (c <= groups) {
    current <- realenrich[[c]]
    realenrich[[c]] <- current[order(rownames(current)), ]
    c <- c + 1
  }
  
  run_one_perm <- function(count) {
    statvec <- setNames(corperms$corStat[, count], rownames(corperms$corP))
    statvec <- na.omit(statvec)
    enrich <- RERconverge::fastwilcoxGMTall(statvec, annotlist, outputGeneVals = FALSE)
    
    res_pvals <- vector("list", groups)
    res_stats <- vector("list", groups)
    names(res_pvals) <- names(realenrich)
    names(res_stats) <- names(realenrich)
    
    for (db in names(realenrich)) {
      current <- enrich[[db]]
      if (!is.null(current) && nrow(current) > 0) {
        current <- current[order(rownames(current)), ]
        ref_names <- rownames(realenrich[[db]])
        current <- current[match(ref_names, rownames(current)), ]
        res_pvals[[db]] <- current$pval
        res_stats[[db]] <- current$stat
      } else {
        res_pvals[[db]] <- rep(NA_real_, nrow(realenrich[[db]]))
        res_stats[[db]] <- rep(NA_real_, nrow(realenrich[[db]]))
      }
    }
    list(pvals = res_pvals, stats = res_stats)
  }
  
  message(sprintf("[FCS] Running pathway permutations parallelized on %d cores...", ncores))
  results <- parallel::mclapply(1:numperms, run_one_perm, mc.cores = ncores)
  
  permenrichP <- vector("list", groups)
  permenrichStat <- vector("list", groups)
  names(permenrichP) <- names(realenrich)
  names(permenrichStat) <- names(realenrich)
  c <- 1
  while (c <= groups) {
    newdf <- data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
    rownames(newdf) <- rownames(realenrich[[c]])
    permenrichP[[c]] <- newdf
    permenrichStat[[c]] <- newdf
    c <- c + 1
  }
  
  for (count in 1:numperms) {
    perm_res <- results[[count]]
    if (is.null(perm_res)) {
      stop(paste("Permutation", count, "failed or returned NULL. Parallel execution might have encountered an error or OOM."))
    }
    for (c in 1:groups) {
      db <- names(realenrich)[c]
      permenrichP[[c]][, count] <- perm_res$pvals[[db]]
      permenrichStat[[c]][, count] <- perm_res$stats[[db]]
    }
  }
  
  data <- vector("list", 5)
  data[[1]] <- corperms$corP
  data[[2]] <- corperms$corRho
  data[[3]] <- corperms$corStat
  data[[4]] <- permenrichP
  data[[5]] <- permenrichStat
  names(data) <- c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  data
}

fcs_permpvalenrich_vectorized <- function(realenrich, permvals) {
  groups <- length(realenrich)
  enrichpvals <- vector("list", groups)
  names(enrichpvals) <- names(realenrich)
  
  for (count in 1:groups) {
    currreal <- realenrich[[count]]
    currenrich <- permvals$enrichStat[[count]]
    
    currreal <- currreal[match(rownames(currenrich), rownames(currreal)), , drop = FALSE]
    mat_enrich <- as.matrix(currenrich)
    obs_stat <- currreal$stat
    
    greater_mat <- abs(mat_enrich) > abs(obs_stat)
    lessnum <- rowSums(greater_mat, na.rm = TRUE)
    denom <- rowSums(!is.na(mat_enrich))
    pvals <- lessnum / denom
    pvals[is.na(obs_stat)] <- NA_real_
    
    names(pvals) <- rownames(currreal)
    enrichpvals[[count]] <- pvals
  }
  enrichpvals
}

# rankings: named list of named-numeric vectors (already zero-floored).
fcs_run_all <- function(rankings, gmts, num_g = 10, perms_file = "NO_FILE") {
  res <- list()
  for (rk in names(rankings)) {
    r <- fcs_run_ranking(rankings[[rk]], gmts, num_g = num_g)
    if (nrow(r) > 0) res[[rk]] <- dplyr::mutate(r, ranking = rk)
  }
  if (length(res) == 0) {
    return(tibble::tibble(ranking = character(), database = character(), pathway = character(),
                          stat = numeric(), p.val = numeric(), p.adj = numeric(), p.perm = numeric(), gene.vals = character()))
  }
  enrich_df <- dplyr::bind_rows(res)
  enrich_df$p.perm <- NA_real_
  
  # Check if perms_file is provided and exists
  if (!is.null(perms_file) && perms_file != "NO_FILE" && file.exists(perms_file)) {
    message("[FCS] Loading null permutations from: ", perms_file)
    corperms <- tryCatch(readRDS(perms_file), error = function(e) NULL)
    
    if (!is.null(corperms) && !is.null(corperms$corStat)) {
      n_perms <- ncol(corperms$corStat)
      message("[FCS] Running pathway permutations (N = ", n_perms, ")...")
      
      ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
      if (is.na(ncores) || ncores < 1) {
        ncores <- parallel::detectCores()
        if (is.na(ncores) || ncores < 1) ncores <- 1
      }
      ncores <- min(ncores, 16)
      
      for (rk in names(rankings)) {
        obs_rk <- enrich_df %>% dplyr::filter(ranking == rk)
        if (nrow(obs_rk) == 0) next
        
        # Reconstruct realenrich list format for getEnrichPerms
        realenrich <- list()
        for (db in unique(obs_rk$database)) {
          db_df <- obs_rk %>% dplyr::filter(database == db) %>% as.data.frame()
          rownames(db_df) <- db_df$pathway
          realenrich[[db]] <- db_df[, c("pval", "stat")]
          colnames(realenrich[[db]]) <- c("pval", "stat")
        }
        
        # Build null statistics matrix scaled/floored for the current ranking
        corperms_rk <- corperms
        if (rk == "global") {
          corperms_rk$corStat <- as.matrix(corperms$corStat)
        } else if (rk == "accelerating") {
          corperms_rk$corStat <- ifelse(as.matrix(corperms$corRho) > 0, as.matrix(corperms$corStat), 0)
        } else if (rk == "decelerating") {
          corperms_rk$corStat <- ifelse(as.matrix(corperms$corRho) < 0, -as.matrix(corperms$corStat), 0)
        } else {
          corperms_rk$corStat <- as.matrix(corperms$corStat)
        }
        
        permvals <- tryCatch({
          fcs_get_enrich_perms_parallel(corperms_rk, realenrich, gmts, ncores = ncores)
        }, error = function(e) {
          message(sprintf("[FCS] fcs_get_enrich_perms_parallel failed for ranking %s: %s", rk, e$message))
          NULL
        })
        
        if (!is.null(permvals)) {
          enrichpvals <- tryCatch({
            fcs_permpvalenrich_vectorized(realenrich, permvals)
          }, error = function(e) {
            message(sprintf("[FCS] fcs_permpvalenrich_vectorized failed for ranking %s: %s", rk, e$message))
            NULL
          })
          
          if (!is.null(enrichpvals)) {
            for (db in names(enrichpvals)) {
              raw_pvals <- enrichpvals[[db]]
              # Apply standard pseudo-count correction (num + 1) / (denom + 1) to pathway permulation p-values
              corrected_pvals <- (raw_pvals * n_perms + 1) / (n_perms + 1)
              
              match_idx <- which(enrich_df$ranking == rk & enrich_df$database == db)
              if (length(match_idx) > 0) {
                p_paths <- enrich_df$pathway[match_idx]
                enrich_df$p.perm[match_idx] <- corrected_pvals[p_paths]
              }
            }
          }
        }
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
