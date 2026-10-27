#!/usr/bin/env Rscript
# term_blocks.R — Thematic blocks of enriched terms (FCS gene sets and AMI module terms).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: 12.FCS_general_report.Rmd and 13.AMI_analysis.Rmd. Defines functions only;
# it needs no package.
#
# Terms are grouped after inference; p/q and calls are never modified. All terms remain
# in the exports. The distance between two terms is 1 - Jaccard of their gene sets and
# the tree is built with complete linkage and cut by height, so every pair inside a block has
# Jaccard >= thr. The number of blocks is not bounded: cutting a complete-linkage tree by count
# isolates the peripheral terms first and leaves one very large, heterogeneous block. Callers
# that need a short list show the first blocks (numbered by their best member) and keep the
# rest in the tables.

# Pairwise Jaccard matrix of a named list of character vectors. Two empty sets have no
# similarity (0): a set without genes is its own block.
term_jaccard <- function(sets) {
  n <- length(sets)
  if (n == 0) return(matrix(numeric(0), 0, 0))
  genes <- unique(unlist(sets, use.names = FALSE))
  inc <- matrix(0, nrow = length(genes), ncol = n)
  for (j in seq_len(n)) inc[match(unique(sets[[j]]), genes), j] <- 1
  inter <- crossprod(inc)
  sz <- diag(inter)
  uni <- outer(sz, sz, "+") - inter
  J <- ifelse(uni > 0, inter / uni, 0)
  diag(J) <- 1
  dimnames(J) <- list(names(sets), names(sets))
  J
}

# sets: named list (name = term key) of gene vectors.
# rank: numeric, lower is better (e.g. q or FDR); ties break by name, so the result does not
#   depend on the order of the input rows.
# thr: minimum Jaccard between any two terms of a block.
# Returns a data.frame in the input order: key, block (integer, 1 = block of the best term),
# is_representative (the best-ranked member of its block), block_size, block_min_jaccard
# (NA for a block of one term).
term_blocks <- function(sets, rank = seq_along(sets), thr = 0.2) {
  n <- length(sets)
  keys <- if (is.null(names(sets))) as.character(seq_len(n)) else names(sets)
  empty <- data.frame(key = character(0), block = integer(0), is_representative = logical(0),
                      block_size = integer(0), block_min_jaccard = numeric(0),
                      stringsAsFactors = FALSE)
  if (n == 0) return(empty)
  stopifnot(length(rank) == n, is.finite(thr), thr >= 0, thr <= 1)
  ord <- order(rank, keys, method = "radix")
  s <- sets[ord]; names(s) <- keys[ord]
  J <- term_jaccard(s)
  if (n == 1) {
    cl <- 1L
  } else {
    hc <- hclust(as.dist(1 - J), method = "complete")
    # The tolerance keeps a pair at exactly thr inside the block despite rounding.
    cl <- cutree(hc, h = 1 - thr + 1e-9)
  }
  # Blocks are numbered by their best member, which is also the representative.
  block <- match(cl, unique(cl))
  rep_flag <- !duplicated(block)
  size <- as.integer(table(block)[as.character(block)])
  mnJ <- vapply(seq_len(max(block)), function(b) {
    ix <- which(block == b)
    if (length(ix) < 2) NA_real_ else min(J[ix, ix][upper.tri(J[ix, ix])])
  }, numeric(1))
  res <- data.frame(key = names(s), block = block, is_representative = rep_flag,
                    block_size = size, block_min_jaccard = mnJ[block], stringsAsFactors = FALSE)
  res[match(keys, res$key), , drop = FALSE]
}

# Adds cluster_id, is_representative, cluster_size and cluster_min_jaccard to the results
# table of fcs_run_all().
# enrich: columns ranking, database, pathway, stat, gene.vals ("gene:rank" of every member of
#   the set in the universe), q.perm, p.perm and selection_col.
# Groups are formed per ranking AND effect sign over the rows with selection_col TRUE, with
# the full gene set of each term (the leading edge describes what drives a term, not what the
# term is about). callers may supply a shared presentation order through order_fun.
# Other rows get NA; p/q and selection calls are never changed.
# cluster_id is "<ranking>_<effect>_C<k>", k following the best-member order.
fcs_add_clusters <- function(enrich, threshold = 0.2, selection_col = "sig_perm", order_fun = NULL) {
  nr <- nrow(enrich)
  enrich$cluster_id <- rep(NA_character_, nr)
  enrich$is_representative <- rep(NA, nr)
  enrich$cluster_size <- rep(NA_integer_, nr)
  enrich$cluster_min_jaccard <- rep(NA_real_, nr)
  if (nr == 0) return(enrich)
  members <- function(x) {
    g <- if (is.na(x) || !nzchar(x)) character(0) else sub(":.*$", "", trimws(strsplit(x, ",\\s*")[[1]]))
    unique(g[nzchar(g)])
  }
  strata <- unique(data.frame(ranking = enrich$ranking, effect = sign(enrich$stat)))
  for (j in seq_len(nrow(strata))) {
    rk <- strata$ranking[j]; effect <- strata$effect[j]
    if (!is.finite(effect)) next
    idx <- which(enrich$ranking == rk & sign(enrich$stat) == effect & enrich[[selection_col]] %in% TRUE)
    if (length(idx) == 0) next
    d <- enrich[idx, , drop = FALSE]
    o <- if (is.null(order_fun)) order(d$q.perm, d$p.perm, -abs(d$stat), d$database, d$pathway) else {
      d$.source_index <- seq_len(nrow(d)); order_fun(d)$.source_index
    }
    idx <- idx[o]
    keys <- paste(enrich$database[idx], enrich$pathway[idx], sep = "\r")
    sets <- lapply(enrich$gene.vals[idx], members); names(sets) <- keys
    # idx is already in presentation order, so the position is the rank.
    bl <- term_blocks(sets, rank = seq_along(idx), thr = threshold)
    enrich$cluster_id[idx] <- sprintf("%s_%s_C%03d", rk,
                                      if (effect > 0) "positive" else if (effect < 0) "negative" else "zero", bl$block)
    enrich$is_representative[idx] <- bl$is_representative
    enrich$cluster_size[idx] <- bl$block_size
    enrich$cluster_min_jaccard[idx] <- bl$block_min_jaccard
  }
  enrich
}
