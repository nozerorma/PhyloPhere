#!/usr/bin/env Rscript
# fcs_cluster.R — Redundancy reduction of the significant gene sets of an FCS ranking.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: 12.FCS_general_report.Rmd. Defines functions only; it needs no package.
#
# The gene sets of the databases overlap (a pathway and its parent term, the same set in two
# databases), so a list of significant sets repeats the same genes. The reduction is applied
# after the significance call and changes no p or q value: it groups the significant sets of
# one ranking by the similarity of their leading edges and names one representative per group.
#
# Assignment is greedy and deterministic. The significant sets are taken in order of
# significance (smallest q.perm first). The first one is a representative; each following set
# joins the first representative whose leading edge is at least `threshold` similar to its
# own, and becomes a representative itself when none is. A set is compared with the
# representatives only, never with the members of a group, so groups do not chain.
# =============================================================================

# Similarity of two leading edges (character vectors of genes).
#   jaccard : |a n b| / |a u b|
#   overlap : |a n b| / min(|a|, |b|); it collapses a set into any larger set that contains it,
#             so it over-merges hierarchical databases (Reactome, HPO).
# Two empty sets have no similarity (0): a set without leading edge is its own group.
fcs_set_similarity <- function(a, b, method = "jaccard") {
  inter <- length(intersect(a, b))
  den <- switch(method,
                jaccard = length(union(a, b)),
                overlap = min(length(a), length(b)),
                stop(sprintf("unknown similarity '%s' (jaccard or overlap)", method)))
  if (den == 0) 0 else inter / den
}

# sets: named list (name = set key) of leading-edge genes, in order of significance (best first).
# Returns a tibble: key, cluster (integer, 1 = group of the first set), is_representative.
fcs_cluster_sets <- function(sets, threshold = 0.5, method = "jaccard") {
  n <- length(sets)
  cluster <- integer(n)
  reps <- integer(0)
  for (i in seq_len(n)) {
    joined <- 0L
    for (k in seq_along(reps)) {
      if (fcs_set_similarity(sets[[i]], sets[[reps[k]]], method) >= threshold) { joined <- k; break }
    }
    if (joined == 0L) { reps <- c(reps, i); joined <- length(reps) }
    cluster[i] <- joined
  }
  data.frame(key = names(sets), cluster = cluster, is_representative = seq_len(n) %in% reps,
             stringsAsFactors = FALSE)
}

# Adds cluster_id and is_representative to the results table of fcs_run_all().
# enrich: the table (columns ranking, database, pathway, q.perm, p.perm, stat, sig_perm).
# leading: long table of leading-edge genes (ranking, database, pathway, gene).
# Groups are formed per ranking over selection_col; callers may supply a shared presentation order.
# Other rows get NA; p/q and selection calls are never changed.
# cluster_id is "<ranking>_C<k>", k following the significance of the representative.
fcs_add_clusters <- function(enrich, leading, threshold = 0.5, method = "jaccard", selection_col = "sig_perm", order_fun = NULL) {
  enrich$cluster_id <- NA_character_
  enrich$is_representative <- NA
  key <- function(d) paste(d$ranking, d$database, d$pathway, sep = "\r")
  keys <- key(enrich)
  le_key <- key(leading)
  genes_by_key <- split(as.character(leading$gene), le_key)
  for (rk in unique(enrich$ranking)) {
    idx <- which(enrich$ranking == rk & enrich[[selection_col]] %in% TRUE)
    if (length(idx) == 0) next
    d <- enrich[idx, , drop = FALSE]
    o <- if (is.null(order_fun)) order(d$q.perm, d$p.perm, -abs(d$stat), d$database, d$pathway) else {
      d$.source_index <- seq_len(nrow(d)); order_fun(d)$.source_index
    }
    idx <- idx[o]
    sets <- lapply(keys[idx], function(k) { g <- genes_by_key[[k]]; if (is.null(g)) character(0) else unique(g) })
    names(sets) <- keys[idx]
    cl <- fcs_cluster_sets(sets, threshold = threshold, method = method)
    enrich$cluster_id[idx] <- sprintf("%s_C%03d", rk, cl$cluster)
    enrich$is_representative[idx] <- cl$is_representative
  }
  enrich
}
