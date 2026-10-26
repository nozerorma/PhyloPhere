#!/usr/bin/env Rscript
# fcs_cluster.R — Redundancy reduction of the significant gene sets of an FCS ranking.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
# =============================================================================
# Sourced by: 12.FCS_general_report.Rmd. Defines functions only; it needs no package.
#
# Group selected sets after inference; p/q and calls are never modified.
# The report selects C>0 and supplies its deterministic C/checks/q/p order.
# Group independently by ranking and effect sign. A term joins the first
# representative with leading-edge Jaccard >= threshold; comparisons against
# group members are forbidden (no chaining). All terms remain in the exports.
# Overlap coefficient remains an explicit alternative, not the default; it
# over-merges hierarchical collections. Pre-test composition deduplication is
# not performed because it would change the inferential family.

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
# Groups are formed per ranking AND effect sign over selection_col; callers may supply a shared presentation order.
# Other rows get NA; p/q and selection calls are never changed.
# cluster_id is "<ranking>_<effect>_C<k>", k following the representative order.
fcs_add_clusters <- function(enrich, leading, threshold = 0.5, method = "jaccard", selection_col = "sig_perm", order_fun = NULL) {
  enrich$cluster_id <- rep(NA_character_,nrow(enrich))
  enrich$is_representative <- rep(NA,nrow(enrich))
  key <- function(d) paste(d$ranking, d$database, d$pathway, sep = "\r")
  keys <- key(enrich)
  le_key <- key(leading)
  genes_by_key <- split(as.character(leading$gene), le_key)
  strata <- unique(data.frame(ranking=enrich$ranking, effect=sign(enrich$stat)))
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
    sets <- lapply(keys[idx], function(k) { g <- genes_by_key[[k]]; if (is.null(g)) character(0) else unique(g) })
    names(sets) <- keys[idx]
    cl <- fcs_cluster_sets(sets, threshold = threshold, method = method)
    enrich$cluster_id[idx] <- sprintf("%s_%s_C%03d", rk, if(effect>0)"positive" else if(effect<0)"negative" else "zero", cl$cluster)
    enrich$is_representative[idx] <- cl$is_representative
  }
  enrich
}
