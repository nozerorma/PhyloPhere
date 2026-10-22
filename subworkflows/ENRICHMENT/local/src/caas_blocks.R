# caas_blocks.R — the CAAS null by blocks of cycles and batches of genes (WP6: transposition and GENE_P per batch).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
#
# Input: per block the matrices.rds of CAAS_CORE_MERGE (list(genes, cycles, gene_stat, stat = list(global, top, bottom))),
# genes x cycles of that block. The universe U (gene order fixed by the caller) and the batch size G define the batches:
# batch k = U[((k-1)G + 1) .. min(kG, |U|)]. Every step works on one block or one batch at a time.
#
#   caas_transpose_block()   block -> one slice per batch (rows of the batch, cycles of the block), zeros for absent genes
#   caas_gene_p_batch()      batch -> reads the slices of every block, GENE_P (gene_p.R) on the full cycle set, writes the
#                            -log10 p slices per block and the observed -log10 p of the batch
#   caas_assemble_block()    block -> genes x cycles matrix of -log10 p, rows in U order (input of the FCS counts)
#
# Slices are RDS lists: list(rows = genes, cols = cycles, values = matrix). The null cycle order inside a block is the
# order of its matrices.rds; the cycle identity is kept by name, so the blocks may come in any order.

caas_batches <- function(universe, batch_size) {
  stopifnot(batch_size >= 1)
  split(universe, ceiling(seq_along(universe) / batch_size))
}

# Slice file of one (ranking, batch, block).
caas_slice_path <- function(dir, rk, k, block_id) file.path(dir, sprintf("%s_b%04d_%s.rds", rk, k, block_id))

# Writes one slice per batch of the block `stat_mat` (genes x cycles, row names = genes; any gene of U may be absent).
caas_transpose_block <- function(stat_mat, batches, rk, block_id, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(rownames(stat_mat)) || anyNA(rownames(stat_mat)) || anyDuplicated(rownames(stat_mat))) stop('Invalid block gene IDs')
  cols <- colnames(stat_mat)
  if (is.null(cols) || anyNA(cols) || anyDuplicated(cols) || any(!nzchar(cols))) stop('Invalid block cycle IDs')
  for (k in seq_along(batches)) {
    rows <- batches[[k]]
    vals <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
    common <- intersect(rows, rownames(stat_mat))
    vals[common, ] <- stat_mat[common, , drop = FALSE]
    saveRDS(list(rows = rows, cols = cols, values = vals), caas_slice_path(dir, rk, k, block_id))
  }
  invisible(length(batches))
}

# GENE_P of batch k of ranking rk. `block_ids` = the blocks of the null; `obs` = observed gene statistic over the batch
# (named vector, the rows of the batch). Writes the -log10 p slices per block and returns the observed -log10 p.
caas_gene_p_batch <- function(k, rk, block_ids, dir, obs, tie_tol = TIE_TOL_GENE_P) {
  parts <- lapply(block_ids, function(b) readRDS(caas_slice_path(dir, rk, k, b)))
  rows <- parts[[1]]$rows
  if (anyDuplicated(block_ids) || !all(vapply(parts, function(p) identical(p$rows, rows) &&
      identical(rownames(p$values), rows) && identical(colnames(p$values), p$cols), logical(1)))) stop('Inconsistent gene batch slices')
  cycles <- unlist(lapply(parts, function(p) p$cols), use.names=FALSE)
  if (anyDuplicated(cycles)) stop('Overlapping null cycles across blocks')
  S <- do.call(cbind, lapply(parts, function(p) p$values))
  cols <- unlist(lapply(parts, function(p) p$cols))
  colnames(S) <- cols
  obs_b <- obs[rows]; obs_b[is.na(obs_b)] <- 0
  g <- gene_p(S, obs_b, tie_tol = tie_tol)
  for (i in seq_along(block_ids)) {
    sel <- match(parts[[i]]$cols, colnames(S))
    out <- list(rows = rows, cols = parts[[i]]$cols, values = g$nlog10_null[, sel, drop = FALSE])
    saveRDS(out, caas_slice_path(dir, paste0(rk, "_gp"), k, block_ids[i]))
  }
  setNames(g$nlog10_obs, rows)
}

# -log10 p matrix of one block (genes x cycles of the block), rows in the order of `universe`.
caas_assemble_block <- function(k_all, rk, block_id, dir, universe) {
  pieces <- lapply(k_all, function(k) readRDS(caas_slice_path(dir, paste0(rk, "_gp"), k, block_id)))
  cols <- pieces[[1]]$cols
  M <- matrix(NA_real_, length(universe), length(cols), dimnames = list(universe, cols))
  for (p in pieces) M[p$rows, ] <- p$values
  M
}
