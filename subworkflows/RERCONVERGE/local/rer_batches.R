# Bounded null generation: one batch of correlations, then gene slices for the
# native median-centred gene p. Public blocks contain only corStat for FCS.
rer_write_batches <- function(res, num_batches, perms_per_batch, seed, run_batch,
                              output, gene_batch_size = 1000L) {
  stopifnot(num_batches >= 1L, perms_per_batch >= 1L, gene_batch_size >= 1L,
            is.finite(seed), seed >= 0, seed + num_batches - 1 <= .Machine$integer.max)
  if (file.exists(output)) stop('RER null output already exists; use a clean task directory')
  dir.create(output, recursive=TRUE)
  scratch <- paste0(output, '.gene_slices')
  dir.create(scratch)
  genes <- rownames(res)
  groups <- split(seq_along(genes), ceiling(seq_along(genes)/gene_batch_size))
  settings <- attr(res, 'rer_settings')
  for (i in seq_len(num_batches)) {
    batch_seed <- as.integer(seed + i - 1L)
    set.seed(batch_seed)
    message(sprintf('[RER] batch %d/%d, seed=%d', i, num_batches, batch_seed))
    batch <- run_batch(perms_per_batch)
    stat <- as.matrix(batch$corStat)
    stopifnot(identical(rownames(stat), genes), ncol(stat) == perms_per_batch)
    colnames(stat) <- sprintf('rer_%06d', (i-1L)*perms_per_batch + seq_len(perms_per_batch))
    block <- list(corStat=stat, batch=i, seed=batch_seed, rer_settings=settings)
    saveRDS(block, file.path(output, sprintf('blk_%06d.rds', i)))
    for (k in seq_along(groups)) {
      saveRDS(as.matrix(batch$corRho[groups[[k]],,drop=FALSE]),
              file.path(scratch, sprintf('genes_%06d_block_%06d.rds', k, i)))
    }
    rm(batch, stat, block)
  }
  res$p.perm <- NA_real_
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    parts <- lapply(seq_len(num_batches), function(i)
      readRDS(file.path(scratch, sprintf('genes_%06d_block_%06d.rds', k, i))))
    rho <- do.call(cbind, parts)
    if (ncol(rho) == 1L) rho <- cbind(rho, NA_real_)
    ppc <- RERconverge::permpvalcor(res[idx,,drop=FALSE], list(corRho=as.data.frame(rho)))
    if (is.data.frame(ppc)) res$p.perm[idx] <- ppc$permpval
    else res$p.perm[idx] <- (as.numeric(ppc)*(num_batches*perms_per_batch)+1)/(num_batches*perms_per_batch+1)
  }
  res$p.perm.adj <- p.adjust(res$p.perm, method='BH')
  attr(res, 'n_perms') <- num_batches * perms_per_batch
  attr(res, 'rer_batch_seed') <- seed
  attr(res, 'rer_perms_per_batch') <- perms_per_batch
  unlink(scratch, recursive=TRUE)
  res
}
