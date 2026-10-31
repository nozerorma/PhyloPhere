# Batches and a content-verified roster. Native median-centred calibration is unchanged.
.rer_batch_script <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if (!exists('rer_file_hash', mode='function')) source(file.path(dirname(normalizePath(.rer_batch_script)), 'rer_contract.R'))
# Bounded null generation: one batch of correlations, then gene slices for the
# native median-centred gene p. Public blocks contain only corStat for FCS.
rer_write_batches <- function(res, num_batches, perms_per_batch, seed, run_batch,
                              output, gene_batch_size = 1000L) {
  num_batches <- rer_integer(num_batches,'num_batches',1)
  perms_per_batch <- rer_integer(perms_per_batch,'perms_per_batch',1)
  gene_batch_size <- rer_integer(gene_batch_size,'gene_batch_size',1)
  seed <- rer_integer(seed,'seed')
  stopifnot(num_batches >= 1L, perms_per_batch >= 1L, gene_batch_size >= 1L,
            is.finite(seed), seed >= 0, seed + num_batches - 1 <= .Machine$integer.max)
  if (file.exists(output)) stop('RER null output already exists; use a clean task directory')
  dir.create(output, recursive=TRUE)
  scratch <- paste0(output, '.gene_slices')
  dir.create(scratch)
  genes <- rownames(res)
  rer_ids(genes, "Null genes")
  roster <- list()
  groups <- split(seq_along(genes), ceiling(seq_along(genes)/gene_batch_size))
  settings <- attr(res, 'rer_settings')
  for (i in seq_len(num_batches)) {
    batch_seed <- as.integer(seed + i - 1L)
    set.seed(batch_seed)
    message(sprintf('[RER] batch %d/%d, seed=%d', i, num_batches, batch_seed))
    batch <- run_batch(perms_per_batch)
    stat <- as.matrix(batch$corStat)
    rho <- as.matrix(batch$corRho)
    if (!is.numeric(rho) || !identical(rownames(rho),genes) || ncol(rho)!=perms_per_batch ||
        any(abs(rho)>1,na.rm=TRUE)) stop('Invalid null rho batch')
    stopifnot(identical(rownames(stat), genes), ncol(stat) == perms_per_batch)
    colnames(stat) <- sprintf('rer_%06d', (i-1L)*perms_per_batch + seq_len(perms_per_batch))
    block <- list(corStat=stat, batch=i, seed=batch_seed, rer_settings=settings)
    saveRDS(block, file.path(output, sprintf('blk_%06d.rds', i)))
    roster[[i]] <- list(file=sprintf('blk_%06d.rds',i),sha256=rer_file_hash(file.path(output,sprintf('blk_%06d.rds',i))),cycles=colnames(stat))
    for (k in seq_along(groups)) {
      saveRDS(as.matrix(batch$corRho[groups[[k]],,drop=FALSE]),
              file.path(scratch, sprintf('genes_%06d_block_%06d.rds', k, i)))
    }
    rm(batch, stat, rho, block)
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
    else stop('Unsupported permpvalcor return schema')
  }
  res$p.perm.adj <- p.adjust(res$p.perm, method='BH')
  attr(res, 'n_perms') <- num_batches * perms_per_batch
  attr(res, 'rer_batch_seed') <- seed
  attr(res, 'rer_perms_per_batch') <- perms_per_batch
  manifest <- list(version='rer_null_v1', genes=genes, n_perms=num_batches*perms_per_batch,
    observed_sha256=rer_hash(list(genes=genes,Rho=res$Rho,P=res$P)),
    settings=settings, blocks=roster)
  jsonlite::write_json(manifest,file.path(output,'rer_null_manifest.json'),auto_unbox=TRUE,null='null',digits=NA)
  attr(res,'rer_null_manifest') <- rer_file_hash(file.path(output,'rer_null_manifest.json'))
  unlink(scratch, recursive=TRUE)
  res
}
