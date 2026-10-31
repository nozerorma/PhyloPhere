# Produce validated machine-readable products independently of HTML rendering.
a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
res <- rer_results(readRDS(a[1]))
nominal <- rer_number(a[4],'nominal threshold'); fdr <- rer_number(a[5],'FDR threshold')
if (nominal < 0 || nominal > 1 || fdr < 0 || fdr > 1) stop('Thresholds must be in [0,1]')
null_path <- a[2]
null_status <- 'unavailable'
manifest <- NULL
if (basename(null_path) != 'NO_FILE') {
  if (!file.exists(null_path)) stop('Missing null input')
  manifest_path <- file.path(null_path,'rer_null_manifest.json')
  if (dir.exists(null_path)) {
    paths <- sort(list.files(null_path,pattern='\\.rds$',full.names=TRUE))
    if (file.exists(manifest_path)) {
      expected <- attr(res,'rer_null_manifest')
      if (is.null(expected) || !identical(expected,rer_file_hash(manifest_path))) stop('Null manifest does not belong to observed results')
      manifest <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
      if (!identical(unlist(manifest$genes,use.names=FALSE),rownames(res)) ||
          !identical(manifest$observed_sha256,rer_hash(list(genes=rownames(res),Rho=res$Rho,P=res$P)))) stop('Null/observed identity mismatch')
      files <- vapply(manifest$blocks,`[[`,'','file')
      if (!identical(basename(paths),files)) stop('Incomplete or extra null blocks')
      for (i in seq_along(paths)) if (!identical(rer_file_hash(paths[i]),manifest$blocks[[i]]$sha256)) stop('Null block checksum mismatch')
      null_status <- 'verified'
    } else {
      if (!is.null(attr(res,'rer_null_manifest'))) stop('Missing required null manifest')
      null_status <- 'legacy_unverified_coverage'
    }
  } else {
    if (!is.null(attr(res,'rer_null_manifest'))) stop('New results require their original manifested null store')
    paths <- null_path
    null_status <- 'legacy_unverified_coverage'
  }
  if (!length(paths)) stop('Empty null input')
  seen <- character()
  for (i in seq_along(paths)) {
    b <- readRDS(paths[i]); m <- as.matrix(b$corStat)
    if (!is.numeric(m) || !ncol(m)) stop('Invalid null corStat')
    rer_ids(rownames(m),'Null gene identifiers')
    if (!identical(rownames(m),rownames(res))) stop('Null genes/order differ from observed results')
    ids <- colnames(m)
    if (length(paths)==1L && is.null(ids) && is.null(manifest)) ids <- paste0('legacy_',seq_len(ncol(m)))
    rer_ids(ids,'Null cycles')
    if (any(ids %in% seen)) stop('Overlapping null cycles')
    if (!is.null(manifest) && !identical(ids,unlist(manifest$blocks[[i]]$cycles,use.names=FALSE))) stop('Null cycle roster mismatch')
    seen <- c(seen,ids)
  }
  recorded <- attr(res,'n_perms')
  if (!is.null(recorded) && length(seen) != recorded) stop('Null cycle count differs from observed results')
}
summary <- rer_summary(res,nominal,fdr)
metadata <- list(version='rer_v1',traitname=a[3],analysis_mode=attr(res,'rer_type'),
  inference=attr(res,'rer_settings'),n_perms=attr(res,'n_perms'),seed=attr(res,'rer_batch_seed'),
  nominal_threshold=nominal,fdr_threshold=fdr,null_status=null_status,
  original_result_sha256=rer_file_hash(a[1]),null_manifest=manifest,
  inference_provenance=if(is.null(attr(res,'rer_settings'))) 'legacy_unknown' else 'recorded')
attr(res,'rer_report_metadata') <- metadata
saveRDS(res,'validated_rer.output')
rer_write_tsv(summary,paste0('rerconverge_summary_',a[3],'.tsv'))
jsonlite::write_json(metadata,'rer_run_metadata.json',auto_unbox=TRUE,null='null',digits=NA)
if (!dir.exists('rer_inputs')) dir.create('rer_inputs')
files <- c('validated_rer.output','rer_run_metadata.json',paste0('rerconverge_summary_',a[3],'.tsv'))
if (!all(file.copy(files,'rer_inputs',overwrite=TRUE))) stop('Cannot preserve RER inputs')
hashes <- lapply(file.path('rer_inputs',files),rer_file_hash)
names(hashes) <- files
jsonlite::write_json(hashes,'rer_inputs/checksums.json',auto_unbox=TRUE)
