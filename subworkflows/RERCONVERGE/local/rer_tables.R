# Produce validated machine-readable products independently of HTML rendering.
a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
res <- rer_results(readRDS(a[1]))
nominal <- rer_number(a[4],'nominal threshold'); fdr_levels <- rer_bh_levels(a[5])
if (nominal < 0 || nominal > 1) stop('Nominal threshold must be in [0,1]')
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
      if (rer_verify_checksums() && (is.null(expected) || !identical(expected,rer_file_hash(manifest_path)))) stop('Null manifest does not belong to observed results')
      manifest <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
      if (!identical(unlist(manifest$genes,use.names=FALSE),rownames(res)) ||
          !identical(manifest$observed_sha256,rer_hash(list(genes=rownames(res),Rho=res$Rho,P=res$P)))) stop('Null/observed identity mismatch')
      files <- vapply(manifest$blocks,`[[`,'','file')
      if (!identical(basename(paths),files)) stop('Incomplete or extra null blocks')
      if (rer_verify_checksums()) for (i in seq_along(paths)) {
        checksum <- manifest$blocks[[i]]$sha256
        if (!is.null(checksum) && !identical(rer_file_hash(paths[i]),checksum)) stop('Null block checksum mismatch')
      }
      null_status <- if(rer_verify_checksums()) 'verified_with_available_checksums' else 'verified_scientific_identity_and_roster'
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
    if (!is.null(manifest) && !identical(b$rer_settings, attr(res,'rer_settings')))
      stop('Null block correlation settings differ from observed results')
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
# Resolve producer IDs only after validating the untouched observed/null pair.
# The CAAS background (or alignment basenames) is the reference convention;
# unmatched RER genes retain their original ID. Numeric gene suffixes are literal.
original_genes <- rownames(res)
reference_path <- if (length(a)>=6L) a[6] else 'NO_FILE'
if (basename(reference_path) != 'NO_FILE') {
  if (!file.exists(reference_path)) stop('Missing gene identifier reference')
  if (length(a)<7L || !file.exists(a[7])) stop('Missing shared gene identifier resolver')
  source(a[7])
  if (dir.exists(reference_path)) {
    fs <- list.files(reference_path, full.names=TRUE)
    fs <- fs[!file.info(fs)$isdir & !grepl('\\.(txt|tsv|csv|log|map)$',fs,ignore.case=FALSE)]
    reference <- tools::file_path_sans_ext(basename(fs))
  } else {
    reference <- readLines(reference_path,warn=FALSE)
    if (length(reference) && reference[1] %in% c('Gene','gene')) reference <- reference[-1]
  }
  rer_ids(reference, 'CAAS reference genes')
  matched <- scoring_match_genes(original_genes, reference)
  canonical <- original_genes
  canonical[!is.na(matched)] <- reference[matched[!is.na(matched)]]
  rer_ids(canonical, 'RER resolved gene identifiers')
} else canonical <- original_genes
changed <- !identical(canonical,original_genes)
original_settings <- attr(res,'rer_settings')
original_null_manifest <- attr(res,'rer_null_manifest')
if (changed) {
  rownames(res) <- canonical
  settings <- original_settings
  if (basename(null_path) == 'NO_FILE') attr(res,'rer_null_manifest') <- NULL
  if (!is.null(settings$matrix_contract$genes)) {
    settings$matrix_contract$genes <- canonical
    attr(res,'rer_settings') <- settings
  }
}
translation <- data.frame(original_gene=original_genes,gene=canonical,stringsAsFactors=FALSE)
rer_write_tsv(translation,'rer_gene_identifiers.tsv')
# Materialize a coherent renamed copy of the already validated null only when
# names changed; the numerical values, cycles and original files are untouched.
if (basename(null_path) != 'NO_FILE') {
  if (!changed) {
    if (!file.symlink(normalizePath(null_path),'normalized_null')) stop('Cannot link RER null')
  } else {
    dir.create('normalized_null')
    for (i in seq_along(paths)) {
      b <- readRDS(paths[i]); m <- as.matrix(b$corStat)
      rownames(m) <- canonical
      if (is.null(colnames(m))) colnames(m) <- paste0('legacy_',seq_len(ncol(m)))
      b$corStat <- m
      if (!is.null(b$rer_settings)) b$rer_settings <- attr(res,'rer_settings')
      target <- file.path('normalized_null',basename(paths[i]))
      saveRDS(b,target)
      if (!is.null(manifest)) manifest$blocks[[i]]$sha256 <- if(rer_verify_checksums()) rer_file_hash(target) else NULL
    }
    if (!is.null(manifest)) {
      manifest$genes <- as.list(canonical)
      manifest$settings <- attr(res,'rer_settings')
      manifest$observed_sha256 <- rer_hash(list(genes=canonical,Rho=res$Rho,P=res$P))
      jsonlite::write_json(manifest,file.path('normalized_null','rer_null_manifest.json'),auto_unbox=TRUE,null='null',digits=NA)
      attr(res,'rer_null_manifest') <- rer_file_hash(file.path('normalized_null','rer_null_manifest.json'))
    }
  }
}
summary <- rer_summary(res,nominal,fdr_levels)
metadata <- list(version='rer_v1',traitname=a[3],analysis_mode=attr(res,'rer_type'),
  inference=attr(res,'rer_settings'),n_perms=attr(res,'n_perms'),seed=attr(res,'rer_batch_seed'),
  nominal_threshold=nominal,fdr_levels=fdr_levels,null_status=null_status,
  original_result_identity=rer_file_identity(a[1]),
  original_result_sha256=if(rer_verify_checksums()) rer_file_hash(a[1]) else NULL,null_manifest=manifest,
  identifier_translation=list(reference=if(basename(reference_path)=='NO_FILE') NULL else as.character(reference_path), changed=sum(original_genes!=canonical), original_settings_sha256=if(changed && !is.null(original_settings)) rer_hash(original_settings) else NULL, original_null_manifest_sha256=if(changed) original_null_manifest else NULL),
  inference_provenance=if(is.null(attr(res,'rer_settings'))) 'legacy_unknown' else 'recorded')
attr(res,'rer_report_metadata') <- metadata
saveRDS(res,'validated_rer.output')
rer_write_tsv(summary,paste0('rerconverge_summary_',a[3],'.tsv'))
jsonlite::write_json(metadata,'rer_run_metadata.json',auto_unbox=TRUE,null='null',digits=NA)
if (!dir.exists('rer_inputs')) dir.create('rer_inputs')
files <- c('validated_rer.output','rer_run_metadata.json','rer_gene_identifiers.tsv',paste0('rerconverge_summary_',a[3],'.tsv'))
if (!all(file.copy(files,'rer_inputs',overwrite=TRUE))) stop('Cannot preserve RER inputs')
identities <- lapply(file.path('rer_inputs',files),rer_file_identity)
names(identities) <- files
jsonlite::write_json(identities,'rer_inputs/file_identity.json',auto_unbox=TRUE)
if (rer_verify_checksums()) {
  hashes <- lapply(file.path('rer_inputs',files),rer_file_hash)
  names(hashes) <- files
  jsonlite::write_json(hashes,'rer_inputs/checksums.json',auto_unbox=TRUE)
}
