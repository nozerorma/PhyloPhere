# Normalize null inputs without loading the live block matrices.
a <- commandArgs(TRUE)
source_type <- a[1]; gene_stat <- a[2]; block_size <- as.integer(a[3])
if (!identical(gene_stat, 'max')) stop('CAAS gene_stat mismatch: only max is supported')
dir.create('blocks')
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
contract <- file.path(dirname(normalizePath(script)), '../../../SCORING/local/src/caas_block_contract.R')
source(contract)
paths <- Sys.glob('live/*/matrices.rds')
stores <- Sys.glob('live/*/blocks')
if (source_type == 'caas' && a[4] != 'disabled' && length(Sys.glob('live/*')) && !length(paths) && !length(stores)) stop('Missing canonical CAAS matrix blocks')
if (source_type == 'caas' && a[4] != 'disabled' && length(stores)) {
  if (length(stores) != 1L || length(paths)) stop('Mixed CAAS stores')
  paths <- Sys.glob(file.path(stores, 'blk_*/matrices.rds'))
  roster_file <- file.path(dirname(stores), 'cycle_roster.txt')
  if (!file.exists(roster_file)) stop('Missing canonical CAAS cycle roster')
  roster <- readLines(roster_file)
  index <- caas_block_index(paths, roster=roster, gene_stat=gene_stat)
  meta <- file.path(dirname(stores), 'null_metadata.json')
  if (file.exists(meta)) {
    saved <- jsonlite::read_json(meta, simplifyVector=FALSE)
    if (!identical(saved$gene_stat, index$gene_stat) || !identical(as.character(unlist(saved$cycles)), index$cycles) ||
        !identical(as.character(unlist(saved$genes)), index$genes)) stop('CAAS metadata disagrees with blocks')
    hashes <- vapply(index$blocks, function(b) b$md5, '')
    old_hashes <- vapply(saved$blocks, function(b) b$md5, '')
    if (!identical(sort(hashes), sort(old_hashes))) stop('CAAS block content differs from metadata')
  }
}
if (a[4] == 'disabled') { paths <- character() } else if (!length(paths) && dir.exists('null_input')) {
  paths <- if (source_type == 'caas') Sys.glob('null_input/blocks/blk_*/matrices.rds') else sort(list.files('null_input', pattern = '\\.rds$', full.names = TRUE, recursive = TRUE))
  if (source_type == 'caas') {
    if (!file.exists('null_input/cycle_roster.txt')) stop('Missing canonical CAAS cycle roster')
    caas_block_index(paths, readLines('null_input/cycle_roster.txt'), gene_stat=gene_stat)
  }
}
if (a[4] != 'disabled' && !length(paths) && file.exists('null_input') && !dir.exists('null_input')) {
  obj <- readRDS('null_input')
  if (source_type == 'caas') {
    if (is.null(obj$gene_stat)) stop('CAAS null lacks gene_stat stamp; rebuild gene matrices with max')
    if (!identical(obj$gene_stat, gene_stat)) stop('CAAS gene_stat mismatch; rebuild the null')
    mats <- obj$caas_corStat_byrank
    if (is.null(mats)) stop('CAAS null lacks caas_corStat_byrank')
    n <- if (length(mats)) ncol(mats[[1]]) else 0L
    for (start in if(n > 0) seq.int(1L, n, by = block_size) else integer()) {
      j <- start:min(n, start + block_size - 1L)
      part <- lapply(mats, function(m) as.matrix(m[, j, drop = FALSE]))
      saveRDS(list(stat = part, gene_stat = gene_stat), sprintf('blocks/blk_%06d.rds', start))
    }
  } else {
    m <- as.matrix(obj$corStat)
    if (!nrow(m) || !ncol(m)) stop('RER null lacks corStat')
    for (start in seq.int(1L, ncol(m), by = block_size)) {
      j <- start:min(ncol(m), start + block_size - 1L)
      x <- m[, j, drop = FALSE]; colnames(x) <- paste0('rer_', j)
      saveRDS(list(corStat = x), sprintf('blocks/blk_%06d.rds', start))
    }
  }
} else if (length(paths)) {
  if (source_type == 'caas') caas_block_index(paths, gene_stat=gene_stat)
  for (i in seq_along(paths)) {
    if (!file.symlink(normalizePath(paths[i]), sprintf('blocks/blk_%06d.rds', i))) stop('Cannot stage block')
  }
}
if (source_type == 'rer' && length(list.files('blocks'))) {
  seen <- character(); settings <- NULL
  for (f in sort(list.files('blocks',full.names=TRUE))) {
    b <- readRDS(f); m <- as.matrix(b$corStat); ids <- colnames(m)
    if (is.null(rownames(m)) || anyNA(rownames(m)) || anyDuplicated(rownames(m))) stop('Invalid RER block gene IDs')
    if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) || any(ids %in% seen)) stop('Invalid or overlapping RER block cycles')
    seen <- c(seen, ids)
    if (!is.null(settings) && !identical(settings,b$rer_settings)) stop('RER block correlation settings differ')
    settings <- b$rer_settings
  }
}
if (!length(list.files('blocks')) && source_type == 'rer') saveRDS(list(corStat=NULL), 'blocks/no_null.rds')
if (!length(list.files('blocks'))) {
  warning('No usable null: CAAS calibrated FCS unavailable; no pathway tests emitted')
  writeLines(paste(c('ranking','database','pathway','stat','pval','p.adj','p.perm','q.perm',
                    'num.genes','gene.vals','sig_perm','direction','redlich_level'), collapse='\t'), 'empty.tsv')
}
