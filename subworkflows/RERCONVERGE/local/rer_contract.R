# Shared input, provenance and tabular contracts for RERconverge.
rer_ids <- function(x, context) {
  if (is.null(x) || anyNA(x) || any(!nzchar(x)) || any(grepl('[\t\r\n]', x)) || anyDuplicated(x))
    stop(context, ': missing, duplicate or invalid identifiers')
  invisible(x)
}

rer_table <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop('Missing table: ', path)
  first <- readLines(path, n=1, warn=FALSE)
  if (!length(first)) stop('Empty table: ', path)
  sep <- if (grepl('\t', first, fixed=TRUE)) '\t' else ','
  # TSV identifiers are literal, including quotes, leading zeros and NA.
  x <- read.table(path, header=TRUE, sep=sep, colClasses='character',
                  na.strings=NULL, quote=if (sep == ',') '"' else '',
                  comment.char='', check.names=FALSE, fill=FALSE)
  rer_ids(names(x), 'Table columns')
  x
}

rer_number <- function(x, context, missing=FALSE) {
  absent <- is.na(x) | trimws(x) %in% c('', 'NA', 'NaN')
  y <- suppressWarnings(as.numeric(x))
  if (any(!absent & !is.finite(y)) || (!missing && any(absent))) stop(context, ': expected finite numbers')
  y[absent] <- NA_real_
  y
}

rer_integer <- function(x, context, minimum=0) {
  y <- rer_number(x, context)
  if (length(y) != 1L || y < minimum || y != floor(y) || y > .Machine$integer.max)
    stop(context, ': expected integer >= ', minimum)
  as.integer(y)
}

rer_species <- function(x) gsub('[[:space:]]+', '_', trimws(x))

rer_read_traits <- function(path, species_col, trait_col=NULL) {
  x <- rer_table(path)
  if (!species_col %in% names(x)) stop('Missing species column: ', species_col)
  x[[species_col]] <- rer_species(x[[species_col]])
  rer_ids(x[[species_col]], 'Trait species')
  if (!is.null(trait_col) && nzchar(trait_col)) {
    if (!trait_col %in% names(x)) stop('Missing trait column: ', trait_col)
    x[[trait_col]] <- rer_number(x[[trait_col]], 'Trait values', missing=TRUE)
    x <- x[!is.na(x[[trait_col]]), , drop=FALSE]
  }
  x
}

rer_trait <- function(path) {
  e <- new.env(parent=emptyenv())
  load(path, envir=e)
  if (!exists('trait_vector', e, inherits=FALSE)) stop('Missing trait_vector')
  v <- e$trait_vector
  rer_ids(names(v), 'Trait species')
  if (!is.numeric(v) || any(!is.finite(v)) || length(unique(v)) < 2L)
    stop('Trait must contain finite, nonconstant numeric values')
  list(values=v, n=e$n_vector, c=e$c_vector, metadata=e$rer_trait_metadata)
}

rer_hash <- function(x) digest::digest(x, algo='sha256', serialize=TRUE)
rer_file_hash <- function(path) digest::digest(file=path, algo='sha256')

rer_verify_snapshot <- function(path) {
  folder <- dirname(path)
  if (basename(folder) != 'rer_inputs') return(invisible(NULL))
  manifest <- file.path(folder,'checksums.json')
  if (!file.exists(manifest)) stop('Missing RER snapshot checksums')
  hashes <- jsonlite::read_json(manifest,simplifyVector=TRUE)
  required <- c('validated_rer.output','rer_run_metadata.json')
  if (!all(required %in% names(hashes)) || !any(grepl('^rerconverge_summary_.*[.]tsv$',names(hashes))))
    stop('Incomplete RER snapshot manifest')
  if (any(basename(names(hashes)) != names(hashes))) stop('Invalid RER snapshot path')
  for (name in names(hashes)) {
    f <- file.path(folder,name)
    if (!file.exists(f) || !identical(rer_file_hash(f),unname(hashes[[name]])))
      stop('RER snapshot checksum mismatch: ',name)
  }
  invisible(NULL)
}

rer_tree_contract <- function(trees) {
  if (is.null(trees$masterTree) || is.null(trees$paths)) stop('Invalid treesObj')
  rer_ids(names(trees$trees), 'Tree genes')
  rer_ids(trees$masterTree$tip.label, 'Master tree tips')
  nodes <- trees$ap$nodeId
  paths <- if (!is.null(nodes)) paste(nodes[[1]],nodes[[2]],sep='->') else colnames(trees$paths)
  rer_ids(paths, 'Tree paths')
  list(version='rer_v1', tree=rer_hash(list(master=trees$masterTree, paths=trees$paths,
        genes=names(trees$trees))), paths=paths, genes=names(trees$trees))
}

rer_check_matrix <- function(mat, trees, species) {
  expected <- rer_tree_contract(trees)
  stamp <- attr(mat, 'rer_contract')
  if (is.null(stamp) || !identical(stamp$version,'rer_v1')) stop('RER matrix lacks verified provenance; rebuild with build_matrix')
  if (!identical(stamp$tree, expected$tree) || !identical(stamp$paths, colnames(mat)) ||
      !identical(colnames(mat), expected$paths) || !identical(stamp$genes, rownames(mat)) ||
      !identical(rownames(mat), expected$genes) || !identical(sort(species), sort(stamp$species)))
    stop('RER matrix/tree/species identity or path order mismatch')
  rer_ids(rownames(mat), 'Matrix genes')
  if (!is.matrix(mat) || !is.numeric(mat) || any(is.infinite(mat))) stop('Invalid RER matrix')
  invisible(stamp)
}

rer_transform <- function(trait, requested) {
  if (!requested %in% c('auto','ha_logit','logit','arcsin','log10','none')) stop('Unknown trait transformation: ', requested)
  v <- trait$values
  effective <- requested
  shift <- NULL
  if (requested == 'auto') {
    if (all(v >= 0 & v <= 1) && !all(v %in% c(0,1))) effective <- 'logit'
    else {
      if (length(v) < 3L || length(v) > 5000L) stop('auto Shapiro-Wilk requires 3..5000 observations; choose an explicit transform')
      effective <- if (stats::shapiro.test(v)$p.value < .05) 'log10' else 'none'
    }
  }
  if (effective == 'ha_logit') {
    if (is.null(trait$n) || is.null(trait$c)) stop('ha_logit requires n_trait and c_trait counts')
    for (a in c('n','c')) {
      rer_ids(names(trait[[a]]), paste('Count',a,'species'))
      if (!setequal(names(v), names(trait[[a]]))) stop('Counts and trait species differ')
    }
    n <- trait$n[names(v)]; c <- trait$c[names(v)]
    if (any(!is.finite(n) | !is.finite(c) | n < 0 | c < 0 | c > n | n != floor(n) | c != floor(c)))
      stop('Counts require finite integers with 0 <= c <= n')
    v <- log((c+.5)/(n-c+.5))
  } else if (effective == 'logit') {
    if (any(v < 0 | v > 1)) stop('logit requires proportions in [0,1]')
    v <- stats::qlogis(pmax(pmin(v,1-1e-4),1e-4))
  } else if (effective == 'arcsin') {
    if (any(v < 0 | v > 1)) stop('arcsin requires proportions in [0,1]')
    v <- asin(sqrt(v))
  } else if (effective == 'log10') {
    if (any(v < 0) || !any(v > 0)) stop('log10 requires nonnegative values and a positive value')
    shift <- min(v[v > 0])/10
    v <- log10(v+shift)
  }
  if (any(!is.finite(v)) || length(unique(v)) < 2L) stop('Transformed trait is nonfinite or constant')
  list(values=v, requested=requested, effective=effective, shift=shift)
}

rer_results <- function(res) {
  if (!is.data.frame(res) || !all(c('Rho','N','P','p.adj') %in% names(res))) stop('Invalid RER result schema')
  rer_ids(rownames(res), 'Result genes')
  expected <- attr(res,'rer_settings')$matrix_contract$genes
  if (!is.null(expected) && !identical(rownames(res),expected)) stop('Result gene coverage differs from matrix roster')
  for (n in intersect(c('Rho','N','P','p.adj','p.perm','p.perm.adj'),names(res))) {
    if (is.logical(res[[n]]) && all(is.na(res[[n]]))) res[[n]] <- as.numeric(res[[n]])
    if (!is.numeric(res[[n]]) || any(is.infinite(res[[n]]))) stop('Invalid result column: ', n)
  }
  for (n in intersect(c('P','p.adj','p.perm','p.perm.adj'), names(res))) {
    if (any(res[[n]] < 0 | res[[n]] > 1, na.rm=TRUE)) stop('Invalid p-value: ', n)
  }
  if (any(abs(res$Rho) > 1, na.rm=TRUE) || any(res$N < 0 | res$N != floor(res$N),na.rm=TRUE)) stop('Invalid Rho or N')
  if ('p.perm' %in% names(res) && any(!is.na(res$p.perm) & (is.na(res$Rho) | is.na(res$P) | res$p.perm == 0))) stop('Empirical evidence without valid association')
  if (!'p.perm' %in% names(res)) res$p.perm <- rep(NA_real_,nrow(res))
  bh <- stats::p.adjust(res$p.perm, 'BH')
  if ('p.perm.adj' %in% names(res) && !isTRUE(all.equal(res$p.perm.adj,bh,check.attributes=FALSE))) stop('Empirical BH disagrees with full result family')
  res$p.perm.adj <- bh
  res
}

rer_summary <- function(res, nominal, fdr) {
  res <- rer_results(res)
  x <- data.frame(gene=rownames(res),res, row.names=NULL,check.names=FALSE)
  x$signif_param <- !is.na(x$p.adj) & x$p.adj < fdr
  x$signif <- !is.na(x$p.perm) & x$p.perm < nominal
  x$signif_perm_fdr <- !is.na(x$p.perm.adj) & x$p.perm.adj < fdr
  x$pos_corr <- x$signif & x$Rho > 0
  x$neg_corr <- x$signif & x$Rho < 0
  x[['Empirical available']] <- !is.na(x$p.perm)
  x[['BH empirical']] <- ifelse(is.na(x$p.perm),NA,x$signif_perm_fdr)
  x$Call <- ifelse(x$pos_corr,'Positively correlated',ifelse(x$neg_corr,'Negatively correlated',ifelse(x$signif,'Significant','')))
  x$rer_stat <- sign(x$Rho)*-log10(x$P)
  x
}

rer_write_tsv <- function(x,path) {
  old <- options(digits=17); on.exit(options(old))
  write.table(x,path,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
}
