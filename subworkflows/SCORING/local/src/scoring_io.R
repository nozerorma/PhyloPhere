# Literal TSV identifiers are never inferred or interpreted as missing values.
scoring_read <- function(path, quote = '', delim = '\t') {
  x <- readr::read_delim(path, delim=delim, quote=quote, col_types=readr::cols(.default=readr::col_character()),
                        na=character(), trim_ws=FALSE, name_repair='minimal', show_col_types=FALSE, progress=FALSE)
  if (nrow(readr::problems(x)) || anyDuplicated(names(x)) || any(!nzchar(names(x)))) stop('Malformed scoring table: ',path)
  ids <- intersect(c('Gene','gene','side','cycle','ranking','family_id'),names(x))
  other <- setdiff(names(x),ids)
  for(n in other) x[[n]] <- type.convert(x[[n]],as.is=TRUE,na.strings=c('','NA'),numerals='allow.loss')
  x
}
scoring_keys <- function(x, keys, label) {
  if(!all(keys %in% names(x))) stop('Missing ',label,' key columns')
  for(k in keys) if(anyNA(x[[k]]) || any(!nzchar(as.character(x[[k]])))) stop('Empty ',label,' key')
  if(anyDuplicated(x[keys])) stop('Duplicate ',label,' keys')
  invisible(x)
}
# Match producer species-labelled IDs to source IDs, preserving exact matches and
# rejecting collisions. Numbered IDs without a species label remain literal.
scoring_gene_aliases <- function(gene) {
  stem <- sub('\\.prot$', '', gene)
  core <- sub('\\.[A-Z][a-z]+_[a-z]+(_[A-Za-z0-9]+)*$', '', stem, perl=TRUE)
  unique(c(gene, stem, core, if (core != stem) sub('\\.[0-9]+$', '', core)))
}
scoring_match_genes <- function(targets, sources) {
  version <- function(gene) {
    stem <- sub('\\.prot$', '', gene)
    stem <- sub('\\.[A-Z][a-z]+_[a-z]+(_[A-Za-z0-9]+)*$', '', stem, perl=TRUE)
    if (grepl('\\.[0-9]+$', stem)) sub('^.*\\.([0-9]+)$', '\\1', stem) else NA_character_
  }
  index <- new.env(hash=TRUE, parent=emptyenv())
  for (i in seq_along(sources)) for (a in scoring_gene_aliases(sources[[i]]))
    index[[a]] <- c(index[[a]], i)
  exact <- match(targets, sources)
  for (i in which(is.na(exact))) {
    hits <- unique(unlist(lapply(scoring_gene_aliases(targets[[i]]), function(a) index[[a]]), use.names=FALSE))
    target_version <- version(targets[[i]])
    if (!is.na(target_version) && length(hits)) {
      source_versions <- vapply(sources[hits], version, character(1))
      hits <- hits[is.na(source_versions) | source_versions == target_version]
    }
    if (length(hits) > 1L) stop('Ambiguous producer gene alias: ', targets[[i]], ': ', paste(sources[hits], collapse=', '))
    if (length(hits)) exact[[i]] <- hits
  }
  used <- exact[!is.na(exact)]
  if (anyDuplicated(used)) {
    hit <- used[duplicated(used)][1]
    stop('Ambiguous producer gene alias: ', sources[[hit]], ' matches ',
         paste(targets[!is.na(exact) & exact == hit], collapse=', '))
  }
  exact
}
scoring_number <- function(x, label, lower=-Inf, upper=Inf, missing=FALSE, integer=FALSE, positive_inf=FALSE) {
  if(is.logical(x) && any(!is.na(x))) stop('Invalid ',label,': logical value')
  if(is.numeric(x) && any(is.nan(x))) stop('Invalid ',label,': NaN')
  v <- suppressWarnings(as.numeric(x))
  absent <- is.na(x) | x %in% c('','NA')
  if(any((!absent & (is.na(v) | !(is.finite(v) | (positive_inf & v==Inf)))) | (absent & !missing)) ||
     any(v < lower | v > upper | (integer & v != floor(v)),na.rm=TRUE)) stop('Invalid ',label)
  v[absent] <- NA_real_
  v
}
scoring_genes <- function(path) {
  x <- readLines(path,warn=FALSE)
  if(length(x) && x[1] %in% c('Gene','gene')) x <- x[-1]
  if(anyNA(x) || any(!nzchar(x)) || any(trimws(x)!=x) || anyDuplicated(x)) stop('Invalid/duplicate universe genes')
  x
}
scoring_write <- function(x,path,...) {
  # Round-trip doubles even across values immediately adjacent to a call threshold.
  for(n in names(x)) if(is.numeric(x[[n]])) x[[n]] <- ifelse(is.na(x[[n]]),'NA',sprintf('%.17g',x[[n]]))
  readr::write_delim(x,path,delim='\t',quote='none',escape='none',na='NA',...)
}
