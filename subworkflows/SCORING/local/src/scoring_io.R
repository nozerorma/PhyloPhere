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
