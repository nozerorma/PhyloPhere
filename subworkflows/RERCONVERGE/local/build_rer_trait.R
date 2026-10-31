# Build a validated trait while preserving the original measurements and counts.
a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
if (length(a) < 5L) stop('Expected trait table, species column, trait column, RData and type output')
x <- rer_read_traits(a[1], a[2], a[3])
raw_trait_vector <- setNames(x[[a[3]]], x[[a[2]]])
if (length(unique(raw_trait_vector)) < 2L) stop('Trait must be nonconstant and nonempty')
mode <- if (length(a) >= 9L) a[9] else 'auto'
if (!mode %in% c('auto','binary','continuous')) stop('Unknown rer_trait_mode')
u <- sort(unique(raw_trait_vector))
trait_type <- if (mode == 'auto') { if (length(u)==2L) 'binary' else 'continuous' } else mode
trait_vector <- raw_trait_vector
if (trait_type == 'binary') {
  if (length(u) != 2L) stop('Binary mode requires exactly two finite values')
  trait_vector[] <- as.numeric(trait_vector == u[2])
}
counts <- function(i) {
  col <- if (length(a) >= i) a[i] else ''
  if (!nzchar(col)) return(NULL)
  if (!col %in% names(x)) stop('Missing count column: ', col)
  setNames(rer_number(x[[col]], paste('Counts', col)), x[[a[2]]])
}
n_vector <- counts(6); c_vector <- counts(7)
if (xor(is.null(n_vector), is.null(c_vector))) stop('Specify both n_trait and c_trait')
if (!is.null(n_vector) && any(n_vector < 0 | c_vector < 0 | c_vector > n_vector |
                            n_vector != floor(n_vector) | c_vector != floor(c_vector)))
  stop('Counts require integers with 0 <= c <= n')
rer_trait_metadata <- list(version='rer_v1', requested_mode=mode, mode=trait_type,
  binary_encoding=if(trait_type=='binary') setNames(c(0,1),as.character(u)) else NULL,
  input_sha256=rer_file_hash(a[1]), species_column=a[2], trait_column=a[3])
save(trait_vector, raw_trait_vector, n_vector, c_vector, rer_trait_metadata, file=a[4])
writeLines(trait_type,a[5])
