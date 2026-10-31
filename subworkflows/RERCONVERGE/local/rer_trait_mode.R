a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
t <- rer_trait(a[1]); mode <- a[2]
if (!mode %in% c('auto','continuous','binary')) stop('Invalid rer_trait_mode')
recorded <- t$metadata$mode
if (mode == 'auto') mode <- if(!is.null(recorded)) recorded else if(length(unique(t$values))==2L) 'binary' else 'continuous'
if (!is.null(recorded) && mode != recorded) stop('Stored trait mode conflicts with requested mode; rebuild trait from raw measurements')
if (mode == 'binary' && !setequal(unique(t$values),c(0,1))) stop('Precomputed binary trait must contain 0/1; rebuild trait')
writeLines(mode,a[3])
