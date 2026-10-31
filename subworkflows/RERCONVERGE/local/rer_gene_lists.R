a <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'rer_contract.R'))
x <- rer_table(a[1])
if (!all(c('gene','Rho','P','p.perm') %in% names(x))) stop('Invalid RER summary schema')
rer_ids(x$gene,'Summary genes')
for (n in c('Rho','P','p.perm')) x[[n]] <- rer_number(x[[n]],n,missing=TRUE)
if (any(abs(x$Rho)>1,na.rm=TRUE) || any(x$P<0 | x$P>1,na.rm=TRUE) || any(x$p.perm<=0 | x$p.perm>1,na.rm=TRUE)) stop('Invalid summary evidence')
threshold <- rer_number(a[3],'Nominal threshold')
if (threshold<0 || threshold>1) stop('Invalid nominal threshold')
sig <- !is.na(x$p.perm) & x$p.perm < threshold
if (any(sig & is.na(x$Rho))) stop('Significant gene without rho')
writeLines(x$gene,'background.txt')
writeLines(x$gene[sig],'rer_significant.txt')
writeLines(x$gene[sig & !is.na(x$Rho) & x$Rho>0],'rer_accelerating.txt')
writeLines(x$gene[sig & !is.na(x$Rho) & x$Rho<0],'rer_decelerating.txt')
stats <- data.frame(gene=x$gene,score_global=sign(x$Rho)*-log10(x$P),
  flag_rer_acc=sig & x$Rho>0,flag_rer_decc=sig & x$Rho<0)
if (basename(a[2]) != 'NO_FILE') {
  gs <- rer_table(a[2])
  candidates <- which(tolower(names(gs))=='gene')
  if (length(candidates)!=1L) stop('Gene scores require exactly one gene column')
  names(gs)[candidates] <- 'gene'
  rer_ids(gs$gene,'Gene-score genes')
  want <- intersect(c('score_top','score_bottom','flag_fade','flag_fade_top','flag_fade_bottom'),names(gs))
  for (n in want) {
    v <- gs[[n]][match(stats$gene,gs$gene)]
    if (startsWith(n,'score_')) v <- rer_number(v,n,missing=TRUE)
    else {
      if (any(!is.na(v) & !v %in% c('TRUE','FALSE','NA',''))) stop('Invalid flag column: ',n)
      v <- as.logical(v)
    }
    stats[[n]] <- v
  }
}
rer_write_tsv(stats,'fcs_stats.tsv')
