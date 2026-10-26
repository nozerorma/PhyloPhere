# Reporting only: preserve producer p/q and their original adjustment families.
scoring_overlap <- function(d, nominal=.05, fdr=.05) {
  if(anyDuplicated(d$Gene)) stop('Duplicate gene evidence keys')
  col <- function(n) if(n %in% names(d)) d[[n]] else rep(NA_real_, nrow(d))
  rows <- list()
  for(rk in c('global','top','bottom')) {
    sf <- if(rk=='global') '' else paste0('_',rk)
    for(tool in c('RER','FADE','Accumulation')) for(level in c('nominal','adjusted')) {
      cp <- col(paste0('caas_',if(level=='nominal') 'p.gen' else 'p.gen.adj',sf))
      ca <- cp < if(level=='nominal') nominal else fdr
      if(tool=='RER') {
        op <- col(if(level=='nominal') 'rer_min_pval' else 'rer_perm_padj')
        measured <- is.finite(op)
        other <- op < if(level=='nominal') nominal else fdr
      } else if(tool=='FADE') {
        sides <- if(rk=='global') c('top','bottom') else rk
        bf <- lapply(sides,function(x) col(paste0('fade_max_bf_',x)))
        calls <- lapply(sides,function(x) col(paste0('fade_significant_',x)))
        # Global FADE is the union of measured sides. A side with no source data
        # stays unavailable; it does not erase evidence from the other side.
        side_measured <- Map(function(b,c) is.finite(b) & !is.na(c), bf, calls)
        other <- Reduce(`|`,lapply(calls,function(x)x %in% TRUE))
        measured <- other | Reduce(`&`, side_measured)
      } else {
        # Only ActualCount>0 genes belong to the accumulation tested BH family.
        # q availability is the producer's marker of that membership.
        measured <- is.finite(col(paste0('accum_fdr',sf))) & !is.na(col(paste0('accum_significant',sf)))
        other <- col(paste0('accum_significant',sf)) %in% TRUE
      }
      ok <- is.finite(cp) & measured
      ca <- ca[ok]; ot <- other[ok]
      cells <- c(sum(ca & ot),sum(ca & !ot),sum(!ca & ot),sum(!ca & !ot))
      pv <- if(length(ca)) stats::fisher.test(matrix(cells,2,byrow=TRUE))$p.value else NA_real_
      rows[[length(rows)+1L]] <- data.frame(ranking=rk,tool=tool,call=level,
        shared_tested_universe=length(ca),CAAS=sum(ca),other=sum(ot),overlap=cells[1],
        both=cells[1],caas_only=cells[2],other_only=cells[3],neither=cells[4],fisher_p=pv)
    }
  }
  out <- do.call(rbind,rows)
  out$fisher_BH <- p.adjust(out$fisher_p,'BH',n=nrow(out))
  out$fisher_family_n <- nrow(out)
  out
}
