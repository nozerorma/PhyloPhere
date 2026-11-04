# Reporting only: preserve producer p/q and their original adjustment families.
# One row per ranking, tool and call: the nominal call and one call per BH level (q < level, labelled "q<0.05", ...).
# FADE has no BH, so only its nominal call is tested. The Fisher p of a row is raw: the BH calls are nested, so their tests
# are not independent and nothing is corrected across rows.
scoring_overlap <- function(d, nominal=.05, levels=c(.05,.1,.15)) {
  if(anyDuplicated(d$Gene)) stop('Duplicate gene evidence keys')
  call_levels <- c(list(list(name='nominal',cut=nominal)),lapply(scoring_bh_levels(levels),function(l)list(name=paste0('q',scoring_bh_labels(l)),cut=l)))
  col <- function(n) if(n %in% names(d)) d[[n]] else rep(NA_real_, nrow(d))
  rows <- list()
  for(rk in c('global','top','bottom')) {
    sf <- if(rk=='global') '' else paste0('_',rk)
    for(tool in c('RER','FADE')) for(cl in call_levels) {
      is_nominal <- cl$name=='nominal'
      if(tool=='FADE' && !is_nominal) next
      cp <- col(paste0('caas_',if(is_nominal) 'p.gen' else 'p.gen.adj',sf))
      ca <- cp < cl$cut
      if(tool=='RER') {
        op <- col(if(is_nominal) 'rer_min_pval' else 'rer_perm_padj')
        measured <- is.finite(op)
        other <- op < cl$cut
      } else if(tool=='FADE') {
        sides <- if(rk=='global') c('top','bottom') else rk
        bf <- lapply(sides,function(x) col(paste0('fade_max_bf_',x)))
        calls <- lapply(sides,function(x) col(paste0('fade_significant_',x)))
        # Global FADE is the union of measured sides. A side with no source data
        # stays unavailable; it does not erase evidence from the other side.
        side_measured <- Map(function(b,c) !is.na(b) & !is.na(c), bf, calls)
        other <- Reduce(`|`,lapply(calls,function(x)x %in% TRUE))
        measured <- other | Reduce(`&`, side_measured)
      }
      ok <- is.finite(cp) & measured & !is.na(ca) & !is.na(other)
      ca <- ca[ok]; ot <- other[ok]
      cells <- c(sum(ca & ot),sum(ca & !ot),sum(!ca & ot),sum(!ca & !ot))
      pv <- if(length(ca)) stats::fisher.test(matrix(cells,2,byrow=TRUE))$p.value else NA_real_
      rows[[length(rows)+1L]] <- data.frame(ranking=rk,tool=tool,call=cl$name,
        shared_tested_universe=length(ca),CAAS=sum(ca),other=sum(ot),overlap=cells[1],
        both=cells[1],caas_only=cells[2],other_only=cells[3],neither=cells[4],fisher_p=pv)
    }
  }
  do.call(rbind,rows)
}
