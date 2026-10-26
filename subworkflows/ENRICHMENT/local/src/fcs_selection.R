# Presentation selection; never recompute p/q or their original BH families.
fcs_selection <- function(d, fdr_thr=.15, p_perm_thr=.025) {
  keys <- c('ranking','database','pathway')
  if(!all(c(keys,'stat','p.adj','p.perm','q.perm') %in% names(d))) stop('Missing FCS selection evidence')
  if(any(vapply(d[keys],function(x) anyNA(x) || any(!nzchar(as.character(x))),logical(1)))) stop('Missing FCS identity')
  if(anyDuplicated(d[keys])) stop('Duplicate FCS ranking/database/pathway identity')
  d$sig_RER <- is.finite(d$q.perm) & d$q.perm < fdr_thr
  for(i in seq_along(c(.05,.10,.15))) {
    level <- c(.05,.10,.15)[i]
    d[[c('R05','R10','R15')[i]]] <- is.finite(d$p.adj) & is.finite(d$p.perm) & d$p.adj < level & d$p.perm < p_perm_thr
  }
  stopifnot(all(!d$R05 | d$R10),all(!d$R10 | d$R15))
  d$checks <- as.integer(d$sig_RER)+as.integer(d$R05)+as.integer(d$R10)+as.integer(d$R15)
  d$C <- (as.integer(d$sig_RER)+(as.integer(d$R05)+as.integer(d$R10)+as.integer(d$R15))/3)/2
  d$selected_C <- d$C > 0
  d$direction <- ifelse(!is.finite(d$stat),NA_character_,ifelse(d$stat>0,'positive',ifelse(d$stat<0,'negative','zero')))
  d
}

# Shared by tables, dotplots, leading-edge composition and group representatives.
fcs_selection_order <- function(d) {
  d[order(-d$C,-d$checks,d$q.perm,d$p.perm,d$database,d$pathway,na.last=TRUE),,drop=FALSE]
}
