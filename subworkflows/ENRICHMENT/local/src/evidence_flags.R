# Report annotations are read from calibrated evidence; no score-derived selection.
evidence_truth <- function(x) !is.na(x) & x %in% c(TRUE,'TRUE','True','true',1,'1')
evidence_column <- function(df,name,default=NA_real_) if(name %in% names(df)) df[[name]] else rep(default,nrow(df))
fcs_evidence_flags <- function(stats) {
 out<-data.frame(gene=as.character(stats$gene),stringsAsFactors=FALSE)
 for(suffix in c('','_top','_bottom'))for(kind in c('nominal','adjusted')) {
  n<-paste0('caas_',kind,suffix)
  v<-evidence_column(stats,n,NA)
  pn<-paste0('caas_',if(kind=='nominal')'p.gen' else 'p.gen.adj',suffix)
  if(pn %in% names(stats)) {p<-suppressWarnings(as.numeric(stats[[pn]]));out[[n]]<-ifelse(is.finite(p),p<.05,NA)}
  else out[[n]]<-ifelse(is.na(v),NA,evidence_truth(v))
 }
 out
}
