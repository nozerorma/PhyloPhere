# Report annotations are read from calibrated evidence; no score-derived selection.
evidence_truth <- function(x) !is.na(x) & x %in% c(TRUE,'TRUE','True','true',1,'1')
evidence_column <- function(df,name,default=NA_real_) if(name %in% names(df)) df[[name]] else rep(default,nrow(df))
fcs_evidence_flags <- function(stats) {
 out<-data.frame(gene=as.character(stats$gene),stringsAsFactors=FALSE)
 for(suffix in c('','_top','_bottom'))for(kind in c('nominal','adjusted')) {
  n<-paste0('caas_',kind,suffix);out[[n]]<-evidence_truth(evidence_column(stats,n,FALSE))
 }
 out
}
