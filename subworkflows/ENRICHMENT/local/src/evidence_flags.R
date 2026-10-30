# Report annotations are read from calibrated evidence; no score-derived selection.
evidence_truth <- function(x) !is.na(x) & x %in% c(TRUE,'TRUE','True','true',1,'1')
evidence_column <- function(df,name,default=NA_real_) if(name %in% names(df)) df[[name]] else rep(default,nrow(df))
# CAAS nominal/BH flags of the stats table. The calibrated columns (caas_nominal, caas_adjusted and their _top/_bottom
# versions) carry the thresholds of the calibration. Only when a flag column is absent is it derived from its p column,
# and then only with the thresholds given (the nominal and fdr_threshold of the calibration metadata); without them it is NA.
fcs_evidence_flags <- function(stats, nominal = NULL, fdr = NULL) {
 out<-data.frame(gene=as.character(stats$gene),stringsAsFactors=FALSE)
 for(suffix in c('','_top','_bottom'))for(kind in c('nominal','adjusted')) {
  n<-paste0('caas_',kind,suffix)
  if(n %in% names(stats)) {
   v<-stats[[n]];out[[n]]<-ifelse(is.na(v),NA,evidence_truth(v));next
  }
  pn<-paste0('caas_',if(kind=='nominal')'p.gen' else 'p.gen.adj',suffix)
  thr<-if(kind=='nominal')nominal else fdr
  out[[n]]<-if(pn %in% names(stats) && !is.null(thr)) {p<-suppressWarnings(as.numeric(stats[[pn]]));ifelse(is.finite(p),p<thr,NA)} else NA
 }
 out
}
