# Report annotations are read from calibrated evidence; no score-derived selection.
evidence_truth <- function(x) !is.na(x) & x %in% c(TRUE,'TRUE','True','true',1,'1')
evidence_column <- function(df,name,default=NA_real_) if(name %in% names(df)) df[[name]] else rep(default,nrow(df))
# CAAS nominal flags of the stats table, and one BH flag per level. The calibrated column caas_nominal[_top|_bottom] carries the
# threshold of the calibration; only when it is absent is it derived from caas_p.gen[_top|_bottom] with the `nominal` given (NA without it).
# With `levels` (the BH levels of the evidence) one flag per level is added, caas_q<level>[_top|_bottom]: the q of the gene is below the
# level, read from the BH class written by scoring (caas_bh_class[_top|_bottom]); NA where the class is absent.
evidence_bh_below <- function(cls, level) ifelse(is.na(cls), NA, startsWith(cls,'<') & suppressWarnings(as.numeric(sub('^<','',cls))) <= level)
fcs_evidence_flags <- function(stats, nominal = NULL, levels = NULL) {
 out<-data.frame(gene=as.character(stats$gene),stringsAsFactors=FALSE)
 for(suffix in c('','_top','_bottom')) {
  n<-paste0('caas_nominal',suffix)
  if(n %in% names(stats)) {
   v<-stats[[n]];out[[n]]<-ifelse(is.na(v),NA,evidence_truth(v))
  } else {
   pn<-paste0('caas_p.gen',suffix)
   out[[n]]<-if(pn %in% names(stats) && !is.null(nominal)) {p<-suppressWarnings(as.numeric(stats[[pn]]));ifelse(is.finite(p),p<nominal,NA)} else NA
  }
  cls<-as.character(evidence_column(stats,paste0('caas_bh_class',suffix),NA_character_))
  for(l in levels) out[[paste0('caas_q',format(l,scientific=FALSE),suffix)]]<-evidence_bh_below(cls,l)
 }
 out
}
