# Verified DOMINO inputs and recurrence display; no recomputation of inference.
ami_module_provenance <- function(path,active_file,network_file) {
 bad <- function(reason)list(available=FALSE,reason=reason)
 if(!file.exists(path))return(bad('missing DOMINO provenance (historical outputs are not certified)'))
 m <- jsonlite::read_json(path,simplifyVector=TRUE)
 if(!identical(as.integer(m$contract_version),2L)||!identical(m$p_value_scale,'raw_upper_tail'))return(bad('old/unknown module p-value scale'))
 genes <- sort(unique(trimws(readLines(active_file,warn=FALSE))));genes<-genes[nzchar(genes)]
 stamp<-paste0(paste(genes,collapse='\n'),if(length(genes))'\n' else '')
 if(!identical(digest::digest(stamp,algo='sha256',serialize=FALSE),m$active_genes_sha256))return(bad('active gene list differs from module inputs'))
 if(is.null(network_file)||!file.exists(network_file)||!identical(digest::digest(file=network_file,algo='sha256'),m$network_sha256))return(bad('network differs from module inputs'))
 list(available=TRUE,reason='verified',metadata=m)
}

ami_recurrence <- function(terms,registry,view) {
 stopifnot(view %in% c('Global','Top','Bottom'))
 plans<-c(caas=tolower(view),fade=paste0('fade_',tolower(view),'_significant'),rer='rer_significant')
 availability<-do.call(rbind,lapply(names(plans),function(tool){
  rs<-registry[[tool]];ix<-which(vapply(rs,function(x)identical(x$list_name,unname(plans[tool])),logical(1)))
  rr<-if(length(ix))rs[[ix[1]]] else NULL
  data.frame(tool=toupper(tool),gene_list=unname(plans[tool]),
   list_label=if(tool=='rer' && view!='Global')'RER: Global (both signs; reference)' else paste(toupper(tool),view,sep=': '),
   available=!is.null(rr)&&isTRUE(rr$available),
   status=if(is.null(rr))'Unavailable: list/method absent' else if(isTRUE(rr$available))'Verified modules and descriptions' else rr$skip_message %||% 'Descriptions unavailable',
   stringsAsFactors=FALSE)
 }))
 required<-c('tool','gene_list','term')
 if(!all(required %in% names(terms)))terms<-data.frame(tool=character(),gene_list=character(),term=character())
 terms<-terms[!is.na(terms$term)&nzchar(terms$term),,drop=FALSE]
 pres<-unique(merge(terms[required],availability[availability$available,c('tool','gene_list')],by=c('tool','gene_list')))
 labels<-sort(unique(terms$term))
 cells<-data.frame(term=character(),tool=character(),gene_list=character(),list_label=character(),present=integer())
 if(length(labels))cells<-do.call(rbind,lapply(labels,function(term){
  d<-availability
  d$term<-term;d$present<-ifelse(d$available,as.integer(paste(d$tool,d$gene_list,sep='\t') %in% paste(pres$tool[pres$term==term],pres$gene_list[pres$term==term],sep='\t')),NA_integer_)
  d[c('term','tool','gene_list','list_label','present')]
 }))
 summary<-data.frame(term=labels,stringsAsFactors=FALSE)
 for(tool in c('CAAS','FADE','RER')) {
  z<-availability[availability$tool==tool,]
  summary[[tool]]<-if(!z$available)rep(NA_integer_,length(labels)) else as.integer(labels %in% pres$term[pres$tool==tool])
 }
 summary$n_evaluable_lists<-rep(sum(availability$available),length(labels))
 summary$n_lists_with_term<-if(length(labels))rowSums(summary[c('CAAS','FADE','RER')],na.rm=TRUE) else integer()
 summary$fraction_evaluable<-ifelse(summary$n_evaluable_lists>0,summary$n_lists_with_term/summary$n_evaluable_lists,NA_real_)
 list(cells=cells,summary=summary,availability=availability)
}

# Legacy filenames are preserved for consumers; visible labels describe Rho.
ami_list_label <- function(x) switch(x,rer_accelerating="RER positive",rer_decelerating="RER negative",rer_significant="RER global (both signs)",x)
