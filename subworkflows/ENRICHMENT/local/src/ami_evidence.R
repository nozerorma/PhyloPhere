# Verified DOMINO inputs and recurrence display; no recomputation of inference.
ami_module_provenance <- function(path,active_file,network_file,background_file=NULL) {
 bad <- function(reason)list(available=FALSE,reason=reason)
 if(!file.exists(path))return(bad('missing DOMINO provenance (historical outputs are not certified)'))
 m <- tryCatch(jsonlite::read_json(path,simplifyVector=TRUE),error=function(e)NULL)
 required <- c('contract_version','p_value_scale','p_adjustment','n_putative','final_modules',
               'active_genes_sha256','network_sha256','slices_sha256','producer_sha256',
               'modules_sha256','stats_sha256','slice_threshold','module_threshold')
 if(is.null(m) || !all(required %in% names(m)) ||
    any(vapply(m[required],length,integer(1))!=1L)) return(bad('incomplete or unreadable DOMINO certificate'))
 if(!identical(as.integer(m$contract_version),2L)||!identical(m$p_value_scale,'raw_upper_tail'))return(bad('old/unknown module p-value scale'))
 if(!identical(m$p_adjustment,'Bonferroni_n_putative'))return(bad('unknown DOMINO adjustment family'))
 counts <- suppressWarnings(as.numeric(c(m$n_putative,m$final_modules)))
 thresholds <- suppressWarnings(as.numeric(c(m$slice_threshold,m$module_threshold)))
 if(any(!is.finite(counts)) || any(counts<0 | counts!=floor(counts)) ||
    any(!is.finite(thresholds)) || any(thresholds<=0 | thresholds>1)) return(bad('invalid DOMINO family counts or thresholds'))
 genes <- sort(unique(trimws(readLines(active_file,warn=FALSE))),method='radix');genes<-genes[nzchar(genes)]
 stamp<-paste0(paste(genes,collapse='\n'),if(length(genes))'\n' else '')
 if(!identical(digest::digest(stamp,algo='sha256',serialize=FALSE),m$active_genes_sha256))return(bad('active gene list differs from module inputs'))
 if(is.null(network_file)||!file.exists(network_file)||!identical(digest::digest(file=network_file,algo='sha256'),m$network_sha256))return(bad('network differs from module inputs'))
 if(!is.null(m$background_sha256) && !is.null(background_file)) {
  bg <- sort(unique(trimws(readLines(background_file,warn=FALSE))),method='radix');bg<-bg[nzchar(bg)]
  bg_stamp<-paste0(paste(bg,collapse='\n'),if(length(bg))'\n' else '')
  if(!identical(digest::digest(bg_stamp,algo='sha256',serialize=FALSE),m$background_sha256))return(bad('background differs from the one the network was built on'))
 }
 producer <- 'src/run_domino_modules.py'
 if(!file.exists(producer) || is.null(m$producer_sha256) || !identical(digest::digest(file=producer,algo='sha256'),m$producer_sha256)) return(bad('DOMINO producer differs or source certificate is absent; regenerate modules'))
 slices <- file.path(dirname(path),'inputs','slices.txt')
 if(!file.exists(slices) || !identical(digest::digest(file=slices,algo='sha256'),m$slices_sha256)) return(bad('DOMINO slices certificate differs or is absent'))
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
   modules_verified=!is.null(rr)&&isTRUE(rr$modules_verified),
   n_modules=if(is.null(rr)||is.null(rr$n_modules))NA_integer_ else rr$n_modules,
   descriptions_available=!is.null(rr)&&isTRUE(rr$descriptions_available),
   available=!is.null(rr)&&isTRUE(rr$available),
   status=if(is.null(rr))'Unavailable: list/method absent' else if(isTRUE(rr$available)) {
     if(identical(rr$n_modules,0L))'Verified: zero modules' else 'Verified modules and descriptions'
   } else rr$skip_message %||% 'Module descriptions unavailable',
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
