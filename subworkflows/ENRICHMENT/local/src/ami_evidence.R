# Verified DOMINO inputs and recurrence display; no recomputation of inference.
ami_reference_symbol <- function(genes) {
 stem <- sub('\\.prot$', '', genes)
 sub('\\.[A-Z][a-z]+_[a-z]+(_[A-Za-z0-9]+)*$', '', stem, perl=TRUE)
}

# The scores table also carries copies of RER/FADE evidence on CAAS IDs.
# Each method's own background defines which rows belong to its seed family.
# Scores carry CAAS IDs (GENE.Species) and a method background may use bare symbols
# (FADE) or suffixed IDs (RER), so rows are matched on the reference symbol and the
# result is returned in the background's own identifier form, as the seed lists are.
ami_method_seed_evidence <- function(scores, hit, background) {
 if (length(hit) != nrow(scores)) stop('Invalid method evidence mask')
 hit_symbols <- unique(ami_reference_symbol(scores$Gene[hit %in% TRUE]))
 background[ami_reference_symbol(background) %in% hit_symbols]
}

ami_module_provenance <- function(path,active_file,network_file,background_file=NULL) {
 bad <- function(reason)list(available=FALSE,reason=reason)
 if(!file.exists(path))return(bad('missing DOMINO provenance (historical outputs are not certified)'))
 m <- tryCatch(jsonlite::read_json(path,simplifyVector=TRUE),error=function(e)NULL)
 required <- c('contract_version','p_value_scale','p_adjustment','n_putative','final_modules',
               'active_genes_sha256','slice_threshold','module_threshold')
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
 strict <- Sys.getenv('PHYLOPHERE_VERIFY_CHECKSUMS') == '1'
 if(is.null(network_file)||!file.exists(network_file)||(strict && !identical(digest::digest(file=network_file,algo='sha256'),m$network_sha256)))return(bad('network differs from module inputs'))
 if(!is.null(m$background_sha256) && !is.null(background_file)) {
  bg <- sort(unique(trimws(readLines(background_file,warn=FALSE))),method='radix');bg<-bg[nzchar(bg)]
  bg_stamp<-paste0(paste(bg,collapse='\n'),if(length(bg))'\n' else '')
  if(!identical(digest::digest(bg_stamp,algo='sha256',serialize=FALSE),m$background_sha256))return(bad('background differs from the one the network was built on'))
 }
 producer <- 'src/run_domino_modules.py'
 if(!file.exists(producer) || (strict && !identical(digest::digest(file=producer,algo='sha256'),m$producer_sha256))) return(bad('DOMINO producer differs or source certificate is absent; regenerate modules'))
 slices <- file.path(dirname(path),'inputs','slices.txt')
 if(!file.exists(slices) || (strict && !identical(digest::digest(file=slices,algo='sha256'),m$slices_sha256))) return(bad('DOMINO slices certificate differs or is absent'))
 list(available=TRUE,reason=if(strict) 'verified' else 'scientific metadata checked; file checksums not verified',metadata=m)
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
 terms<-terms[paste(terms$tool,terms$gene_list) %in% paste(availability$tool,availability$gene_list),,drop=FALSE]  # nominal lists only: BH level lists are separate analyses, not part of the cross-method recurrence
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

# Displayed recurrence: the terms found in at least `min_lists` evaluable lists, most recurrent first (radix order, so
# ties do not depend on the locale). A term in a single list says nothing about recurrence; the exports keep the full summary.
ami_recurrent_terms <- function(summary, min_lists = 2L) {
 if (!nrow(summary) || max(summary$n_evaluable_lists) < min_lists) return(summary[0, , drop = FALSE])
 shown <- summary[summary$n_lists_with_term >= min_lists, , drop = FALSE]
 shown[order(-shown$n_lists_with_term, shown$term, method = 'radix'), , drop = FALSE]
}

# One line above the displayed table: what the filter kept, or why recurrence cannot be evaluated.
ami_recurrence_note <- function(summary, min_lists = 2L) {
 if (!nrow(summary)) return('_No verified module descriptions available in this view._\n\n')
 if (max(summary$n_evaluable_lists) < min_lists)
  return(sprintf('_Fewer than %d evaluable lists in this view: recurrence is undefined._\n\n', min_lists))
 sprintf('_%d of %d terms occur in at least %d evaluable lists; the exported tables keep every term._\n\n',
         nrow(ami_recurrent_terms(summary, min_lists)), nrow(summary), min_lists)
}

# BH classes in the AMI report. The class labels ("<0.05", ..., "ns") are written by scoring; `levels` are the BH levels.
# The colour ladder is that of report_bh_colours() in TRAIT_ANALYSIS/local/src/report_style.R, which the AMI report does not read.
ami_bh_token <- function(levels) vapply(levels, format, '', scientific = FALSE)
ami_bh_labels <- function(levels) paste0('<', ami_bh_token(sort(levels)))   # strictest first
ami_bh_colours <- function(labels) {
 ladder <- if (length(labels) == 1L) '#E15759' else grDevices::colorRampPalette(c('#4E79A7', '#F28E2B', '#E15759'))(length(labels))
 c(ns = '#BAB0AC', stats::setNames(rev(ladder), labels))
}
# TRUE where the q of the class is below `level` (the class is the strictest level below q, so also below every looser one).
ami_bh_below <- function(cls, level) ifelse(is.na(cls), NA, startsWith(cls, '<') & suppressWarnings(as.numeric(sub('^<', '', cls))) <= level)

# Level lists of the tabs: "global_q0.05" is the BH level 0.05 of the list "global"; any other name is the nominal list.
ami_level_value <- function(name) ifelse(grepl('_q[0-9.]+$', name), suppressWarnings(as.numeric(sub('^.*_q', '', name))), 0)
ami_level_title <- function(name) ifelse(grepl('_q[0-9.]+$', name), paste0('q<', sub('^.*_q', '', name)), 'Nominal')

# Legacy filenames are preserved for consumers; visible labels describe Rho.
ami_list_label <- function(x) switch(x,rer_accelerating="RER positive",rer_decelerating="RER negative",rer_significant="RER global (both signs)",x)
