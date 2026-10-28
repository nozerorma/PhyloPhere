# Shared presentation helpers. Raw tables and inferential families are unchanged.
report_db_dictionary <- function() data.frame(
 database=c('geneontology_Biological_Process_noRedundant','geneontology_Molecular_Function_noRedundant','geneontology_Cellular_Component_noRedundant','pathway_KEGG','pathway_Reactome','pathway_Wikipathway','pathway_Wikipathway_cancer','network_CORUMA','network_Transcription_Factor_target','network_Kinase_target','phenotype_Human_Phenotype_Ontology','disease_Disgenet','disease_OMIM','TCGA.PanCancer.OncoVar','h','c6'),
 abbreviation=c('GO-BP','GO-MF','GO-CC','KEGG','REAC','WP','WP-C','CORUM','TF','KIN','HPO','DisGeNET','OMIM','OncoVar','HALLMARK','C6'),stringsAsFactors=FALSE)
report_db_abbrev <- function(x) {
 d <- report_db_dictionary(); a <- d$abbreviation[match(x,d$database)]
 ifelse(is.na(a),as.character(x),a)
}
report_collapsible <- function(x, limit=85L) vapply(as.character(x),function(s) {
 if(is.na(s)) return(NA_character_)
 esc <- function(z) as.character(htmltools::htmlEscape(z))
 if(nchar(s)<=limit) return(esc(s))
 paste0('<details><summary>',esc(paste0(substr(s,1,limit),'…')),'</summary><div style="max-width:65em;white-space:normal">',esc(s),'</div></details>')
},character(1),USE.NAMES=FALSE)
report_table <- function(df, caption=NULL, sci_cols=NULL, ...) {
 df <- as.data.frame(df,check.names=FALSE)
 if('database' %in% names(df) && !'abbreviation' %in% names(df)) df$database <- report_db_abbrev(df$database)
 # Provenance is enforced and exported by producers, rather than shown as an annotation.
 df <- df[,!names(df) %in% c('gene_stat','flag_gene_stat'),drop=FALSE]
 for(n in names(df)) if(is.factor(df[[n]])) df[[n]] <- as.character(df[[n]])
 text_cols <- unname(which(vapply(df,is.character,logical(1)))-1L)
 render_text <- DT::JS("function(data,type,row,meta){if(type !== 'display' || typeof data !== 'string') return data; function esc(s){return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;');} if(data.length <= 85) return esc(data); return '<details><summary>'+esc(data.slice(0,85))+'…</summary><div style=\"max-width:65em;white-space:normal\">'+esc(data)+'</div></details>';}" )
 w <- DT::datatable(df,caption=caption,rownames=FALSE,filter='top',escape=FALSE,extensions='Buttons',
   options=list(pageLength=15,lengthMenu=c(10,15,25,50,100),scrollX=TRUE,autoWidth=FALSE,deferRender=TRUE,processing=TRUE,order=list(),dom='Blfrtip',
     buttons=list(list(extend='csv',exportOptions=list(orthogonal='export')),list(extend='excel',exportOptions=list(orthogonal='export'))),
     columnDefs=if(length(text_cols))list(list(targets=text_cols,render=render_text)) else list()))
 nums <- names(df)[vapply(df,is.numeric,logical(1)) & !vapply(df,is.integer,logical(1))]
 # Counts remain integers; continuous evidence always has four decimal places.
 nums <- nums[!grepl('^(n_|N |Rows$|checks$|num.genes$|Genes$|Module$|Position$|.*count$)',nums,ignore.case=TRUE)]
 # sci_cols (p-values and similar quantities that span many orders of magnitude) keep three significant digits.
 sci_cols <- intersect(sci_cols,names(df)); nums <- setdiff(nums,sci_cols)
 if(length(nums)) w <- DT::formatRound(w,nums,digits=4,mark='')
 if(length(sci_cols)) w <- DT::formatSignif(w,sci_cols,digits=3,mark='')
 w
}
# Inert templates keep hidden/remote sections from initializing large tables.
# Each widget is bound by HTMLWidgets when its section reaches the viewport.
.report_lazy_script_emitted <- FALSE
report_emit <- function(w) {
 r<-knitr::knit_print(w);knitr::knit_meta_add(attr(r,'knit_meta'))
 if(!.report_lazy_script_emitted) {
  cat('<script>(function(){function start(){var pending=document.querySelectorAll(".report-widget-pending");function mount(el){if(el.dataset.mounted)return;el.dataset.mounted="1";var t=el.querySelector("template");if(!t)return;el.replaceChildren(t.content);if(window.HTMLWidgets)HTMLWidgets.staticRender();}if(window.IntersectionObserver){var observer=new IntersectionObserver(function(items){items.forEach(function(item){if(item.isIntersecting){mount(item.target);observer.unobserve(item.target);}});},{rootMargin:"0px 0px 200px 0px"});pending.forEach(function(el){observer.observe(el);});}else{pending.forEach(mount);}}if(document.readyState==="loading")document.addEventListener("DOMContentLoaded",start);else start();})();</script>\n')
  .report_lazy_script_emitted <<- TRUE
 }
 cat('<div class="report-widget-pending" style="min-height:100px"><template>',r,'</template></div>\n\n',sep='')
 invisible(NULL)
}

report_representatives <- function(d) {
 if(!'is_representative' %in% names(d)) return(d)
 # Ungrouped rows have no defensible representative and remain visible.
 d[is.na(d$is_representative) | d$is_representative %in% TRUE,,drop=FALSE]
}
report_term_views <- function(d, level=3L, signed=FALSE) {
 if(signed) {
  for(s in c('positive','negative','zero','unavailable')) {
   v<-if('direction' %in% names(d)) as.character(d$direction) else rep(NA_character_,nrow(d))
   sub<-d[if(s=='unavailable')is.na(v) else !is.na(v)&v==s,,drop=FALSE]
   if(!nrow(sub)) next
   cat('\n',strrep('#',level),' ',tools::toTitleCase(s),' {.tabset}\n\n',sep='')
   report_term_views(sub,level+1L,FALSE)
  }
 } else {
  cat('\n',strrep('#',level),' Representatives\n\n',sep='')
  report_emit(report_table(report_representatives(d)))
  cat('\n',strrep('#',level),' All terms\n\n',sep='')
  report_emit(report_table(d))
 }
 invisible(NULL)
}
report_term_reduction_note <- function() cat('Representatives summarize blocks of terms with overlapping genes; All terms preserves every member. FCS blocks group the selected sets of a ranking and effect sign by complete linkage on the Jaccard similarity of their full gene sets (saved threshold). AMI blocks group the terms of one module by complete linkage on the Jaccard similarity of the module genes they annotate (saved threshold); the module label and the block summary show the first blocks only (saved maximum), and the tables keep every block. Grouping never changes p/q values or significance families.\n\n')
