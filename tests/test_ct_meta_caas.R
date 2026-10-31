# Actual report renders and scientific/reporting invariants on synthetic inputs.
root <- normalizePath('.')
work <- tempfile('ct-meta-caas-'); dir.create(work)
cat('ARTIFACT_DIR:', work, '\n')
source(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'))
set_report_style()
stopifnot(length(report_palette(0))==0, length(unique(report_palette(24)))==24,
          identical(report_theme()$text$family, 'sans'))
plot <- ggplot2::ggplot(data.frame(x=c('a','b'),y=1:2),ggplot2::aes(x,y,fill=x))+ggplot2::geom_col()
stopifnot(identical(ggplot2::ggplot_build(plot)$data[[1]]$fill,report_palette(2)))
new_case <- function(name) {
 d<-file.path(work,name);dir.create(d);dir.create(file.path(d,'src'))
 file.copy(list.files(file.path(root,'subworkflows/CT_META_CAAS/local'),pattern='Rmd$',full.names=TRUE),d)
 file.copy(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'),file.path(d,'src'))
 file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src/report_presentation.R'),file.path(d,'src'))
 d
}
render <- function(d,name,p) {
 e<-new.env(parent=globalenv())
 rmarkdown::render(file.path(d,name),params=p,quiet=TRUE,envir=e)
 stopifnot(file.exists(file.path(d,sub('Rmd$','html',name))))
 html<-paste(readLines(file.path(d,sub('Rmd$','html',name)),warn=FALSE),collapse='\n')
 stopifnot(grepl('Helvetica Neue',html,fixed=TRUE),grepl('max-width: 95%',html,fixed=TRUE))
 e
}
for(case in c('discovery','empty','legacy_background')) {
 d<-new_case(paste0('pattern_',case))
 disc<-data.frame(gene=c('101','101','002'),position=c(0,0,1),trait=c('H1','H2','H1'),
                  pattern=c(1,1,2),is_conserved_meta=c(FALSE,TRUE,FALSE),conserved_pair='0:',shared_residue_count=c(0,1,0))
 if(case=='empty') disc<-disc[FALSE,]
 readr::write_tsv(disc,file.path(d,'discovery.tab'))
 bg<-if(case=='empty')character()else c('101','002','101','')
 if(case=='legacy_background')bg<-c('Gene',bg)
 writeLines(bg,file.path(d,'background_genes.output'))
 e<-render(d,'7.CAAS_pattern_annotation.Rmd',list(discovery_input='discovery.tab',background_input='background_genes.output'))
 stopifnot(nrow(e$background_df)==if(case=='empty')0 else 2,
           e$n_positions==if(case=='empty')0 else 2)
 if(case!='empty')stopifnot(e$pattern_summary$Discovery_rows[1]==2,e$pattern_summary$Unique_positions[1]==1,
                           identical(e$discovery_df$Gene,c('101','101','002')),
                           e$discovery_df$is_conserved_meta[2]=='TRUE',e$discovery_df$conserved_pair[2]=='0:',
                           e$discovery_df$shared_residue_count[2]=='1')
 cat('PASS pattern:',case,'\n')
}
for(case in c('baseline','top_only','bottom_only','empty','missing_q','missing_p')) {
 d<-new_case(paste0('significance_',case))
 pos<-data.frame(Gene=c('101','101','002'),Position=c(0,0,1),side=c('top','bottom','top'),
                 CAAS_score=c(.8,.2,.4),p.emp=c(.01,.01,NA_real_),p.adj_bh=c(.02,.02,NA_real_))
 genes<-data.frame(Gene=c('101','002'),gene_caas_score=c(.8,.4),gene_caas_score_top=c(.8,.4),
                   gene_caas_score_bottom=c(.2,0),caas_p.gen=c(.01,NA_real_),caas_p.gen.adj=c(.02,NA_real_))
 for(sf in c('_top','_bottom'))for(n in c('caas_p.gen','caas_p.gen.adj')) genes[[paste0(n,sf)]]<-genes[[n]]
 if(case=='top_only')pos<-pos[pos$side=='top',]
 if(case=='bottom_only')pos<-pos[pos$side=='bottom',]
 if(case=='empty'){pos<-pos[FALSE,];genes<-genes[FALSE,]}
 if(case=='missing_q')pos$p.adj_bh<-NA_real_
 if(case=='missing_p'){pos$p.emp<-NA_real_;pos$p.adj_bh<-NA_real_}
 readr::write_tsv(pos,file.path(d,'position_scores.tsv'));readr::write_tsv(genes,file.path(d,'gene_scores.tsv'))
 readr::write_tsv(data.frame(gene=rep('101',3),ranking=c('global','top','bottom'),N_valid=49),file.path(d,'gene_evidence.tsv'))
 jsonlite::write_json(list(gene_stat='max',nominal_threshold=.05,fdr_threshold=.05,N=99,genes=2),
                      file.path(d,'gene_evidence_metadata.json'),auto_unbox=TRUE)
 e<-render(d,'15.CAAS_significance_report.Rmd',list(position_scores_input='position_scores.tsv',gene_scores_input='gene_scores.tsv'))
 stopifnot(nrow(e$views_pos$global)==nrow(unique(pos[c('Gene','Position')])),
           all(is.na(e$resolution$N_roster[1])),!file.exists(file.path(d,'significance/pattern_significance_summary.tsv')))
 if(nrow(pos)) {
  stopifnot(all(e$views_pos$global$Gene %in% c('101','002')),
            isTRUE(all.equal(e$views_pos$global$p.emp,pos$p.emp[match(e$views_pos$global$Gene,pos$Gene)])))
 }
 if(case=='baseline')stopifnot(e$views_pos$global$CAAS_score[1]==.8,
                              sum(e$views_pos$global$BH %in% TRUE)==1,
                              e$views_gene$global$N_valid[1]==49,
                              e$views_gene$global$min_resolvable_p[1]==1/50)
 if(case=='missing_p')stopifnot(all(as.character(e$views_pos$global$evidence_status)=='p unavailable'))
 if(case=='missing_q')stopifnot(as.character(e$views_pos$global$evidence_status[1])=='p available; q unavailable')
 stopifnot(length(list.files(file.path(d,'significance'),pattern='tsv$'))==10)
 cat('PASS significance:',case,'\n')
}
# All seven curated reports retain the same page design and load paletteer directly
# or through the shared style. Special tree/heatmap dimensions remain figure-specific.
reports<-c(list.files(file.path(root,'subworkflows/TRAIT_ANALYSIS/local'),pattern='^[0-4]\\..*Rmd$',full.names=TRUE),
           list.files(file.path(root,'subworkflows/CT_META_CAAS/local'),pattern='Rmd$',full.names=TRUE))
stopifnot(length(reports)==7)
for(f in reports) {
 s<-paste(readLines(f),collapse='\n')
 stopifnot(grepl('theme: flatly',s,fixed=TRUE),grepl('Helvetica Neue',s,fixed=TRUE),
           grepl('max-width: 95%',s,fixed=TRUE),grepl('font-size: 38px',s,fixed=TRUE),
           grepl('font-size: 24px',s,fixed=TRUE),grepl('font-size: 16px',s,fixed=TRUE))
}
cat('PASS shared paletteer defaults, plot fonts and seven-report HTML style audit\n')
