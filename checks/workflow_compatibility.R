# Run serially in the phylophere environment; temporary artifacts stay outside the checkout.
root <- normalizePath(commandArgs(TRUE)[1])
work <- normalizePath(commandArgs(TRUE)[2])
fail <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
source(file.path(root,'subworkflows/RERCONVERGE/local/rer_contract.R'))
source(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/contrast_inputs.R'))
source(file.path(root,'subworkflows/SCORING/local/src/caas_block_contract.R'))
source(file.path(root,'subworkflows/SCORING/local/src/scoring_io.R'))
Sys.unsetenv('PHYLOPHERE_VERIFY_CHECKSUMS')
setwd(work)
# Parse edited R and report chunks before the integration checks.
paths <- c(list.files(file.path(root,'subworkflows/RERCONVERGE/local'),pattern='[.]R$',full.names=TRUE),
 list.files(file.path(root,'subworkflows/SCORING/local/src'),pattern='[.]R$',full.names=TRUE),
 file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/contrast_inputs.R'),
 file.path(root,'subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R'))
for(p in paths) invisible(parse(p))
for(p in c('CT_META_CAAS/local/15.CAAS_significance_report.Rmd','CT_POSTPROC/local/8.Characterization_report.Rmd',
 'SCORING/local/10.Scoring_report.Rmd','TRAIT_ANALYSIS/local/3.CI-composition.Rmd','TRAIT_ANALYSIS/local/4.Independent_contrasts.Rmd','FADE/local/6.FADE_report.Rmd')) {
 lines <- readLines(file.path(root,'subworkflows',p)); code <- character(); active <- FALSE
 for(line in lines) {
   if(startsWith(line,'```{r')) {active<-TRUE;next}
   if(active && startsWith(line,'```')) {active<-FALSE;next}
   if(active) code<-c(code,line)
 }
 invisible(parse(text=code))
}
cat('PASS R and report syntax\n')
# Native RER estimator: legacy values are not trusted or silently labelled weighted.
suppressPackageStartupMessages({library(ape);library(RERconverge)})
set.seed(61);tr<-rtree(16,rooted=FALSE,tip.label=LETTERS[1:16])
writeLines(vapply(1:24,function(i) {tr$edge.length<-runif(nrow(tr$edge),.2,2);paste0('g.',i,'\t',write.tree(tr))},''),'trees.nwk')
trees<-readTrees('trees.nwk');species<-LETTERS[1:16]
reference<-getAllResiduals(trees,useSpecies=species,transform='sqrt',use.weights=TRUE,norm='scale',min.sp=3)
legacy<-reference; legacy[]<-0
recovered<-suppressWarnings(rer_compatible_matrix(legacy,trees,species,3))
stopifnot(isTRUE(all.equal(as.numeric(recovered),as.numeric(reference))),
 attr(recovered,'rer_contract')$provenance_status=='recomputed_from_legacy_inputs')
wrong<-recovered;attr(wrong,'rer_contract')$tree<-'another-tree'
fail(rer_compatible_matrix(wrong,trees,species,3))
wrong_axis<-legacy;colnames(wrong_axis)<-rev(rer_tree_contract(trees)$paths)
fail(rer_compatible_matrix(wrong_axis,trees,species,3))
cat('PASS legacy RER reconstruction with native weighted estimator and mismatched axes rejection\n')
# Exercise the actual inference and table CLIs with a legacy matrix.
saveRDS(trees,'trees.rds');saveRDS(legacy,'legacy.rds')
trait_vector<-setNames(seq_along(species),species);n_vector<-c_vector<-NULL
save(trait_vector,n_vector,c_vector,file='traits.RData')
rs<-file.path(R.home('bin'),'Rscript')
command<-c(file.path(root,'subworkflows/RERCONVERGE/local/continuous_rer.R'),
 'traits.RData','trees.rds','paths.rds','legacy.rds','continuous.output','3','0','0','0','2','cc','none','79')
status<-system2(rs,shQuote(command),stdout='continuous.log',stderr='continuous.log')
if(status!=0L)stop(paste(readLines('continuous.log'),collapse='\n'))
result<-readRDS('continuous.output')
stopifnot(attr(result,'rer_settings')$matrix_contract$provenance_status=='recomputed_from_legacy_inputs')
command<-c(file.path(root,'subworkflows/RERCONVERGE/local/rer_tables.R'),'continuous.output','NO_FILE','trait','.05','.05')
status<-system2(rs,shQuote(command),stdout='tables.log',stderr='tables.log')
if(status!=0L)stop(paste(readLines('tables.log'),collapse='\n'))
stopifnot(file.exists('rer_inputs/file_identity.json'),!file.exists('rer_inputs/checksums.json'))
cat('PASS legacy continuous RER inference and table CLIs\n')

# Scientific block checks remain active with no byte hashes.
mm<-matrix(c(.2,.4),1,dimnames=list('NA',c('b_1','b_2')))
b<-list(genes='NA',cycles=colnames(mm),gene_stat='max',n_units=2,stat=list(global=mm,top=mm,bottom=mm*0))
saveRDS(b,'matrices.rds');idx<-caas_block_index('matrices.rds',c('b_1','b_2'))
stopifnot(is.null(idx$blocks[[1]]$md5))
Sys.setenv(PHYLOPHERE_VERIFY_CHECKSUMS='1');stopifnot(nzchar(caas_block_index('matrices.rds')$blocks[[1]]$md5));Sys.unsetenv('PHYLOPHERE_VERIFY_CHECKSUMS')
b$gene_stat<-'sum';saveRDS(b,'bad.rds');fail(caas_block_index('bad.rds'))
cat('PASS scientific block checks and optional checksum mode\n')
# Reconstructed contrast bundles are produced from current inputs, never guessed.
df<-data.frame(species=c('A','B','C'),value=c(.1,.5,.9),n=c(10,20,30),k=c(1,10,27))
ci<-contrast_trait_ci(df,'value','n','k',FALSE)
expected<-binom::binom.confint(df$k,df$n,methods='bayes',prior.shape1=.5,prior.shape2=.5,conf.level=.95)
stopifnot(identical(ci$trait,df$value),identical(ci$trait_ci_lb,expected$lower),identical(ci$category,c('bottom','intermediate','top')))
pv<-list(version=1L,tree='current',seed=1L)
rebuilds<-0L
rebuild<-function(){rebuilds<<-rebuilds+1L;list(trait_ci=ci,ctx=list(selected_model='BM'))}
x<-suppressWarnings(read_observed_contrasts('selection.rds',pv,rebuild));y<-read_observed_contrasts('selection.rds',pv,rebuild)
stopifnot(rebuilds==1L,identical(x,y),x$provenance_status=='recomputed_from_current_inputs')
fail(read_observed_contrasts('selection.rds',list(tree='different'),rebuild))
cat('PASS contrast reconstruction, interval equivalence and provenance mismatch rejection\n')
# Render full significance reports with literal IDs, with and without metadata.
for(known in c(TRUE,FALSE)) {
 d<-file.path(work,paste0('meta_',known));dir.create(d);dir.create(file.path(d,'src'))
 file.copy(file.path(root,'subworkflows/CT_META_CAAS/local/15.CAAS_significance_report.Rmd'),d)
 file.copy(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'),file.path(d,'src'))
 file.copy(file.path(root,'subworkflows/ENRICHMENT/local/src/report_presentation.R'),file.path(d,'src'))
 ids<-c('NA','001','"Q"')
 pos<-data.frame(Gene=ids,Position=c(0,1,2),side='top',CAAS_score=c(.8,.5,.3),p.emp=c(.01,.1,NA),p.adj_bh=c(.03,.3,NA))
 genes<-data.frame(Gene=ids,gene_caas_score=pos$CAAS_score,gene_caas_score_top=pos$CAAS_score,gene_caas_score_bottom=0,
 caas_p.gen=pos$p.emp,caas_p.gen.adj=pos$p.adj_bh)
 for(sf in c('_top','_bottom'))for(n in c('caas_p.gen','caas_p.gen.adj'))genes[[paste0(n,sf)]]<-genes[[n]]
 scoring_write(pos,file.path(d,'position_scores.tsv'));scoring_write(genes,file.path(d,'gene_scores.tsv'))
 if(known) jsonlite::write_json(list(gene_stat='max',nominal_threshold=.05,fdr_threshold=.05),file.path(d,'gene_evidence_metadata.json'),auto_unbox=TRUE)
 e<-new.env(parent=globalenv())
 rmarkdown::render(file.path(d,'15.CAAS_significance_report.Rmd'),quiet=TRUE,envir=e,
 params=list(position_scores_input='position_scores.tsv',gene_scores_input='gene_scores.tsv'))
 stopifnot(setequal(e$pos$Gene,ids),identical(e$pos$p.emp,pos$p.emp),identical(e$pos$p.adj_bh,pos$p.adj_bh))
 if(!known)stopifnot(all(is.na(e$views_pos$global$BH)),all(is.na(e$summary$BH)))
 else stopifnot(sum(e$views_pos$global$BH %in% TRUE)==1L)
 cat('PASS complete significance report, original metadata:',known,'\n')
}
# Render CT_POSTPROC with actual removals but no original settings or universe ledger.
d<-file.path(work,'postproc');dir.create(d);dir.create(file.path(d,'src'));dir.create(file.path(d,'sweep'))
file.copy(file.path(root,'subworkflows/CT_POSTPROC/local/8.Characterization_report.Rmd'),d)
file.copy(file.path(root,'subworkflows/CT_POSTPROC/local/src/report_io.R'),file.path(d,'src'))
file.copy(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'),file.path(d,'src'))
disc<-data.frame(Gene=c('NA','001'),Position=c(0,1),side=c('top','none'),asr_path_score=c(.8,0))
scoring_write(disc,file.path(d,'input.tsv'));scoring_write(disc[2,,drop=FALSE],file.path(d,'retained.tsv'))
scoring_write(data.frame(Gene=disc$Gene,length=10,length_source='map',n_CAAS=1,n_CAAS_per_length=10,category='Normal'),file.path(d,'stats.tsv'))
scoring_write(data.frame(Gene=character(),category=character()),file.path(d,'removed.tsv'))
scoring_write(data.frame(Parameter=character(),Minlen=numeric(),Maxcaas=numeric(),DiscardedCount=integer()),file.path(d,'summary.tsv'))
scoring_write(data.frame(Gene=disc$Gene,Position=disc$Position,clustering_flag='Good'),file.path(d,'selected.tsv'))
scoring_write(data.frame(gene=disc$Gene,length=10),file.path(d,'annotation.tsv'))
e<-new.env(parent=globalenv())
rmarkdown::render(file.path(d,'8.Characterization_report.Rmd'),quiet=TRUE,envir=e,params=list(discovery_file='input.tsv',filtered_discovery_file='retained.tsv',
 gene_stats_file='stats.tsv',removed_genes_file='removed.tsv',filter_summary_file='summary.tsv',selected_cluster_file='selected.tsv',gene_ensembl_file='annotation.tsv',filter_dir='sweep'))
stopifnot(nrow(e$removal_events)==1,e$removal_events$Gene=='NA',e$removal_events$removal_stage=='unknown',!e$known_universe)
cat('PASS complete legacy postproc report without invented removal causes\n')
# Full SCORING report remains descriptive when original calibration metadata are absent.
d<-file.path(work,'scoring');dir.create(d);dir.create(file.path(d,'src'))
file.copy(file.path(root,'subworkflows/SCORING/local/10.Scoring_report.Rmd'),d)
for(f in c('scoring_io.R','report_evidence.R'))file.copy(file.path(root,'subworkflows/SCORING/local/src',f),file.path(d,'src'))
file.copy(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'),file.path(d,'src'))
genes<-data.frame(Gene=c('NA','001','g3'),n_positions=c(2L,1L,0L),gene_caas_score=c(.8,.4,0),
 gene_caas_score_top=c(.8,.4,0),gene_caas_score_bottom=c(.2,0,0),caas_p.gen=c(.01,.5,NA_real_))
positions<-data.frame(Gene=c('NA','NA','001'),Position=c(0L,1L,0L),side=c('top','bottom','top'),
 CAAS_score=c(.8,.2,.4),p.emp=c(.01,.5,.2),p.adj_bh=c(.04,.5,.4))
scoring_write(genes,file.path(d,'gene_scores.tsv'));scoring_write(positions,file.path(d,'position_scores.tsv'))
scoring_write(data.frame(Gene=character(),correlation=numeric()),file.path(d,'correlations.tsv'))
e<-new.env(parent=globalenv())
rmarkdown::render(file.path(d,'10.Scoring_report.Rmd'),quiet=TRUE,envir=e,params=list(position_scores_file='position_scores.tsv',gene_scores_file='gene_scores.tsv',gene_corr_file='correlations.tsv',output_dir=d))
stopifnot(!e$.known_metadata,identical(e$pos_df$p.emp,positions$p.emp))
cat('PASS complete legacy SCORING report\n')
# An opt-in audit can read metadata produced without checksums: missing is not mismatch.
jsonlite::write_json(list(gene_stat='max',nominal_threshold=.05,fdr_threshold=.05,table_md5=list()),file.path(d,'gene_evidence_metadata.json'),auto_unbox=TRUE)
Sys.setenv(PHYLOPHERE_VERIFY_CHECKSUMS='1')
e<-new.env(parent=globalenv())
rmarkdown::render(file.path(d,'10.Scoring_report.Rmd'),quiet=TRUE,envir=e,params=list(position_scores_file='position_scores.tsv',gene_scores_file='gene_scores.tsv',gene_corr_file='correlations.tsv',output_dir=d))
Sys.unsetenv('PHYLOPHERE_VERIFY_CHECKSUMS')
cat('PASS SCORING report integrity opt-in with unavailable historical checksums\n')

# FCS compatibility with historical hash manifests preserves block-axis checks.
d<-file.path(work,'fcs');dir.create(d);setwd(d)
dir.create('live/saved/blocks/blk_000001',recursive=TRUE)
b<-list(genes='NA',cycles=c('b_1','b_2'),gene_stat='max',n_units=2,stat=list(global=mm,top=mm,bottom=mm*0))
saveRDS(b,'live/saved/blocks/blk_000001/matrices.rds');writeLines(b$cycles,'live/saved/cycle_roster.txt')
idx<-caas_block_index('live/saved/blocks/blk_000001/matrices.rds',b$cycles)
idx$blocks[[1]]$md5<-'old-checksum'
jsonlite::write_json(idx,'live/saved/null_metadata.json',auto_unbox=TRUE)
command<-c(file.path(root,'subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R'),'caas','max','2','enabled')
status<-system2(rs,shQuote(command),stdout='prepare.log',stderr='prepare.log')
if(status!=0L)stop(paste(readLines('prepare.log'),collapse='\n'))
unlink('blocks',recursive=TRUE)
idx$blocks[[1]]$md5<-NULL
jsonlite::write_json(idx,'live/saved/null_metadata.json',auto_unbox=TRUE)
Sys.setenv(PHYLOPHERE_VERIFY_CHECKSUMS='1')
status<-system2(rs,shQuote(command),stdout='missing-checksum.log',stderr='missing-checksum.log')
if(status!=0L)stop(paste(readLines('missing-checksum.log'),collapse='\n'))
idx$blocks[[1]]$md5<-'old-checksum'
jsonlite::write_json(idx,'live/saved/null_metadata.json',auto_unbox=TRUE)
status<-system2(rs,shQuote(command),stdout='audit.log',stderr='audit.log');stopifnot(status!=0L)
Sys.unsetenv('PHYLOPHERE_VERIFY_CHECKSUMS')
idx$blocks[[1]]$genes<-'another-gene'
jsonlite::write_json(idx,'live/saved/null_metadata.json',auto_unbox=TRUE)
status<-system2(rs,shQuote(command),stdout='bad-axis.log',stderr='bad-axis.log');stopifnot(status!=0L)
cat('PASS FCS historical manifests, integrity opt-in and semantic mismatch rejection\n')
