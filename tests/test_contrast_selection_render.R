# Actual report renders, with complete export checks, in an isolated directory.
options(phylo_debug=FALSE)
repo <- normalizePath('.')
local <- file.path(repo,'subworkflows/TRAIT_ANALYSIS/local')
work <- tempfile('contrast-render-');dir.create(work)
invisible(file.copy(list.files(local,full.names=TRUE),work,recursive=TRUE))
invisible(file.copy(file.path(repo,'subworkflows/CT/local/scripts',c('pss_core.R','lean_contrast_selector.R')),file.path(work,'src')))
setwd(work);cat('ARTIFACT_DIR:',work,'\n')
tr <- ape::read.tree(text='(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):2);')
ape::write.tree(tr,'tree.nwk')
cases <- commandArgs(trailingOnly=TRUE)
if (!length(cases))cases <- c('one','two','four','six','zero','count','ordinal','headers')
for(case in cases) {
 n <- switch(case,one=1,two=2,four=4,6)
 score <- if(case=='zero')rep(.5,n) else if(case=='ordinal')c(0,2,1,2,0,1) else c(.1,.9,.2,.8,.3,.7)[seq_len(n)]
 dat <- data.frame(species=tr$tip.label[seq_len(n)],family='Unknown, family',score=score,n=100,
                   c=if(case=='zero')rep(50,n) else round(score*100))
 primary <- if(case=='headers') 'score rate' else 'score'
 if(case=='headers')names(dat)[names(dat)=='score'] <- primary
 write.csv(dat,paste0(case,'.csv'),row.names=FALSE)
 out <- paste0(case,'_results');sd <- file.path(out,'1.Data-exploration/1.Species_distribution');dir.create(sd,recursive=TRUE)
 dat$g_median <- median(score);write.csv(dat,file.path(sd,'trait_stats.csv'),row.names=FALSE)
 counts <- case %in% c('zero','count')
 p <- list(trait_file=paste0(case,'.csv'),tree_file='tree.nwk',output_dir=out,seed='1998',clade_name='primates',taxon_of_interest='family',sp_colname='species',traitname=primary,n_trait=if(counts)'n' else '',c_trait=if(counts)'c' else '',secondary_trait='',branch_trait='',pss_top_pct='0',perm_strategy=if(case=='headers')'auto' else 'bm',trait_type=if(case=='ordinal')'ordinal' else 'continuous')
 for(report in c('3.CI-composition.Rmd','4.Independent_contrasts.Rmd')) {
  pp <- p;if(startsWith(report,'4'))pp <- c(pp,list(max_contrasts='0',max_fop='3'))
  rmarkdown::render(report,params=pp,output_file=paste0(case,'-',substr(report,1,1),'.html'),envir=new.env(),quiet=TRUE)
  cat('PASS actual render:',case,report,'\n')
 }
 ci_dir <- file.path(out,'1.Data-exploration/5.CI_overlaps')
 bundle <- readRDS(file.path(ci_dir,'observed_selection.rds'))
 pairs <- read.csv(file.path(ci_dir,'pairwise_data.csv'))
 exported <- pairs[!pairs$trait_overlap,c('species1','species2')]
 cc <- bundle$candidates$cand_df
 stopifnot(nrow(exported)==if(is.null(cc))0 else nrow(cc))
 traitfile <- file.path(out,'2.CT/1.Traitfiles/traitfile.tab')
 if(case %in% c('one','two','zero'))stopifnot(file.info(traitfile)$size==0)
 else stopifnot(file.info(traitfile)$size>0)
 stopifnot(file.exists(file.path(out,'2.CT/1.Traitfiles/observed_selection.rds')),
           file.exists(file.path(out,'2.CT/2.Permulation_traitfiles/permulation_traitfile.tab')),
           file.exists(file.path(ci_dir,'candidate_species_audit.tsv')))
}
cat('Passed',2L*length(cases),'actual contrast report renders\n')
