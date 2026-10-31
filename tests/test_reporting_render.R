# Real Rmd smoke checks for small cohorts; run from the repository root.
options(phylo_debug=FALSE)
repo<-normalizePath('.')
local<-file.path(repo,'subworkflows/TRAIT_ANALYSIS/local')
work<-tempfile('reporting-render-');dir.create(work);invisible(file.copy(list.files(local,full.names=TRUE),work,recursive=TRUE))
setwd(work);cat('ARTIFACT_DIR:',work,'\n')
ape::write.tree(ape::read.tree(text='((S1:1,S2:1):1,(S3:1,S4:1):1);'),'tree.nwk')
for(n in c(1,2,3)) {
 write.csv(data.frame(species=paste0('S',seq_len(n)),family='Hominidae',score=seq_len(n)),'traits.csv',row.names=FALSE)
 p<-list(trait_file='traits.csv',tree_file='tree.nwk',output_dir=paste0('results',n),seed='1998',clade_name='primates',taxon_of_interest='family',sp_colname='species',traitname='score',n_trait='',c_trait='',secondary_trait='',branch_trait='',pss_top_pct='0.05',perm_strategy='best_model',trait_type='continuous')
 rmarkdown::render('1.Dataset_exploration.Rmd',params=p,output_file=paste0('dataset',n,'.html'),envir=new.env(),quiet=TRUE)
 stats<-read.csv(file.path(p$output_dir,'1.Data-exploration/1.Species_distribution/trait_stats.csv'))
 stopifnot(nrow(stats)==n,all(stats$g_median==median(seq_len(n))))
 if(n<3) {
  rmarkdown::render('2.Phenotype_exploration.Rmd',params=p,output_file=paste0('phenotype',n,'.html'),envir=new.env(),quiet=TRUE)
  stopifnot(file.exists(paste0('phenotype',n,'.html')))
 }
 cat('PASS actual dataset render N=',n,if(n<3)' + phenotype' else '', '\n')
}
