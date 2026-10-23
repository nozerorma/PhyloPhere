a<-commandArgs(TRUE);root<-normalizePath(a[1]);scratch<-a[2];dir.create(scratch,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'subworkflows/ENRICHMENT/local/src/gene_p.R'))
write<-function(d,f)write.table(d,f,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
g<-c('g1','g2','g3','g4');S<-matrix(c(.9,.2,.8,0, .1,.3,.8,0, .4,.1,.8,0, .6,.3,.8,0),4,dimnames=list(g,paste0('b_',1:4)))
obs<-c(.95,.2,.8,0)
for(batch in c(1,3)) {
 d<-file.path(scratch,paste0('batch',batch));dir.create(d);dir.create(file.path(d,'blocks'));setwd(d)
 write(data.frame(Gene=g,gene_caas_score=obs),'genes.tsv');write(data.frame(Gene='g1',Position=1,side='top',p.emp=.01,p.adj_bh=.2),'positions.tsv')
 write(data.frame(gene=g,score_global=obs,score_top=obs,score_bottom=obs),'stats.tsv');writeLines(g,'universe.txt')
 for(i in 1:2)saveRDS(list(gene_stat='size_adj_max',stat=setNames(rep(list(S[,(i*2-1):(i*2),drop=FALSE]),3),c('global','top','bottom'))),sprintf('blocks/b%d.rds',i))
 result<-system2(file.path(R.home('bin'),'Rscript'),c(file.path(root,'subworkflows/SCORING/local/src/scoring_evidence.R'),'--genes genes.tsv --positions positions.tsv --stats stats.tsv --universe universe.txt --blocks blocks --batch',batch),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(result,'status')))stop(paste(result,collapse='\n'))
 tab<-read.delim('gene_evidence.tsv');expected<-gene_p(S,obs)
 stopifnot(isTRUE(all.equal(tab$p.gen[tab$ranking=='global'],expected$p_obs,check.attributes=FALSE)),isTRUE(all.equal(tab$p.gen.adj[tab$ranking=='global'],p.adjust(expected$p_obs,'BH'),check.attributes=FALSE)))
 # Reassemble the null of each rank from batch slices: exact LOO identity.
 z<-matrix(NA_real_,4,4,dimnames=dimnames(S))
 for(f in list.files('calibration',pattern='^global_gp_',full.names=TRUE)) {b<-readRDS(f);z[b$rows,b$cols]<-b$values}
 stopifnot(isTRUE(all.equal(z,expected$nlog10_null)))
 stopifnot(all(tab$N_roster==4),!any(tab$adjusted & !tab$nominal))
 if(batch==1)reference<-tab else stopifnot(identical(tab,reference))
}
cat('PASS shared gene evidence: observed, BH full universe, exact LOO, batch invariance\n')
# Boundary, missing-null, empty-universe, NA-position and duplicate-side cases.
for(case in c('no_null','one_cycle','boundary','empty_universe','na_position')) {
 d<-file.path(scratch,case);dir.create(d);dir.create(file.path(d,'blocks'));setwd(d)
 write(data.frame(Gene='g1',gene_caas_score=1),'genes.tsv')
 write(data.frame(Gene=c('g1','g1'),Position=c(1,1),side=c('top','bottom'),p.emp=if(case=='na_position') NA else .05,p.adj_bh=if(case=='na_position') NA else .05),'positions.tsv')
 write(data.frame(gene='g1',score_global=1,score_top=1,score_bottom=1),'stats.tsv')
 writeLines(if(case=='empty_universe')character() else 'g1','universe.txt')
 N<-if(case=='boundary')19L else if(case=='one_cycle')1L else 0L
 if(N) {
  m<-matrix(0,1,N,dimnames=list('g1',paste0('b_',seq_len(N))))
  saveRDS(list(gene_stat='size_adj_max',stat=setNames(rep(list(m),3),c('global','top','bottom'))),'blocks/b.rds')
 }
 result<-system2(file.path(R.home('bin'),'Rscript'),c(file.path(root,'subworkflows/SCORING/local/src/scoring_evidence.R'),'--genes genes.tsv --positions positions.tsv --stats stats.tsv --universe universe.txt --blocks blocks'),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(result,'status')))stop(paste(result,collapse='\n'))
 tab<-read.delim('gene_evidence.tsv')
 stopifnot(!any(tab$nominal),!any(tab$adjusted),nrow(read.delim('position_lists/active_global.tsv'))==0)
 if(case=='na_position')stopifnot(nrow(read.delim('position_lists/active_top.tsv'))==0)
 if(case=='boundary')stopifnot(all(tab$p.gen==.05),all(tab$p.gen.adj==.05)) else stopifnot(file.exists('calibration/empty.tsv'))
 if(case=='empty_universe')stopifnot(nrow(tab)==0) else if(case!='boundary')stopifnot(all(is.na(tab$p.gen)))
}
cat('PASS: strict boundary, no null, insufficient cycles, empty universe, NA position, position-side consistency\n')
# Incompatible null: a block calibrated with a different gene statistic fails informatively.
d<-file.path(scratch,'stat_mismatch');dir.create(d);dir.create(file.path(d,'blocks'));setwd(d)
write(data.frame(Gene='g1',gene_caas_score=1),'genes.tsv')
write(data.frame(Gene='g1',Position=1,side='top',p.emp=.01,p.adj_bh=.2),'positions.tsv')
write(data.frame(gene='g1',score_global=1,score_top=1,score_bottom=1),'stats.tsv')
writeLines('g1','universe.txt')
saveRDS(list(gene_stat='max',stat=setNames(rep(list(matrix(0,1,2,dimnames=list('g1',paste0('b_',1:2)))),3),c('global','top','bottom'))),'blocks/b.rds')
result<-system2(file.path(R.home('bin'),'Rscript'),c(file.path(root,'subworkflows/SCORING/local/src/scoring_evidence.R'),'--genes genes.tsv --positions positions.tsv --stats stats.tsv --universe universe.txt --blocks blocks'),stdout=TRUE,stderr=TRUE)
stopifnot(!is.null(attr(result,'status')),any(grepl('Gene statistic mismatch',result)))
cat('PASS: incompatible null rejected with diagnostic\n')
