# Independent WP9 variant, inclusive-LOO and native AUC oracles; no scientific default selected.
suppressPackageStartupMessages({library(dplyr);library(RERconverge)})
source('adhoc_scripts/wp9_gene_stat_control.R')
source('subworkflows/ENRICHMENT/local/src/evidence_flags.R')
source('subworkflows/ENRICHMENT/local/src/fcs_enrich.R')
source('subworkflows/ENRICHMENT/local/src/gene_p.R')
source('subworkflows/SCORING/local/src/caas_block_contract.R')
same<-function(a,b)stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
reject<-function(x)stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
set.seed(91);genes<-paste0('g',1:40);cycles<-paste0('b_',1:11)
n<-matrix(sample(0:4,440,TRUE),40,dimnames=list(genes,cycles))
f<-matrix(sample(c(.1,.25,.5,1),440,TRUE),40,dimnames=list(genes,cycles));f[n==0]<-0
m<-matrix(sample(c(.125,.25,.5,1),440,TRUE),40,dimnames=list(genes,cycles));m[n==0]<-0
parts<-list(max=m,fmax=f,n=n)
for(v in wp9_variants) {
 x<-wp9_stat(parts,v)
 beta<-switch(v,fmax=0,beta025=.25,beta05=.5,beta075=.75,size_adj_max=1)
 for(i in 1:40)for(j in 1:11) {
  expected<-if(n[i,j]==0)0 else if(v=='max') m[i,j] else f[i,j]^(n[i,j]^beta)
  same(x[i,j],expected)
 }
}
stopifnot(identical(wp9_stat(parts,'fmax'),f))
same(wp9_stat(parts,'size_adj_max')[n>0],f[n>0]^n[n>0])
reject(wp9_stat(parts,'unknown'))
bad<-parts;bad$n[1,1]<-1.5;reject(wp9_stat(bad,'fmax'))
bad<-parts;rownames(bad$fmax)<-rev(genes);reject(wp9_stat(bad,'fmax'))
# Near ties: gene tolerance retained, pathways exact; native AUC checked per column.
obs<-setNames(sample(c(0,.25,.5,1),40,TRUE),genes)
gmts<-list(db=list(genesets=list(A=genes[1:12],B=genes[10:28],C=genes[25:40]),geneset.names=c('A','B','C')))
x<-gene_p(wp9_stat(parts,'fmax'),obs)
got<-wp9_fcs(x$nlog10_null,x$nlog10_obs,genes,gmts)
null_native<-matrix(NA_real_,nrow(got$observed),11)
for(j in 1:11) {
 invisible(capture.output(native<-RERconverge::fastwilcoxGMTall(setNames(x$nlog10_null[,j],genes),gmts,outputGeneVals=FALSE,num.g=10)))
 null_native[,j]<-native$db$stat[match(got$observed$pathway,rownames(native$db))]
}
same(got$stats$db[got$observed$pathway,,drop=FALSE],null_native)
expected<-vapply(seq_len(nrow(null_native)),function(i)(1+sum(null_native[i,]>=got$observed$stat[i]))/12,numeric(1))
same(got$observed$p.perm,expected)
same(got$observed$q.perm,p.adjust(expected,'BH'))
null_p<-matrix(NA_real_,nrow(null_native),11)
for(i in seq_len(nrow(null_p)))for(j in 1:11) null_p[i,j]<-(1+sum(null_native[i,-j]>=null_native[i,j]))/11
same(got$cycles$p05,colSums(null_p<=.05))
same(got$cycles$q15,apply(null_p,2,function(p)sum(p.adjust(p,'BH')<.15)))
# A pathway difference below gene tolerance still compares exactly at the pathway stage.
tiny<-matrix(c(.1,.1+5e-13,.2),1,dimnames=list('A',c('a','b','c')))
same(gene_p(tiny,.1+5e-13,tie_tol=0)$p_obs,3/4)
# Disjoint references: exhaustive query comparisons and independent native evaluation AUC.
reference<-wp9_stat(parts,'fmax')[,1:5];query<-wp9_stat(parts,'fmax')[,6:11]
pquery<-wp9_reference_p(reference,query,1e-12)
for(g in 1:40)for(j in 1:6)same(pquery[g,j],(1+sum(reference[g,]>=query[g,j]-1e-12))/6)
real<-list(db=as.data.frame(got$observed));rownames(real$db)<-real$db$pathway
check<-wp9_holdout(wp9_stat(parts,'fmax'),genes,gmts,real)
stopifnot(setequal(check$cycle,cycles),!anyDuplicated(check$cycle),all(check$reference_cycles%in%c(5,6)))
raw<-wp9_stat(parts,'fmax')
for(k in 1:2) {
 ev<-check$cycle[check$fold==k];ref<-setdiff(cycles,ev)
 rscore<-matrix(0,40,length(ref),dimnames=list(genes,ref));escore<-matrix(0,40,length(ev),dimnames=list(genes,ev))
 for(g in 1:40) {
  for(j in seq_along(ref))rscore[g,j]<--log10((1+sum(raw[g,ref[-j]]>=raw[g,ref[j]]-1e-12))/length(ref))
  for(j in seq_along(ev))escore[g,j]<--log10((1+sum(raw[g,ref]>=raw[g,ev[j]]-1e-12))/(length(ref)+1))
 }
 native_stats<-function(scores) {
  out<-matrix(NA_real_,3,ncol(scores),dimnames=list(c('A','B','C'),colnames(scores)))
  for(j in seq_len(ncol(scores))) {
   invisible(capture.output(z<-RERconverge::fastwilcoxGMTall(setNames(scores[,j],genes),gmts,outputGeneVals=FALSE,num.g=10)))
   out[,j]<-z$db$stat[match(rownames(out),rownames(z$db))]
  };out
 }
 ra<-native_stats(rscore);ea<-native_stats(escore)
 for(j in seq_along(ev)) {
  p<-vapply(1:3,function(i)(1+sum(ra[i,]>=ea[i,j]))/(length(ref)+1),numeric(1))
  same(check$p05[check$cycle==ev[j]],sum(p<=.05))
  same(check$q15[check$cycle==ev[j]],sum(p.adjust(p,'BH')<.15))
 }
}

# Canonical all-zero genes/cycles, roster order and b0 must survive reconstruction.
dir<-tempfile('wp9_store_');dir.create(dir);dir.create(file.path(dir,'blocks'));dir.create(file.path(dir,'b0'))
writeLines(cycles,file.path(dir,'cycle_roster.txt'))
for(k in 1:2) {
 sel<-if(k==1)1:6 else 7:11
 dest<-file.path(dir,'blocks',sprintf('blk_%04d',k));dir.create(dest)
 b<-list(genes=genes,cycles=cycles[sel],gene_stat='size_adj_max',n_units=1,
   parts=lapply(parts,function(p)setNames(rep(list(p[,sel,drop=FALSE]),3),c('global','top','bottom'))))
 b$stat<-setNames(rep(list(wp9_stat(lapply(parts,function(p)p[,sel,drop=FALSE]),'size_adj_max')),3),c('global','top','bottom'))
 saveRDS(b,file.path(dest,'matrices.rds'));writeLines(cycles[sel],file.path(dest,'cycles.txt'))
}
b0<-b;b0$cycles<-'b_0';b0$parts<-lapply(parts,function(p){v<-p[,1,drop=FALSE];colnames(v)<-'b_0';setNames(rep(list(v),3),c('global','top','bottom'))})
b0$stat<-setNames(lapply(c('global','top','bottom'),function(rk)wp9_stat(lapply(b0$parts,function(p)p[[rk]]),'size_adj_max')),c('global','top','bottom'))
saveRDS(b0,file.path(dir,'b0/matrices.rds'))
data<-wp9_store(dir,genes,'global');same(data$parts$global$fmax,f[sort(genes),])
subset<-wp9_store(dir,genes[-1],'global');stopifnot(genes[1]%in%subset$excluded_store_genes);same(subset$parts$global$fmax,f[sort(genes[-1]),])
reject(wp9_store(dir,genes,'global',1))
jsonlite::write_json(list(gene_stat='max'),file.path(dir,'null_metadata.json'));reject(wp9_store(dir,genes,'global'));unlink(file.path(dir,'null_metadata.json'))
file.rename(file.path(dir,'blocks/blk_0002/matrices.rds'),file.path(dir,'removed.rds'))
reject(wp9_store(dir,genes,'global'))
unlink(dir,recursive=TRUE)
# Legacy pools include scored positions outside the pathway universe; zero cycles come from the RDS roster.
legacy<-tempfile('wp9_legacy_');dir.create(legacy)
positions<-data.frame(Gene=c('g1','g2','g1','g2'),Position=1L,side=c('top','top','bottom','top'),
 cycle=c('b_1','b_1','b_3','b_3'),s=c(.2,.8,.3,.9))
saveRDS(positions,file.path(legacy,'pos.rds'))
write.table(data.frame(Gene=c('g1','g2'),Position=1L,side='top',caas_score=c(.4,.7)),file.path(legacy,'b0.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
reference<-matrix(c(.5,1,0,0,.5,1),2,dimnames=list(c('g1','g2'),c('b_1','b_2','b_3')))
saveRDS(list(caas_corStat_byrank=list(global=reference)),file.path(legacy,'null.rds'))
parsed<-wp9_legacy(file.path(legacy,'pos.rds'),file.path(legacy,'b0.tsv'),file.path(legacy,'null.rds'),'g1','global')
same(as.numeric(parsed$parts$global$fmax),c(.5,0,.5));stopifnot('g2'%in%parsed$excluded_store_genes)
unlink(legacy,recursive=TRUE)
cat('PASS WP9 variants, zero detection, native pathway AUC, exact pathway ties, exhaustive LOO/BH, canonical reconstruction and incomplete-store rejection\n')
