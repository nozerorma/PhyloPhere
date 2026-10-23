# RER missingness is evaluated independently in every null cycle.
suppressPackageStartupMessages({library(readr);library(dplyr)})
src<-normalizePath('subworkflows/ENRICHMENT/local/src')
source(file.path(src,'evidence_flags.R'));source(file.path(src,'fcs_enrich.R'))
wd<-tempfile('rer_fcs_cli_');dir.create(wd);dir.create(file.path(wd,'gmts'))
genes<-paste0('g',1:40);set.seed(11)
vals<-setNames(rnorm(40),genes);vals[1:10]<-NA
m<-matrix(rnorm(40*25),40,dimnames=list(genes,paste0('p',1:25)));m[11:15,1:7]<-NA
sets<-list(A=genes[1:20],B=genes[15:40])
writeLines(vapply(names(sets),function(n)paste(c(n,'description',sets[[n]]),collapse='\t'),''),file.path(wd,'gmts/db.gmt'))
write_tsv(data.frame(gene=genes,score_global=vals),file.path(wd,'stats.tsv'));writeLines(genes,file.path(wd,'universe.txt'))
saveRDS(list(corStat=m),file.path(wd,'null.rds'))
gmts<-fcs_load_gmts(file.path(wd,'gmts'))
ref<-fcs_run_all(list(global=fcs_build_vals(vals,genes,score_source='rer')),gmts,perms_file=file.path(wd,'null.rds'),score_source='rer')
cli<-function(...) {
 a<-c(file.path(src,'fcs_blocks.R'),...)
 log<-system2(file.path(R.home('bin'),'Rscript'),a,stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(log,'status')))stop(paste(log,collapse='\n'))
}
cli('rer_observed','--obs',file.path(wd,'stats.tsv'),'--universe',file.path(wd,'universe.txt'),
    '--gmt-dir',file.path(wd,'gmts'),'--out',file.path(wd,'observed.rds'))
cli('rer_null','--matrices',file.path(wd,'null.rds'),'--observed',file.path(wd,'observed.rds'),
    '--obs',file.path(wd,'stats.tsv'),'--universe',file.path(wd,'universe.txt'),'--gmt-dir',file.path(wd,'gmts'),'--out',file.path(wd,'counts.tsv'))
cli('close','--observed',file.path(wd,'observed.rds'),'--counts',file.path(wd,'counts.tsv'),'--out',file.path(wd,'result.tsv'))
got<-read_tsv(file.path(wd,'result.tsv'),show_col_types=FALSE);got<-got[match(ref$pathway,got$pathway),]
stopifnot(isTRUE(all.equal(ref$p.perm,got$p.perm,tolerance=1e-12)),identical(ref$sig_perm,got$sig_perm))
cat('PASS RER CLI preserves independently evaluable null genes, bilateral AUC and native counts\n')
unlink(wd,recursive=TRUE)
