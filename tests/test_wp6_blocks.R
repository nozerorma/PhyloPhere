# Independent native AUC oracle plus exhaustive per-gene calibration.
suppressPackageStartupMessages({library(dplyr); library(RERconverge)})
source('subworkflows/ENRICHMENT/local/src/percentile_flags.R')
source('subworkflows/ENRICHMENT/local/src/fcs_enrich.R')
source('subworkflows/ENRICHMENT/local/src/gene_p.R')
source('subworkflows/ENRICHMENT/local/src/caas_blocks.R')
same <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
reject <- function(x) stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
set.seed(72)
genes <- paste0('g',1:400); cycles <- paste0('b_',1:25)
raw <- matrix(sample(c(0,1,2,NA),400*25,replace=TRUE),400,dimnames=list(genes,cycles))
observed <- setNames(sample(0:3,400,replace=TRUE),genes)
gp <- gene_p(raw,observed)
# This oracle counts each comparison separately and never uses sorted ranks.
for(g in seq_len(nrow(raw))) {
 x <- raw[g,]; ok <- !is.na(x)
 same(gp$p_obs[g],(1+sum(x>=observed[g]-1e-12,na.rm=TRUE))/(1+sum(ok)))
 for(j in which(ok)) same(gp$p_null[g,j],(1+sum(x[-j]>=x[j]-1e-12,na.rm=TRUE))/sum(ok))
}
gmts <- list(db=list(genesets=lapply(1:30,function(i) sample(genes,20)),geneset.names=paste0('set',1:30)))
for(source_type in c('caas','rer')) {
 if(source_type=='caas') {vals<-setNames(gp$nlog10_obs,genes); m<-gp$nlog10_null; alt<-'greater'}
 else {vals<-setNames(rnorm(400),genes);m<-matrix(rnorm(400*25),400,dimnames=list(genes,cycles));m[1,1]<-Inf;m[2,1]<--Inf;m[3,2]<-NA;alt<-'two.sided'}
 obs <- fcs_observed(list(global=vals),gmts,num_g=10,score_source=source_type)$enrich
 ref <- matrix(NA_real_,nrow(obs),ncol(m))
 for(j in seq_len(ncol(m))) {
   invisible(capture.output(native <- fastwilcoxGMTall(setNames(m[,j],genes),gmts,outputGeneVals=FALSE,num.g=10)))
   ref[,j] <- native$db$stat[match(obs$pathway,rownames(native$db))]
 }
 target <- if(alt=='greater') ref>=obs$stat else abs(ref)>=abs(obs$stat)
 expected <- (1+rowSums(target,na.rm=TRUE))/(1+rowSums(!is.na(ref)))
 for(B in c(1,7,25)) {
   counts <- fcs_null_counts_by_block('global',obs,m,gmts,alt,num_g=10,block_cols=B)
   got <- fcs_close_counts(obs,counts)
   same(got$p.perm,expected)
 }
 cat('PASS native AUC oracle:',source_type,'blocks 1,7,N\n')
}
# Physical gene slices, including absent genes, overlap rejection and ID contracts.
for(G in c(1,7,400)) {
 dir<-tempfile('wp6_');dir.create(dir); batches<-caas_batches(genes,G)
 blocks<-split(seq_len(25),ceiling(seq_len(25)/7))
 for(b in names(blocks)) caas_transpose_block(raw[,blocks[[b]],drop=FALSE],batches,'global',b,dir)
 got<-unlist(lapply(seq_along(batches),function(k) caas_gene_p_batch(k,'global',names(blocks),dir,observed)))
 same(unname(got),unname(gp$nlog10_obs))
 for(b in names(blocks)) same(caas_assemble_block(seq_along(batches),'global',b,dir,genes),gp$nlog10_null[,blocks[[b]],drop=FALSE])
 reject(caas_gene_p_batch(1,'global',rep(names(blocks)[1],2),dir,observed))
 unlink(dir,recursive=TRUE)
}
reject(caas_transpose_block(matrix(1,2,2,dimnames=list(c('g','g'),c('b1','b2'))),list('g'),'global','1',tempdir()))
cat('PASS exhaustive gene p, physical slices at 1/7/G, overlaps and IDs\n')
