suppressPackageStartupMessages(library(RERconverge))
source('subworkflows/RERCONVERGE/local/rer_batches.R')
same <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
res <- data.frame(Rho=c(-.7,0,.5,NA,.9), P=c(.01,1,.1,NA,0), row.names=paste0('g',1:5))
attr(res,'rer_settings') <- list(mode='binary',min.sp=15,min.pos=3,clade='terminal',winR=2)
make_batch <- function(n) {
 rho<-matrix(runif(5*n,-1,1),5,dimnames=list(rownames(res),NULL));rho[4,]<-NA
 p<-abs(rho);p[5,1]<-0
 list(corRho=as.data.frame(rho),corStat=as.data.frame(sign(rho)*-log10(p)))
}
parts <- lapply(1:3,function(i) {set.seed(1998+i-1);make_batch(7)})
full <- list(corRho=do.call(cbind,lapply(parts,`[[`,'corRho')))
expected <- permpvalcor(res,full)$permpval
for(G in c(1,2,5)) {
 dir<-tempfile('rerblocks_')
 got<-rer_write_batches(res,3,7,1998,make_batch,dir,gene_batch_size=G)
 same(got$p.perm,expected);same(got$p.perm.adj,p.adjust(expected,'BH'))
 stopifnot(attr(got,'n_perms')==21)
 for(i in 1:3) {
  b<-readRDS(file.path(dir,sprintf('blk_%06d.rds',i)))
  same(unname(b$corStat),unname(as.matrix(parts[[i]]$corStat)))
  stopifnot(is.null(b$corP),is.null(b$corRho),b$seed==1998+i-1,
            identical(b$rer_settings,attr(res,'rer_settings')))
 }
 stopifnot(!dir.exists(paste0(dir,'.gene_slices')))
 unlink(dir,recursive=TRUE)
}
cat('PASS RER seeded batches reproduce native gene p; gene batches 1/2/G; unbounded corStat; settings and cleanup\n')
