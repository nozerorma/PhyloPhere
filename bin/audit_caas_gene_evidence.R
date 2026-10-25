#!/usr/bin/env Rscript
# Validation only: native gene p and an explicitly experimental detection x
# pooled-conditional-tail candidate. Does not update production evidence/lists.
argval <- function(a,k,default=NULL) {i<-match(k,a);if(is.na(i))default else a[i+1L]}
classes_for <- function(b,nd) findInterval(pmax(nd,1),b,left.open=TRUE)+1L
fit_fact_gene <- function(S,K=20L) {
  nd<-rowSums(S>0,na.rm=TRUE);v<-sort(unique(nd[nd>0]))
  if(!length(v))b<-numeric() else {
    w<-v*tabulate(match(nd,v),nbins=length(v));cl<-pmin(K,ceiling(K*cumsum(w)/sum(w)-1e-9))
    last<-c(which(diff(cl)!=0),length(cl));b<-head(v[last],-1)
  }
  cls<-classes_for(b,nd)
  pools<-lapply(seq_len(length(b)+1L),function(k){x<-as.numeric(S[cls==k,,drop=FALSE]);sort(x[!is.na(x)&x>0])})
  list(nd=nd,breaks=b,pools=pools,N=ncol(S))
}
fact_gene_p <- function(value,nd,fit) {
  cls<-classes_for(fit$breaks,nd+1);p<-(nd+1)/(fit$N+1)
  for(k in seq_along(fit$pools)) {
    at<-which(cls==k & !is.na(value));pool<-fit$pools[[k]]
    if(length(at)&&length(pool))p[at]<-p[at]*(1+length(pool)-findInterval(value[at]-1e-12,pool,left.open=TRUE))/(1+length(pool))
  }
  p[value<=1e-12 & !is.na(value)]<-1;p[is.na(value)]<-NA_real_
  pmin(1,p)
}
write_tsv <- function(x,p) write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE)
main <- function(a) {
  repo<-argval(a,'--repo');null<-argval(a,'--null-dir');out<-argval(a,'--output-dir')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  source(file.path(repo,'subworkflows/ENRICHMENT/local/src/gene_p.R'))
  genes<-trimws(readLines(argval(a,'--universe')));genes<-genes[nzchar(genes)&!genes%in%c('Gene','gene')];stopifnot(!anyDuplicated(genes))
  obs<-read.delim(argval(a,'--core-genes'),check.names=FALSE)
  roster<-readLines(file.path(null,'cycle_roster.txt'));stopifnot(!anyDuplicated(roster),!'b_0'%in%roster)
  files<-sort(list.files(file.path(null,'blocks'),pattern='^blk_',full.names=TRUE));N<-length(roster);m<-length(genes)
  ranks<-c('global','top','bottom');columns<-c('gene_caas_score','gene_caas_score_top','gene_caas_score_bottom')
  S<-setNames(lapply(ranks,function(r)matrix(0,m,N,dimnames=list(genes,roster))),ranks)
  block_cycles<-list();seen<-character()
  for(i in seq_along(files)) {
    b<-readRDS(file.path(files[i],'matrices.rds'));stopifnot(identical(b$gene_stat,'size_adj_max'))
    cy<-b$cycles;stopifnot(!any(cy%in%seen),all(cy%in%roster));seen<-c(seen,cy);block_cycles[[i]]<-cy
    for(r in ranks){at<-match(genes,rownames(b$stat[[r]]));ok<-!is.na(at);S[[r]][ok,cy]<-b$stat[[r]][at[ok],,drop=FALSE]}
    cat('read',basename(files[i]),'\n');flush.console()
  }
  stopifnot(setequal(seen,roster));rm(b);gc()
  b0<-readRDS(file.path(null,'b0','matrices.rds'))
  allobs<-list();allcal<-list();allstrata<-list();b0comp<-list();settings<-list()
  set.seed(1998)
  for(ri in seq_along(ranks)) {
    rk<-ranks[ri];mat<-S[[rk]];o<-obs[[columns[ri]]][match(genes,obs$Gene)];o[is.na(o)]<-0
    obs0<-b0$stat[[rk]][match(genes,rownames(b0$stat[[rk]])),1];obs0[is.na(obs0)]<-0
    b0comp[[rk]]<-list(different=sum(abs(o-obs0)>1e-12),max_error=max(abs(o-obs0)))
    p<-valid<-numeric(m)
    batches<-split(seq_len(m),ceiling(seq_len(m)/500))
    for(ix in batches){g<-gene_p(mat[ix,,drop=FALSE],o[ix]);p[ix]<-g$p_obs;valid[ix]<-rowSums(!is.na(mat[ix,,drop=FALSE]))}
    stopifnot(all(valid==N))
    fit<-fit_fact_gene(mat);pf<-fact_gene_p(o,fit$nd,fit)
    tab<-data.frame(gene=genes,ranking=rk,statistic=o,N_valid=valid,p.gen=p,p.gen.adj=p.adjust(p,'BH'),experimental_p.gen_fact=pf,experimental_q.gen_fact=p.adjust(pf,'BH'),nd=fit$nd,class_observed=classes_for(fit$breaks,fit$nd+1))
    allobs[[rk]]<-tab
    settings[[rk]]<-list(class_breaks=fit$breaks,pool_sizes=lengths(fit$pools))
    rm(fit);gc()
    for(fold in 0:1) {
      train_blocks<-which((seq_along(files)-1)%%2==fold);test_blocks<-which((seq_along(files)-1)%%2!=fold)
      tr<-unlist(block_cycles[train_blocks],use.names=FALSE);te<-unlist(lapply(block_cycles[test_blocks],function(cy)sample(cy,min(10,length(cy)))),use.names=FALSE)
      train<-mat[,tr,drop=FALSE];test<-mat[,te,drop=FALSE];fit<-fit_fact_gene(train)
      pe<-pf<-matrix(1,m,length(te));nd<-fit$nd
      for(g in seq_len(m)) {
        v<-test[g,];ok<-which(v>1e-12);sx<-sort(train[g,])
        if(length(ok)) {
          pe[g,ok]<-(1+length(sx)-findInterval(v[ok]-1e-12,sx,left.open=TRUE))/(length(sx)+1)
          pf[g,ok]<-fact_gene_p(v[ok],rep(nd[g],length(ok)),fit)
        }
      }
      grp<-findInterval(nd,c(0,1,2,6,21,101,501,2501,5001),rightmost.closed=TRUE)
      for(method in c('p.gen','experimental_p.gen_fact')) {
        P<-if(method=='p.gen')pe else pf
        for(j in seq_along(te)) {
          q<-p.adjust(P[,j],'BH')
          for(alpha in c(.0001,.001,.005,.01,.025,.05))allcal[[length(allcal)+1L]]<-data.frame(ranking=rk,fold=fold,cycle=te[j],method=method,alpha=alpha,rejected=sum(P[,j]<=alpha),family=m,rate=mean(P[,j]<=alpha),bh_005=sum(q<.05))
        }
        for(k in sort(unique(grp)))for(alpha in c(.0001,.001,.005,.01,.025,.05)) {
          Q<-P[grp==k,,drop=FALSE];rates<-colMeans(Q<=alpha)
          allstrata[[length(allstrata)+1L]]<-data.frame(ranking=rk,fold=fold,group=k,method=method,alpha=alpha,genes=sum(grp==k),cycles=length(te),rate=mean(rates),sd_cycles=sd(rates),rejected=sum(Q<=alpha),pairs=length(Q))
        }
      }
      rm(train,test,fit,pe,pf);gc();cat(rk,'fold',fold,'done\n');flush.console()
    }
    S[[rk]]<-NULL;rm(mat);gc()
  }
  z<-do.call(rbind,allobs);write_tsv(z,file.path(out,'gene_evidence.tsv'));write_tsv(do.call(rbind,allcal),file.path(out,'gene_calibration_cycles.tsv'));write_tsv(do.call(rbind,allstrata),file.path(out,'gene_calibration_strata.tsv'))
  published<-argval(a,'--published-evidence')
  errors<-NULL
  if(!is.null(published)) {
    old<-read.delim(published,check.names=FALSE);j<-merge(z,old,by=c('gene','ranking'),suffixes=c('_audit','_published'));stopifnot(nrow(j)==nrow(z),nrow(old)==nrow(z))
    errors<-sapply(c('statistic','p.gen','p.gen.adj'),function(k)max(abs(j[[paste0(k,'_audit')]]-j[[paste0(k,'_published')]])))
    stopifnot(max(errors)<1e-12)
  }
  jsonlite::write_json(list(N=N,genes=m,family='separate BH by ranking over supplied cleaned universe',gene_stat='size_adj_max',nominal_threshold=.05,fdr_threshold=.05,b0_comparison=b0comp,published_max_error=errors,experimental=list(estimand='P(size_adj_max > 0) * P(size_adj_max >= observed | size_adj_max > 0)',class_fit=settings,production_selector_changed=FALSE)),file.path(out,'gene_summary.json'),pretty=TRUE,auto_unbox=TRUE)
}
if(sys.nframe()==0L)main(commandArgs(TRUE))
