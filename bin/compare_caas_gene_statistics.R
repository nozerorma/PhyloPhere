#!/usr/bin/env Rscript
# Read-only scientific comparison: existing canonical nulls, no simulation or selection changes.
# Run on Slurm: Rscript bin/compare_caas_gene_statistics.R --repo REPO --store CAAS_STORE --universe CLEAN_IDS --out NEW_DIR
# Requires canonical max/fmax/n components and 10000 cycles. No FCS computation.
arg <- function(k) commandArgs(TRUE)[match(k,commandArgs(TRUE))+1L]
repo<-arg('--repo');store<-arg('--store');out<-arg('--out');uf<-arg('--universe')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(jsonlite))
wp9_variants <- c('max','fmax','beta025','beta05','beta075','size_adj_max')
wp9_stat <- function(parts, variant) {
 if(!variant %in% wp9_variants) stop('Unknown WP9 variant')
 if(!all(c('max','fmax','n') %in% names(parts))) stop('Missing WP9 parts')
 if(!identical(dimnames(parts$max),dimnames(parts$fmax)) || !identical(dimnames(parts$n),dimnames(parts$fmax))) stop('WP9 parts IDs differ')
 n<-parts$n; f<-parts$fmax
 if(any(!is.finite(n)|n<0|n!=floor(n)) || any(!is.finite(f)|f<0|f>1) || any(!is.finite(parts$max)|parts$max<0)) stop('Invalid WP9 parts')
 if(any((n==0)!=(f==0)) || any((n==0)!=(parts$max==0))) stop('WP9 parts detection differs')
 if(variant=='max') return(parts$max)
 beta<-c(fmax=0,beta025=.25,beta05=.5,beta075=.75,size_adj_max=1)[[variant]]
 # Preserve exact canonical endpoint operations and keep absent genes at zero, including n^0.
 x<-f
 ok<-n>0
 x[ok]<-f[ok]^(n[ok]^beta)
 x
}

wp9_cor <- function(a,b) {
 ok<-is.finite(a)&is.finite(b)
 if(sum(ok)<3 || length(unique(a[ok]))<2 || length(unique(b[ok]))<2) return(NA_real_)
 cor(a[ok],b[ok],method='spearman')
}

wp9_store <- function(store,universe,rankings,max_cells=5e7) {
 files<-Sys.glob(file.path(store,'blocks/blk_*/matrices.rds'))
 roster<-readLines(file.path(store,'cycle_roster.txt'))
 index<-caas_block_index(files,roster)
 meta_path<-file.path(store,'null_metadata.json')
 if(file.exists(meta_path)) {
  saved<-jsonlite::read_json(meta_path,simplifyVector=FALSE)
  if(!identical(saved$gene_stat,index$gene_stat)||!identical(as.character(unlist(saved$cycles)),index$cycles)||
     !identical(as.character(unlist(saved$genes)),index$genes))stop('WP9 store metadata disagrees with blocks')
  hashes<-vapply(index$blocks,function(b)b$md5,'');old<-vapply(saved$blocks,function(b)b$md5,'')
  if(!identical(sort(hashes),sort(old)))stop('WP9 store content differs from metadata')
 }
 genes<-sort(unique(universe))
 # The productive CAAS transposition selects the fixed supplied universe from the broader block store.
 excluded<-setdiff(index$genes,genes)
 if(length(genes)*length(roster)*length(rankings)>max_cells) stop('Whole-memory validation exceeds --max-cells; use one ranking per job or explicitly budget more memory')
 make<-function(cy) matrix(0,length(genes),length(cy),dimnames=list(genes,cy))
 parts<-setNames(lapply(rankings,function(rk)setNames(lapply(c('max','fmax','n'),function(p)make(roster)),c('max','fmax','n'))),rankings)
 for(f in files) {
  b<-readRDS(f)
  if(is.null(b$parts)) stop('Canonical blocks lack WP9 parts; do not infer them from stat')
  common<-intersect(b$genes,genes)
  for(rk in rankings)for(p in names(parts[[rk]]))parts[[rk]][[p]][common,b$cycles]<-b$parts[[p]][[rk]][common,,drop=FALSE]
 }
 b0path<-file.path(store,'b0/matrices.rds');b0<-readRDS(b0path)
 if(is.null(b0$parts)||!identical(b0$cycles,'b_0')||anyDuplicated(b0$genes)) stop('Invalid observed b0 parts/IDs')
 if(!identical(b0$gene_stat,index$gene_stat))stop('Observed b0 gene_stat differs from null')
 for(rk in rankings) {
  bp<-lapply(b0$parts,function(p)p[[rk]])
  expected<-wp9_stat(bp,b0$gene_stat)
  if(!isTRUE(all.equal(b0$stat[[rk]],expected,tolerance=1e-12)))stop('Observed b0 stat disagrees with parts')
 }
 observed<-setNames(lapply(rankings,function(rk)setNames(lapply(c('max','fmax','n'),function(p){
  m<-make('b_0');v<-b0$parts[[p]][[rk]]
  if(!identical(rownames(v),b0$genes)||!identical(colnames(v),'b_0')) stop('Observed b0 parts dimensions/IDs differ')
  common<-intersect(b0$genes,genes);m[common,'b_0']<-v[common,'b_0'];m
 }),c('max','fmax','n'))),rankings)
 list(parts=parts,observed=observed,genes=genes,cycles=roster,source='canonical_blocks',excluded_store_genes=excluded,inputs=c(files,b0path,file.path(store,'cycle_roster.txt'),if(file.exists(meta_path))meta_path))
}

source(file.path(repo,'subworkflows/SCORING/local/src/caas_block_contract.R'))
genes<-sort(unique(trimws(readLines(uf))));stopifnot(length(genes)>1000,!any(!nzchar(genes)))
write_tsv<-function(d,f)write.table(d,file.path(out,f),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
summary<-pairs<-cal<-levels<-list();meta<-list()
for(rk in c('global','top','bottom')) {
 cat('loading',rk,'\n');flush.console()
 d<-wp9_store(store,genes,rk,max_cells=2e8);P<-d$parts[[rk]];O<-d$observed[[rk]]
 N<-length(d$cycles);m<-length(genes);stopifnot(N==10000)
 # Verify reconstructed endpoint against every canonical saved block.
 files<-Sys.glob(file.path(store,'blocks/blk_*/matrices.rds'))
 blockcy<-lapply(files,function(f)readRDS(f)$cycles)
 for(f in files){b<-readRDS(f);ix<-match(genes,b$genes);ok<-!is.na(ix);S0<-wp9_stat(lapply(P,function(x)x[ok,b$cycles,drop=FALSE]),'size_adj_max');stopifnot(max(abs(S0-b$stat[[rk]][ix[ok],,drop=FALSE]))<1e-12)}
 rm(b,S0);gc()
 set.seed(1998);design<-lapply(0:1,function(fold){list(ref=unlist(blockcy[which((seq_along(files)-1)%%2==fold)],use.names=FALSE),test=unlist(lapply(blockcy[which((seq_along(files)-1)%%2!=fold)],function(cy)sample(cy,10)),use.names=FALSE))})
 level<-colSums(P$n>0);obslevel<-sum(O$n>0);near<-level>=.8*obslevel & level<=1.5*obslevel
 levels[[rk]]<-data.frame(ranking=rk,cycle=d$cycles,detected_genes=level,observed_detected=obslevel,observed_like=near)
 observed<-list()
 for(v in c('max','fmax','size_adj_max')) {
  S<-wp9_stat(P,v);o<-as.numeric(wp9_stat(O,v));p<-rep(NA_real_,m)
  for(g in seq_len(m)){sx<-sort(S[g,]);stopifnot(length(sx)==N);p[g]<-(1+N-findInterval(o[g]-1e-12,sx,left.open=TRUE))/(N+1)}
  q<-p.adjust(p,'BH');observed[[v]]<-p
  tab<-data.frame(gene=genes,ranking=rk,variant=v,statistic=o,n_positions=as.numeric(O$n),p.gen=p,p.gen.adj=q)
  write_tsv(tab,paste0(rk,'_',v,'_genes.tsv'))
  null_rhon<-sapply(seq(1,N,by=50),function(j){ok<-P$n[,j]>0;wp9_cor(S[ok,j],P$n[ok,j])})
  summary[[length(summary)+1]]<-data.frame(ranking=rk,variant=v,genes=m,cycles=N,observed_detected=obslevel,observed_like_cycles=sum(near),nominal=sum(p<.05),bh=sum(q<.05),min_p=min(p),min_q=min(q),floor_genes=sum(p==1/(N+1)),rho_score_n=wp9_cor(o[o>0],as.numeric(O$n)[o>0]),null_mean_rho_n=mean(null_rhon,na.rm=TRUE))
  for(fi in 1:2){de<-design[[fi]];te<-match(de$test,d$cycles);tr<-match(de$ref,d$cycles);pp<-matrix(NA_real_,m,length(te))
   for(g in seq_len(m)){sx<-sort(S[g,tr]);pp[g,]<-(1+length(sx)-findInterval(S[g,te]-1e-12,sx,left.open=TRUE))/(length(sx)+1)}
   for(j in seq_along(te))cal[[length(cal)+1]]<-data.frame(ranking=rk,variant=v,fold=fi-1,cycle=d$cycles[te[j]],detected_genes=level[te[j]],observed_like=near[te[j]],reference_cycles=length(tr),family=m,p001=sum(pp[,j]<=.001),p01=sum(pp[,j]<=.01),p05=sum(pp[,j]<=.05),bh05=sum(p.adjust(pp[,j],'BH')<.05))
   rm(pp)
  }
  cat('done',rk,v,'nominal',sum(p<.05),'BH',sum(q<.05),'qmin',min(q),'\n');flush.console();rm(S);gc()
 }
 for(vv in list(c('max','fmax'),c('fmax','size_adj_max'),c('max','size_adj_max'))){a<-observed[[vv[1]]];b<-observed[[vv[2]]];ok<-as.numeric(O$n)>0
  pairs[[length(pairs)+1]]<-data.frame(ranking=rk,a=vv[1],b=vv[2],spearman_all=wp9_cor(a,b),spearman_detected=wp9_cor(a[ok],b[ok]),common_nominal=sum(a<.05&b<.05),a_only=sum(a<.05&b>=.05),b_only=sum(b<.05&a>=.05),delta_log10_gt1=sum(abs(log10(a)-log10(b))>1),delta_log10_gt05=sum(abs(log10(a)-log10(b))>.5))
 }
 meta[[rk]]<-list(inputs=as.list(tools::md5sum(c(d$inputs,uf))),excluded_store_genes=d$excluded_store_genes,design=design)
 rm(d,P,O,observed);gc()
}
write_tsv(do.call(rbind,summary),'summary.tsv');write_tsv(do.call(rbind,pairs),'pairwise.tsv');write_tsv(do.call(rbind,cal),'calibration_cycles.tsv');write_tsv(do.call(rbind,levels),'detection_levels.tsv')
write_json(list(status='complete',method='inclusive own-gene upper tail; tolerance 1e-12; add-one; 10000 cycles; BH separately by ranking/variant; 200 disjoint evaluation cycles, 5000 reference per fold',seed=1998,input_metadata=meta),file.path(out,'metadata.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
