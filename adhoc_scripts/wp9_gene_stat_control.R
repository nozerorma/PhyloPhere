#!/usr/bin/env Rscript
# WP9 validation-only comparison. Reconstructs whole gene/cycle matrices deliberately; not a productive pipeline stage.
wp9_root <- function() {
 f<-grep('^--file=',commandArgs(FALSE),value=TRUE)
 if(length(f)) dirname(dirname(normalizePath(sub('^--file=','',f[1])))) else getwd()
}
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
wp9_ks <- function(x) {
 x<-sort(x[is.finite(x)&x>0]);n<-length(x)
 if(!n) return(NA_real_)
 max(seq_len(n)/n-x,x-(seq_len(n)-1)/n)
}
wp9_cor <- function(a,b) {
 ok<-is.finite(a)&is.finite(b)
 if(sum(ok)<3 || length(unique(a[ok]))<2 || length(unique(b[ok]))<2) return(NA_real_)
 cor(a[ok],b[ok],method='spearman')
}
wp9_mean <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
wp9_icc <- function(S,n) {
 k<-rowSums(n>0);ok<-k>=20
 if(sum(ok)<2) return(NA_real_)
 m<-rowSums(S)[ok]/k[ok]
 variance<-(rowSums(S*S)[ok]-k[ok]*m*m)/(k[ok]-1)
 within<-mean(variance);between<-max(var(m)-within/mean(k[ok]),0)
 if(between+within==0) NA_real_ else between/(between+within)
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
wp9_legacy <- function(positions,observed_file,null_file,universe,rankings,max_cells=5e7) {
 # The historical RDS supplies the full roster, including cycles with zero detections.
 reference_object<-readRDS(null_file)
 if(!is.null(reference_object$gene_stat)&&reference_object$gene_stat!='size_adj_max')stop('Legacy reconstruction requires the historical size_adj_max reference')
 reference<-reference_object$caas_corStat_byrank
 rm(reference_object)
 roster<-colnames(reference$global)
 if(length(roster)<2 || anyDuplicated(roster)) stop('Invalid historical roster')
 genes<-sort(unique(universe))
 if(length(genes)*length(roster)*length(rankings)>max_cells) stop('Whole-memory validation exceeds --max-cells')
 d<-data.table::as.data.table(readRDS(positions));o<-data.table::fread(observed_file,showProgress=FALSE)
 o<-o[is.finite(caas_score)&caas_score>0,.(Gene,Position,side,cycle='b_0',s=caas_score)]
 d<-rbind(d[,.(Gene,Position,side,cycle,s)],o)
 if(any(!d$cycle %in% c('b_0',roster))) stop('Position cycles outside historical roster')
 # Rebuild each historical cycle pool before selecting the fixed pathway universe.
 excluded<-setdiff(unique(d$Gene),genes)
 parts<-list();observed<-list();cy<-c('b_0',roster)
 for(rk in rankings) {
  a<-if(rk=='global') d else d[side==rk]
  a<-a[is.finite(s)&s>0,.(s=max(s)),by=.(cycle,Gene,Position)]
  g<-a[,.(n=.N,m=max(s)),by=.(cycle,Gene)]
  pools<-lapply(split(a$s,a$cycle),sort)
  g[,f:=mapply(function(v,c) findInterval(v+1e-12,pools[[c]])/length(pools[[c]]),m,cycle)]
  g<-g[Gene%in%genes]
  make<-function(v){m<-matrix(0,length(genes),length(cy),dimnames=list(genes,cy));m[cbind(match(g$Gene,genes),match(g$cycle,cy))]<-v;m}
  allparts<-list(max=make(g$m),fmax=make(g$f),n=make(g$n))
  parts[[rk]]<-lapply(allparts,function(m)m[,-1,drop=FALSE]);observed[[rk]]<-lapply(allparts,function(m)m[,1,drop=FALSE])
  ref<-reference[[rk]];reconstructed<-wp9_stat(parts[[rk]],'size_adj_max')
  common<-intersect(genes,rownames(ref))
  if(!identical(unname(reconstructed[common,]>0),unname(ref[common,roster]>0))) stop('Legacy detection differs from saved null')
  delta<-max(abs(reconstructed[common,]-ref[common,roster]))
  if(!is.finite(delta)||delta>1e-12)stop('Legacy reconstructed statistic differs from saved null')
  cat(sprintf('legacy %s size_adj_max max delta %.17g (no extra rounding)\n',rk,delta))
 }
 list(parts=parts,observed=observed,genes=genes,cycles=roster,source='historical_positions_no_extra_rounding',excluded_store_genes=excluded,inputs=c(positions,observed_file,null_file))
}
wp9_fcs <- function(S,obs,genes,gmts,num_g=10,max_g=0) {
 observed<-fcs_observed(list(global=setNames(obs,genes)),gmts,num_g,max_g,score_source='caas')$enrich
 if(is.null(observed)||!nrow(observed)) return(list(observed=data.frame(),cycles=data.frame(),stats=list()))
 real<-lapply(split(as.data.frame(observed),observed$database),function(d){rownames(d)<-d$pathway;d[,c('pval','stat'),drop=FALSE]})
 stats<-fcs_null_enrichstat_vectorized(S,gmts,real,num_g,max_g)
 observed$p.perm<-NA_real_;cycles<-list()
 for(db in names(stats)) {
  m<-stats[[db]];idx<-which(observed$database==db);o<-observed$stat[idx][match(rownames(m),observed$pathway[idx])]
  # For pathways, tolerance is exactly zero; do not import the tolerance of gene_p.
  p<-gene_p(m,o,tie_tol=0)
  observed$p.perm[idx]<-p$p_obs[match(observed$pathway[idx],rownames(m))]
  q<-apply(p$p_null,2,p.adjust,method='BH');if(is.null(dim(q))) q<-matrix(q,nrow=nrow(m))
  cycles[[db]]<-data.frame(database=db,cycle=colnames(S),valid=colSums(!is.na(p$p_null)),
   p01=colSums(p$p_null<=.01,na.rm=TRUE),p05=colSums(p$p_null<=.05,na.rm=TRUE),
   q15=colSums(q<.15,na.rm=TRUE))
 }
 observed<-fcs_classify_significance(observed,.15,.025)
 list(observed=observed,cycles=dplyr::bind_rows(cycles),stats=stats)
}
# Calibrate query columns only against the disjoint reference, with inclusive comparisons.
wp9_reference_p <- function(reference,query,tie_tol=0) {
 if(!identical(rownames(reference),rownames(query)))stop('Reference/query row IDs differ')
 p<-matrix(NA_real_,nrow(query),ncol(query),dimnames=dimnames(query))
 for(g in seq_len(nrow(query))) {
  ref<-sort(reference[g,!is.na(reference[g,])]);ok<-!is.na(query[g,])
  if(!length(ref))next
  p[g,ok]<-(1+length(ref)-findInterval(query[g,ok]-tie_tol,ref,left.open=TRUE))/(length(ref)+1)
 }
 p
}
wp9_holdout <- function(S,genes,gmts,real,num_g=10,max_g=0,seed=1998) {
 if(ncol(S)<4)stop('Holdout needs at least four cycles')
 set.seed(seed);permutation<-sample(seq_len(ncol(S)));folds<-split(permutation,rep(1:2,length.out=ncol(S)))
 rows<-list()
 for(k in 1:2) {
  evaluation<-folds[[k]];reference<-setdiff(seq_len(ncol(S)),evaluation)
  ref<-S[,reference,drop=FALSE]
  ref_gp<-gene_p(ref,rep(0,nrow(ref)))$nlog10_null
  eval_gp<--log10(wp9_reference_p(ref,S[,evaluation,drop=FALSE],1e-12))
  ref_auc<-fcs_null_enrichstat_vectorized(ref_gp,gmts,real,num_g,max_g)
  eval_auc<-fcs_null_enrichstat_vectorized(eval_gp,gmts,real,num_g,max_g)
  for(db in names(real)) {
   p<-wp9_reference_p(ref_auc[[db]],eval_auc[[db]],0)
   q<-apply(p,2,p.adjust,method='BH');if(is.null(dim(q)))q<-matrix(q,nrow=nrow(p))
   rows[[length(rows)+1]]<-data.frame(database=db,fold=k,cycle=colnames(S)[evaluation],reference_cycles=length(reference),
    valid=colSums(!is.na(p)),p01=colSums(p<=.01,na.rm=TRUE),p05=colSums(p<=.05,na.rm=TRUE),q15=colSums(q<.15,na.rm=TRUE))
  }
 }
 dplyr::bind_rows(rows)
}
wp9_run <- function(data,gmts,variants=wp9_variants,num_g=10,max_g=0,holdout=TRUE) {
 gene_rows<-list();gene_cycles<-list();path_rows<-list();path_cycles<-list();changes<-list();holdouts<-list()
 for(rk in names(data$parts)) {
  parts<-data$parts[[rk]];n<-parts$n;level<-colSums(n>0);obslevel<-sum(data$observed[[rk]]$n>0)
  near<-level>=.8*obslevel & level<=1.5*obslevel
  observed_ps<-list()
  for(v in variants) {
   t0<-Sys.time();S<-wp9_stat(parts,v);obs<-as.numeric(wp9_stat(data$observed[[rk]],v));gp<-gene_p(S,obs)
   observed_ps[[v]]<-gp$nlog10_obs
   g01<-colSums(gp$p_null<=.01,na.rm=TRUE);g05<-colSums(gp$p_null<=.05,na.rm=TRUE)
   prop<-rowMeans(n>0);mean_score<-rowSums(S)/pmax(1,rowSums(n>0));seen<-rowSums(n>0)>=.02*ncol(n)
   gene_q_null<-apply(gp$p_null,2,function(p)sum(p.adjust(p,'BH')<.15,na.rm=TRUE))
   gene_rows[[length(gene_rows)+1]]<-data.frame(ranking=rk,variant=v,genes=nrow(S),null_cycles=ncol(S),obs_scored=obslevel,matched_cycles=sum(near),
    ks=if(v=='max') NA_real_ else wp9_mean(apply(S,2,wp9_ks)),rho_n=wp9_mean(vapply(seq_len(ncol(S)),function(j){ok<-n[,j]>0;wp9_cor(S[ok,j],n[ok,j])},numeric(1))),
    rho_prop=wp9_cor(prop[seen],mean_score[seen]),icc=wp9_icc(S,n),obs_p01=sum(gp$p_obs<=.01),matched_p01=wp9_mean(g01[near]),share_p01=wp9_mean(g01[near]>=sum(gp$p_obs<=.01)),
    obs_p05=sum(gp$p_obs<=.05),matched_p05=wp9_mean(g05[near]),share_p05=wp9_mean(g05[near]>=sum(gp$p_obs<=.05)),rho_level=wp9_cor(level,g05),
    gene_q15=sum(p.adjust(gp$p_obs,'BH')<.15))
   gene_cycles[[length(gene_cycles)+1]]<-data.frame(ranking=rk,variant=v,cycle=colnames(S),level=level,observed_like=near,p01=g01,p05=g05,q15=gene_q_null)
   f<-wp9_fcs(gp$nlog10_null,gp$nlog10_obs,data$genes,gmts,num_g,max_g)
   if(holdout && nrow(f$observed)) {
    real<-lapply(split(as.data.frame(f$observed),f$observed$database),function(d){rownames(d)<-d$pathway;d[,c('pval','stat'),drop=FALSE]})
    h<-wp9_holdout(S,data$genes,gmts,real,num_g,max_g)
    h$ranking<-rk;h$variant<-v;h$level<-level[match(h$cycle,colnames(S))];h$observed_like<-near[match(h$cycle,colnames(S))]
    holdouts[[length(holdouts)+1]]<-h
   }
   if(nrow(f$observed)) {f$observed$variant<-v;f$observed$ranking<-rk;path_rows[[length(path_rows)+1]]<-f$observed}
   if(nrow(f$cycles)) {
    f$cycles$variant<-v;f$cycles$ranking<-rk;f$cycles$level<-level[match(f$cycles$cycle,colnames(S))]
    f$cycles$observed_like<-near[match(f$cycles$cycle,colnames(S))];path_cycles[[length(path_cycles)+1]]<-f$cycles
   }
   cat(sprintf('done %s %s: %d genes, %d cycles, %d pathways, %.1f s\n',rk,v,nrow(S),ncol(S),nrow(f$observed),as.numeric(difftime(Sys.time(),t0,units='secs'))))
  }
  if(all(c('fmax','size_adj_max')%in%names(observed_ps))) {
   ok<-as.numeric(data$observed[[rk]]$n)>0;delta<-abs(observed_ps$fmax-observed_ps$size_adj_max)
   changes[[rk]]<-data.frame(ranking=rk,obs_scored=sum(ok),n_gt8=sum(data$observed[[rk]]$n>8),
      fraction_delta_gt1=wp9_mean(delta[ok]>1),spearman=wp9_cor(observed_ps$fmax[ok],observed_ps$size_adj_max[ok]))
  }
 }
 list(gene_summary=dplyr::bind_rows(gene_rows),gene_cycles=dplyr::bind_rows(gene_cycles),pathways=dplyr::bind_rows(path_rows),pathway_cycles=dplyr::bind_rows(path_cycles),gene_changes=dplyr::bind_rows(changes),holdout_cycles=dplyr::bind_rows(holdouts))
}
wp9_main <- function(args=commandArgs(TRUE)) {
 if('--help'%in%args) {
  cat('WP9 validation-only, whole-memory control. Required: --out NEW_DIR --universe IDS\n',
      'Canonical: --store caas_permulation/\n',
      'Historical: --legacy-positions pos.rds --legacy-observed b0.tsv.gz --legacy-null null.rds\n',
      'Optional: --rankings global,top,bottom --variants max,fmax,beta025,beta05,beta075,size_adj_max\n',
      '--databases h,pathway_KEGG --min-genes 10 --max-genes 0 --max-cells 50000000\n',
      'Outputs: gene summaries/cycles, observed pathways, pathway LOO cycles, disjoint holdout cycles, metadata and sessionInfo.\n')
  return(invisible(NULL))
 }
 get<-function(k,default=NULL){i<-match(k,args);if(is.na(i)) default else {if(i==length(args))stop(paste('Missing value',k));args[i+1]}}
 root<-wp9_root();src<-file.path(root,'subworkflows/ENRICHMENT/local/src')
 suppressPackageStartupMessages({library(dplyr);library(RERconverge);library(data.table)})
 source(file.path(src,'percentile_flags.R'));source(file.path(src,'fcs_enrich.R'));source(file.path(src,'gene_p.R'))
 source(file.path(root,'subworkflows/SCORING/local/src/caas_block_contract.R'))
 output<-get('--out');if(is.null(output)||dir.exists(output))stop('--out must name a new directory')
 ufile<-get('--universe');universe<-trimws(readLines(ufile));universe<-universe[nzchar(universe)]
 if(anyDuplicated(universe))stop('Duplicate universe IDs')
 rankings<-strsplit(get('--rankings','global,top,bottom'),',',fixed=TRUE)[[1]]
 if(!length(rankings)||anyDuplicated(rankings)||any(!rankings%in%c('global','top','bottom')))stop('Invalid rankings')
 variants<-strsplit(get('--variants',paste(wp9_variants,collapse=',')),',',fixed=TRUE)[[1]]
 if(anyDuplicated(variants)||any(!variants%in%wp9_variants))stop('Invalid variants')
 gmtdir<-get('--gmt-dir',file.path(root,'subworkflows/ENRICHMENT/dat'));gmts<-fcs_load_gmts(gmtdir)
 dbs<-get('--databases');if(!is.null(dbs)){dbs<-strsplit(dbs,',',fixed=TRUE)[[1]];if(any(!dbs%in%names(gmts)))stop('Unknown database');gmts<-gmts[dbs]}
 max_cells<-as.numeric(get('--max-cells','50000000'));if(!is.finite(max_cells)||max_cells<1)stop('Invalid max-cells')
 store<-get('--store');data<-if(!is.null(store))wp9_store(store,universe,rankings,max_cells) else
   wp9_legacy(get('--legacy-positions'),get('--legacy-observed'),get('--legacy-null'),universe,rankings,max_cells)
 expected<-get('--expected-cycles')
 if(!is.null(expected)&&(!is.finite(as.numeric(expected))||as.numeric(expected)!=length(data$cycles)))stop('Cycle count differs from --expected-cycles')
 dir.create(output,recursive=TRUE)
 num_g<-as.integer(get('--min-genes','10'));max_g<-as.integer(get('--max-genes','0'))
 if(is.na(num_g)||num_g<1||is.na(max_g)||max_g<0)stop('Invalid gene-set bounds')
 metadata<-list(status='running',source=data$source,excluded_store_genes=data$excluded_store_genes,genes=length(data$genes),cycles=data$cycles,rankings=rankings,variants=variants,
  databases=names(gmts),min_genes=num_g,max_genes=max_g,gene_tie_tol=1e-12,pathway_tie_tol=0,
  holdout='two disjoint reference/evaluation folds, swapped, seed 1998; reference gene ranks LOO, evaluation ranks against full reference',
  comparators='gene LOO across all N cycles; pathway LOO across all N calibrated columns, inclusive upper AUC tail; diagnostic, not independent replicates',
  code_md5=as.list(tools::md5sum(c(file.path(root,'adhoc_scripts/wp9_gene_stat_control.R'),file.path(src,c('fcs_enrich.R','gene_p.R'))))),
  inputs=as.list(tools::md5sum(c(data$inputs,ufile,list.files(gmtdir,pattern='[.]gmt$',full.names=TRUE)))),args=args)
 saveRDS(metadata,file.path(output,'metadata.rds'));writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
 result<-wp9_run(data,gmts,variants,num_g,max_g)
 for(n in names(result))write.table(result[[n]],file.path(output,paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
 metadata$status<-'complete'
 saveRDS(metadata,file.path(output,'metadata.rds'))
}
if(sys.nframe()==0L)wp9_main()
