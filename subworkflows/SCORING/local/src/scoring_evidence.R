#!/usr/bin/env Rscript
# Shared gene calibration and activity tables. Memory bounded by gene batch x N.
a <- commandArgs(TRUE)
arg <- function(k, default=NULL) { i<-match(k,a); if(is.na(i)) default else a[i+1L] }
script <- sub('^--file=', '', grep('^--file=',commandArgs(),value=TRUE)[1])
src <- dirname(normalizePath(script))
source(file.path(src,'../../../ENRICHMENT/local/src/gene_p.R'))
source(file.path(src,'../../../ENRICHMENT/local/src/caas_blocks.R'))
read <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
write <- function(x,p) write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
if (!identical(arg('--gene-stat','max'), 'max')) stop('Only gene_stat=max is supported')
nom <- as.numeric(arg('--nominal','0.05')); fdr <- as.numeric(arg('--fdr','0.05'))
stopifnot(is.finite(nom), nom>0,nom<=1,is.finite(fdr),fdr>0,fdr<=1)
genes <- read(arg('--genes')); positions <- read(arg('--positions')); stats <- read(arg('--stats'))
u <- if(file.exists(arg('--universe'))) trimws(readLines(arg('--universe'))) else stats$gene; u<-u[nzchar(u)&!u%in%c('gene','Gene')]
if(anyDuplicated(u)) stop('Duplicate genes in calibration universe')
if(anyDuplicated(genes$Gene)||anyDuplicated(stats$gene))stop('Duplicate observed genes')
files <- sort(list.files(arg('--blocks'),pattern='\\.rds$',full.names=TRUE))
dir.create('calibration',showWarnings=FALSE);dir.create('gene_lists',showWarnings=FALSE);dir.create('position_lists',showWarnings=FALSE)
batch_size<-as.integer(arg('--batch','1000'));stopifnot(batch_size>=1)
ranks<-c('global','top','bottom'); obs_cols<-c(global='score_global',top='score_top',bottom='score_bottom')
cycles<-character();block_ids<-sprintf('blk_%06d',seq_along(files))
for(f in files) {
 b<-readRDS(f);cy<-colnames(b$stat$global)
 if(is.null(cy)||anyNA(cy)||anyDuplicated(cy)||any(cy%in%cycles))stop('Invalid/overlapping null cycles')
 if(!identical(b$gene_stat,arg('--gene-stat','max')))stop('Gene statistic mismatch')
 for(r in ranks)if(!identical(colnames(b$stat[[r]]),cy))stop('Ranking cycles differ')
 cycles<-c(cycles,cy)
}
N<-length(cycles);available<-N>=2L && length(u)>0
batches<-caas_batches(u,batch_size)
# Read each large block once; subsequent gene batches read only their small slices.
if(available) {
 dir.create('raw_slices',showWarnings=FALSE)
 for(i in seq_along(files)) {
  b<-readRDS(files[i])
  for(rk in ranks)caas_transpose_block(b$stat[[rk]],batches,rk,block_ids[i],'raw_slices')
 }
}
# Keep the complete tested gene universe, including genes without observed signal.
missing<-setdiff(u,genes$Gene)
if(length(missing)) {
 extra<-genes[rep(NA_integer_,length(missing)),,drop=FALSE];extra$Gene<-missing
 for(n in intersect(c('n_positions','n_positions_top','n_positions_bottom','gene_caas_score','gene_caas_score_top','gene_caas_score_bottom'),names(extra)))extra[[n]]<-0
 genes<-rbind(genes,extra)
}
missing<-setdiff(u,stats$gene)
if(length(missing)) {
 extra<-stats[rep(NA_integer_,length(missing)),,drop=FALSE];extra$gene<-missing
 for(n in intersect(c('score_global','score_top','score_bottom'),names(extra)))extra[[n]]<-0
 stats<-rbind(stats,extra)
}
all_evidence<-list()
for(rk in ranks) {
 obs<-setNames(rep(0,length(u)),u); ix<-match(u,stats$gene);v<-as.numeric(stats[[obs_cols[[rk]]]][ix]);obs[!is.na(v)]<-v[!is.na(v)]
 p<-setNames(rep(NA_real_,length(u)),u);valid<-setNames(rep(NA_integer_,length(u)),u)
 if(available)for(k in seq_along(batches)) {
  rows<-batches[[k]];S<-matrix(0,length(rows),N,dimnames=list(rows,cycles))
  block_cycles<-list()
  for(i in seq_along(files)) {
   b<-readRDS(caas_slice_path('raw_slices',rk,k,block_ids[i]));cy<-b$cols;S[,cy]<-b$values;block_cycles[[i]]<-cy
  }
  g<-gene_p(S,obs[rows]);p[rows]<-g$p_obs;valid[rows]<-rowSums(!is.na(S))
  write(data.frame(gene=rows,nlog10_obs=g$nlog10_obs),file.path('calibration',sprintf('obs_%s_%04d.tsv',rk,k)))
  for(i in seq_along(files)) {
   cy<-block_cycles[[i]];saveRDS(list(rows=rows,cols=cy,values=g$nlog10_null[,cy,drop=FALSE]),caas_slice_path('calibration',paste0(rk,'_gp'),k,block_ids[i]))
  }
 }
 q<-p.adjust(p,'BH');tab<-data.frame(gene=u,ranking=rep(rk,length(u)),statistic=unname(obs),gene_stat=rep(arg('--gene-stat','max'),length(u)),p.gen=unname(p),p.gen.adj=unname(q),nominal=!is.na(p)&p<nom,adjusted=!is.na(q)&q<fdr,N_roster=rep(N,length(u)),N_valid=unname(valid),family_id=rep(paste0('CAAS:',rk),length(u)),contract_version=rep(1L,length(u)))
 all_evidence[[rk]]<-tab
 suffix<-if(rk=='global')'' else paste0('_',rk)
 for(nm in c('p.gen','p.gen.adj','nominal','adjusted'))genes[[paste0('caas_',nm,suffix)]]<-tab[[nm]][match(genes$Gene,u)]
 for(nm in c('p.gen','p.gen.adj','nominal','adjusted'))stats[[paste0('caas_',nm,suffix)]]<-tab[[nm]][match(stats$gene,u)]
 # Header + tab-separated rows: prepare_ami_lists.py reads these columns by name (base write() cannot write a data.frame).
 active<-tab[tab$nominal,,drop=FALSE];write.table(data.frame(Gene=active$gene,active[,-1,drop=FALSE]),file.path('gene_lists',paste0('active_',rk,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
write(do.call(rbind,all_evidence),'gene_evidence.tsv')
# Existing position p-values and BH already use the pooled unique-position family.
positions$caas_nominal<-!is.na(positions$p.emp)&positions$p.emp<nom
positions$caas_adjusted<-!is.na(positions$p.adj_bh)&positions$p.adj_bh<fdr
for(rk in ranks) {
 x<-positions[positions$caas_nominal & (rk=='global' | positions$side==rk),,drop=FALSE]
 if(rk=='global')x<-x[!duplicated(paste(x$Gene,x$Position,sep='\t')),,drop=FALSE]
 write(x,file.path('position_lists',paste0('active_',rk,'.tsv')))
}
for(n in c("rer_min_pval","rer_perm_padj")) if(n %in% names(genes)) stats[[if(n=="rer_min_pval")"rer_p.emp" else "rer_p.adj"]] <- genes[[n]][match(stats$gene,genes$Gene)]
write(genes,'gene_scores.tsv');write(positions,'position_scores.tsv');write(stats,'fcs_stats.tsv')
jsonlite::write_json(list(contract_version=1L,available=available,N=N,genes=length(u),blocks=block_ids,nominal_threshold=nom,fdr_threshold=fdr,gene_stat=arg('--gene-stat','max')), 'calibration/metadata.json',auto_unbox=TRUE,pretty=TRUE)
file.copy('calibration/metadata.json','gene_evidence_metadata.json',overwrite=TRUE)
if(!available)write(data.frame(ranking=character(),database=character(),pathway=character(),stat=numeric(),pval=numeric(),p.adj=numeric(),p.perm=numeric(),q.perm=numeric(),num.genes=integer(),gene.vals=character(),sig_perm=logical(),direction=character(),redlich_level=character()),'calibration/empty.tsv')
